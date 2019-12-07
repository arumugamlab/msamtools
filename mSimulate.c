#include <limits.h>
#include <math.h>
#include "mSimulate.h"
#define VAR_FEATURES (0.00)
#define VAR_FRAGMENTS (0.10)

static const num_t mPI    = 3.14159265359; 
static const num_t mTwoPI = 6.28318530718;
static num_t mIntMaxPlusOneInverse   =   1.0/4294967296.0; /* this is UINT_MAX + 1 */
static num_t mIntMaxInverse = 1.0/4294967295.0; /* UINT_MAX */
mt_struct        **mtss;
static mt_struct  *rng1;
static mt_struct  *rng2;

/* For a request of N generators, we will actually have 2N RNGs. This is so that there *
 * are two uniform RNGs per instance so we can use two per instance to generate normal *
 * random variables as well using the same setup.                                      *
 * Therefore, generator itself is labelled as 2*self, and the second one is            *
 * generator+1 for the normal variate.
 */

void mInitRand(int wordsize, int exponent, uint32_t seed, int generators) {
	int i;
	int gen_count;
	mtss = get_mt_parameters_st(wordsize, exponent, 0, 2*generators-1, seed, &gen_count);
	if (gen_count != 2*generators)
		mDie("RNG init failed: %d out of %d generated!\n", gen_count, 2*generators);
	for (i=0; i<gen_count; i++)
		sgenrand_mt(i, mtss[i]);
}

void mFreeRand(int count) {
	free_mt_struct_array(mtss, 2*count);
}

void mSetGenerator(int gen) {
	rng1 = mtss[2*gen];
	rng2 = mtss[2*gen+1];
}

uint32_t mGenerateRandomNumber() { return genrand_mt(rng1); }

/* scaled by (1/(1+UINT_MAX)) so that the returned value is [0,1). */

num_t mGenerateUniformRandomFraction() { return (num_t) genrand_mt(rng1)*mIntMaxPlusOneInverse; }

num_t mGenerateUniformRandomValueInInterval(num_t a, num_t b) {
	return a+(b-a)*genrand_mt(rng1)*mIntMaxInverse;
}


/* since Box-Muller transform requires uniform random number in (0,1], add 1 and scale */

num_t mGenerateNormalRandomValue(num_t mu, num_t sigma) {
	num_t U1 = (1.0+genrand_mt(rng1))*mIntMaxPlusOneInverse;
	num_t U2 = (1.0+genrand_mt(rng2))*mIntMaxPlusOneInverse;
	return mu+sigma*cos(mTwoPI*U2)*sqrt(-2*log(U1));
}

/*
float mGenerateUniformRandomFraction() { return rand() / mRandMax; }
*/


void mInitSample(mSample* p, int size, int groups) {
	int i;
	mIVector *pop;
	mFVector *abundance;
	p->size = size;
	p->groups = groups;
	pop = (mIVector*) mMalloc(sizeof(mIVector));
	mInitIVector(pop, groups);
	for (i=0; i<groups; i++)
		mPushIVector(pop, 0);
	p->pop = pop;
	abundance = (mFVector*) mMalloc(sizeof(mFVector));
	mInitFVector(abundance, groups);
	for (i=0; i<groups; i++)
		mPushFVector(abundance, 0.0);
	p->abundance = abundance;
}

void mFreeSample(mSample* p) {
	mFreeFVector(p->abundance);
	mFree(p->abundance);
	mFreeIVector(p->pop);
	mFree(p->pop);
}

void mWriteSample(FILE* stream, mSample* s) {
	fprintf(stream, "<Sample size=\"%d\" features=\"%d\">\n", s->size, s->groups);
	fprintf(stream, "<Counts>\n");
	mWriteIVector(stream, s->pop, ",");
	fprintf(stream, "</Counts>\n");
	fprintf(stream, "<Abundance>\n");
	mWriteFVector(stream, s->abundance, ",");
	fprintf(stream, "</Abundance>\n");
	fprintf(stream, "</Sample>\n");
}

void mInitSimulation(mSimulation* s, int nsamples, int nfragments) {
	int i;
	s->nsamples = nsamples;
	s->nfragments = nfragments;
	s->sample   = (mSample**) mMalloc(nsamples*sizeof(mSample*));
	for (i=0; i<nsamples; i++)
		s->sample[i] = (mSample*) mMalloc(sizeof(mSample));
}

void mFreeSimulation(mSimulation* s) {
	int i;
	for (i=0; i<s->nsamples; i++) {
		mFreeSample(s->sample[i]);
		mFree(s->sample[i]);
	}
	mFree(s->sample);
}

void mWriteSimulation(FILE* stream, mSimulation* s) {
	int i;
	fprintf(stream, "<Simulation samples=\"%d\" features=\"%d\" fragments=\"%d\" pseudocount=\"%f\">\n", s->nsamples, s->nfeatures, s->nfragments, s->pseudocount);
	for (i=0; i<s->nsamples; i++) {
		fprintf(stream, "<Id index=\"%d\">\n", i);
		mWriteSample(stream, s->sample[i]);
		fprintf(stream, "</Id>\n");
	}
	fprintf(stream, "</Simulation>\n");
}

void init_experiment(mSimulation* simulation, mMatrix* dist, int samples, int features, int fragments) {

	int i;

	mInitSimulation(simulation, samples, fragments);
	mInitMatrix(dist, samples, samples);
	simulation->pseudocount = 0.000001;
	simulation->mode        = UNIFORM;

	/* Create an array of feature vectors, and fill them with 0 */

	simulation->nfeatures = features + (int) (VAR_FEATURES*features);    /* allocate max possible */
	for (i=0; i<samples; i++) {
		mInitSample(simulation->sample[i], fragments, simulation->nfeatures);
	}
}

int simulate(mSimulation* simulation, mMatrix* dist, int samples, int features, int fragments) {
	int   i;
	num_t sim_assigned = 1.0;                    /* assignment rate for this sample, default 100% */
	int   skip = (simulation->unassigned != 0);  /* should we skip any feature in distance calculation? like unassigned? */

	for (i=0; i<samples; i++) {
		if (simulation->unassigned) sim_assigned  = mGenerateUniformRandomFraction(); /* will never hit if unassigned == 0 */
		mSimulateSample(simulation->mode, simulation->sample[i], sim_assigned);
	}
/*	mWriteSimulation(stdout, simulation); */
	mNormalizeAbundance(simulation);
	mCalculatePopulation(simulation);
	mAddPseudoCount(simulation);
	mCalculateDistance(skip, simulation, dist);

	return 1;
}

void mSimulateSample(int mode, mSample* p, num_t assignment_rate) {
	int i;
	int sim_feature_dev;
	int sim_fragment_dev;
	int features  = p->groups;
	int fragments = p->size;
	int actual_fragments = 0;  /* reads in single-end, or inserts in paired-end */
	int actual_features  = 0;       /* number of features simulated based on original features */
	int range_max;
	int *pop;
	num_t *abundance;

	/* reset counts */

	mFillIVector(p->pop, 0);

	/* Recommended number of features is <features>, but simulator can randomly choose
	   in the interval [features - sim_feature_dev , features + sim_feature_dev]
	   Same for fragments as well.
	*/

	sim_feature_dev  = (int) (VAR_FEATURES*features);
	sim_fragment_dev = (int) (VAR_FRAGMENTS*fragments);

	actual_features  = features  - sim_feature_dev  + 2 * sim_feature_dev  * mGenerateUniformRandomFraction();
	actual_fragments = fragments - sim_fragment_dev + 2 * sim_fragment_dev * mGenerateUniformRandomFraction();

	switch (mode) {
		case MANI:

			/* you actually generate fragments, so need exact number of fragments */

			/* number of groups to populate, 1 reserved for unassigned */
			/* since range*frac() can be 0<=x<range, we use actual_features-1 and floor */
			range_max = actual_features - 1;
			pop    = p->pop->elem+1;       /* offset 1 to reduce number of additions */

			for (i=0; i<actual_fragments; i++) {
				if (mGenerateUniformRandomFraction() <= assignment_rate) {
					int f = (int) floor(range_max * mGenerateUniformRandomFraction());
					pop[f]++;
				} else {
					pop[-1]++;  /* because this is offset by 1 already */
				}
			}
			p->size   = actual_fragments;
			p->groups = features + sim_feature_dev;
			break;

		case NORMAL:
			abundance = p->abundance->elem;
			for (i=0; i<actual_features; i++) {
				if ((abundance[i] = mGenerateNormalRandomValue(0.2, 0.3)) < 0.0) abundance[i] = 0.0;
			}
			p->groups = features + sim_feature_dev;
			break;

		case UNIFORM:
			abundance = p->abundance->elem;
			for (i=0; i<actual_features; i++) {
				abundance[i] = mGenerateUniformRandomFraction();
			}
			p->groups = features + sim_feature_dev;
			break;

		case TEMPLATE:
			mDie("Template-based simulation shouldnt be here!\n");
		default:
			break;
	}

	/* shuffling not required if there are no missing values */
	/* and i am not sure if shuffling is a good idea anyway! */
	/* -1 for unassigned */
	/* mShuffleFisherYates(p->groups - 1, pop); */
}

void mSimulateSample_old(mSample* p, num_t assignment_rate) {
	int i;
	int size = p->size;
	int groups = p->groups - 1; /* 1 reserved for unassigned */
	for (i=0; i<size; i++) {
		if (mGenerateUniformRandomFraction() <= assignment_rate) {
			int f = (int) (groups * mGenerateUniformRandomFraction());
			p->pop->elem[1+f]++;
		} else {
			p->pop->elem[0]++;
		}
	}
}

void mShuffleFisherYates(int n, int* array) {
	int i;
	for (i=n; i>1; i--) {
		int j = (int) (i * mGenerateUniformRandomFraction());
		int tmp = array[j];
		array[j] = array[i-1];
		array[i-1] = tmp;
	}
}

void mNormalizeAbundance(mSimulation* sim) {
	int   i, j;
	int   nsamples    = sim->nsamples;
	int   nfeatures   = sim->nfeatures;
	for (i=0; i<nsamples; i++) {
		num_t *abundance = sim->sample[i]->abundance->elem;
		num_t sum = 0;
		for (j=0; j<nfeatures; j++) {
			sum += abundance[j];
		}
		for (j=0; j<nfeatures; j++)
			abundance[j] /= sum;
	}
}

void mAddPseudoCount(mSimulation *sim) {
	int i, j;
	int   nsamples    = sim->nsamples;
	int   nfeatures   = sim->nfeatures;
	num_t pseudocount = sim->pseudocount;
	num_t denominator = 1 + sim->nfeatures*pseudocount;
	for (i=0; i<nsamples; i++) {
		num_t   *abundance = sim->sample[i]->abundance->elem;
		for (j=0; j<nfeatures; j++) {
			abundance[j] += pseudocount;
			abundance[j] /= denominator;
		}
	}
}

void mCalculateAbundance(mSimulation* sim) {
	int i, j;
	int   nsamples    = sim->nsamples;
	num_t pseudocount = sim->pseudocount;
	num_t denominator = 1 + sim->nfeatures*pseudocount;
	for (i=0; i<nsamples; i++) {
		mSample *sample = sim->sample[i];
		int size = sample->size;
		mIVector *pop = sample->pop;
		mFVector *abundance = sample->abundance;
		for (j=0; j<sim->nfeatures; j++) {
			abundance->elem[j] = (1.0*pop->elem[j]/size + pseudocount)/denominator;
if (abundance->elem[j] <= 0) {
	fprintf(stdout, "sample %d has 0 from elem=%d, size=%d\n", i, pop->elem[j], size);
}
		}
	}
}

/* update population size */

void mCalculatePopulation(mSimulation *sim) {
	int samples = sim->nsamples;
	int size    = sim->nfragments;
	int features= sim->nfeatures;
	int k, j;
	for (k=0; k<samples; k++) {
		num_t sum = 0;
		mSample *s         = sim->sample[k];
		num_t   *abundance = s->abundance->elem;
		int     *pop       = s->pop->elem;
		for (j=0; j<features; j++) {
			pop[j] = size*abundance[j];
			sum += pop[j];
		}
		s->size = (int) sum;
	}

}

void mCalculateDistance(int skip, mSimulation* sim, mMatrix *dist) {
	int i, j;
	int samples = sim->nsamples;
	for (i=0; i<samples; i++) {
		dist->elem[i][i] = 0.0;
		for (j=i+1; j<samples; j++) {
			/* offset by 1, since unassigned should not be used in distance calculation */
			dist->elem[j][i] = mFJensenShannonDistance(sim->sample[i]->abundance, sim->sample[j]->abundance, skip); 
			if (isnan(dist->elem[j][i])) {
				dist->elem[j][i] = mFJensenShannonDistance(sim->sample[i]->abundance, sim->sample[j]->abundance, skip); 
				mWriteIVector(stdout, sim->sample[i]->pop, ",");
				mWriteFVector(stdout, sim->sample[i]->abundance, ",");
				mWriteIVector(stdout, sim->sample[j]->pop, ",");
				mWriteFVector(stdout, sim->sample[j]->abundance, ",");
				mDie("JSD between samples %d and %d is NaN\n", i, j);
			}
			dist->elem[i][j] = dist->elem[j][i];
		}
	}
}
