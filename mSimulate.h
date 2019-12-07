#ifndef M_SIMULATE_H
#define M_SIMULATE_H

#include "mDefinitions.h"
#include "mCommon.h"
#include "mMatrix.h"
#include "dc.h"

/* distances */

#define EUCLIDEAN  (1)
#define JSD        (2)
#define KLD        (3)
#define NONE      (-1)

#define NORMAL   (9)
#define MANI     (10)
#define UNIFORM  (11)
#define TEMPLATE (12)
#define TEMPLATE_NORMAL (13)
#define TEMPLATE_UNIF   (14)

struct mSample {
	int size;
	int groups;
	mIVector *pop;
	mFVector *abundance;
};
typedef struct mSample mSample;

void mInitSample(mSample* p, int size, int groups);
void mFreeSample(mSample* p);
void mWriteSample(FILE* stream, mSample* s);

struct mSimulation {
	int       nsamples;
	int       nfeatures;
	int       nfragments;
	num_t     pseudocount;
	int       mode;
	int       unassigned;
	mSample **sample;
	mMatrix  *template;
};
typedef struct mSimulation mSimulation;

void mInitSimulation(mSimulation* s, int nsamples, int fragments);
void mWriteSimulation(FILE* stream, mSimulation* sim);
void mFreeSimulation(mSimulation* s);

/* general */

void mCalculateDistance(int skip, mSimulation* sim, mMatrix *dist);
void mSimulateSample(int mode, mSample* p, num_t assignment_rate);
void mShuffleFisherYates(int n, int* array);
void mCalculateAbundance(mSimulation* sim);
void mNormalizeAbundance(mSimulation* sim);
void mCalculatePopulation(mSimulation *sim);
void mAddPseudoCount(mSimulation *sim);

int simulate(mSimulation* simulation, mMatrix* dist, int samples, int features, int fragments);
void init_experiment(mSimulation* simulation, mMatrix* dist, int samples, int features, int fragments);

/* random number generation */

void mInitRand(int wordsize, int exponent, uint32_t seed, int generators);
void mFreeRand(int count);
void mSetGenerator(int gen);
num_t mGenerateUniformRandomFraction();
num_t mGenerateNormalRandomValue(num_t mu, num_t sigma);
num_t mGenerateUniformRandomValueInInterval(num_t a, num_t b);
uint32_t mGenerateRandomNumber();

#endif
