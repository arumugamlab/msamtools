#include "msam.h"
#include "zoeTools.h"

#define MULTI_ADD_ALL (1)
#define MULTI_SHARE_EQUAL (2)
#define MULTI_SHARE_PROPORTIONAL (3)

#define UNIT_REL (1)
#define UNIT_FPKM (2)
#define UNIT_TPM (3)
#define UNIT_ABN (4)

void mInitInsertCounts(int share_type);
int  mEstimateInsertCountOnFile(samfile_t *input, int share_type);
void mWriteCompressedSeqAbundance();
mMatrix* mInsertCountToAbundanceMatrix(int row, const char* rowname, int share_type);

void mInitInsertCounts(int share_type) {
	int n_features  = global->n_features;

	global->ui_insert_count = (uint32_t*) mCalloc(n_features, sizeof(uint32_t));
	global->ub_target_hit   = (uint8_t*) mCalloc(n_features, sizeof(uint8_t));

	if (share_type == MULTI_SHARE_EQUAL) {
		int i;
		double *f = (double*) mMalloc(n_features*sizeof(double));
		for (i=0; i<n_features; i++) f[i] = 0.0f;
		global->d_insert_count = f;
	}

	if (share_type == MULTI_SHARE_PROPORTIONAL) {
		global->multi_mappers = (mVector*) mMalloc(sizeof(mVector));
		mInitVector(global->multi_mappers, 65536);
	}
}

void mFreeInsertCounts(int share_type) {
	if (share_type == MULTI_SHARE_PROPORTIONAL) {
		int i;
		mVector *mm = global->multi_mappers;
		for (i=0; i<mm->size; i++) {
			mIVector *vec = (mIVector*) mm->elem[i];
			mFreeIVector(vec);
			mFree(vec);
		}
		mFreeVector(mm);
		mFree(mm);
	}

	if (share_type == MULTI_SHARE_EQUAL) {
		mFree(global->d_insert_count);
	}

	mFree(global->ui_insert_count);
	mFree(global->ub_target_hit);
}

/* This mode makes sure that the total insert counts are preserved in the end */
/* It shares a total count of 1 among the seen sequences */
void mEstimateInsertCountOnPool(mBamPool *pool, int share_type) {
	bam1_t **elem  = pool->elem;
	int     *fmap  = global->fmap;
	int      size  = pool->size;
	int      i;

/* Note the multiple return's in this function for avoiding too many conditionals.
 * Maintaining might be difficult, but I chose this over if-else-if-else-etc */

	switch(size) {
		case 1: /* With one hit, you got to be unique hit, so add 2 no matter what share_type */
			global->ui_insert_count[fmap[elem[0]->core.tid]] += 2;
			global->uniq_mapper_count++;
			return;

		case 2: /* could be 2 mate-pairs, or 1 mate mapped twice*/
		{
			/* Doing this separately so that unique-mapping count can be correctly estimated */
			/* Otherwise it is similar to what happens in "default:" */
			int fid0  = fmap[elem[0]->core.tid];
			int fid1  = fmap[elem[1]->core.tid];

			if (fid0 == fid1) { /* Hits twice the same target, so add 2 independent of share_type and mate status */
				global->ui_insert_count[fid0] += 2;
				global->uniq_mapper_count++;
				return;
			}

			/* If we are here, fid0 and fid1 are different, so multi-mappers */
	/* TO DO: Should I differentiate between mate-types? Probably not! */
			global->multi_mapper_count++;
			switch(share_type) {
				case MULTI_ADD_ALL:     /* For 2 distinct hits, MULTI_ADD_ALL adds 2 to each independent of mate status */
					global->ui_insert_count[fid0]+=2;
					global->ui_insert_count[fid1]+=2;
					break;
				case MULTI_SHARE_EQUAL:     /* For 2 distinct hits, MULTI_SHARE_EQUAL adds 1 to each independent of mate status */
					global->ui_insert_count[fid0]++;
					global->ui_insert_count[fid1]++;
					break;
				case MULTI_SHARE_PROPORTIONAL:
				{
					/******************************************************************************
					 * 12.07.2021:                                                                *
					 * I used to treat two-ends mapping to different sequences as unique hits.    *
					 * But after some thought, I decided to change this behavior to treating them *
					 * as ambiguous hits.                                                         *
					 ******************************************************************************/
					mIVector *mappers = (mIVector*) mMalloc(sizeof(mIVector)); /* Will be freed by main */
					mInitIVector(mappers, 2);
					mPushIVector(mappers, fid0);
					mPushIVector(mappers, fid1);
					mPushVector(global->multi_mappers, mappers);
					break;
				}
				default:
					mDie("Do not understand share_type=%d", share_type);
					break;
			}
			return;
		}

		default:
		{
			uint8_t  *ub_target_hit = global->ub_target_hit;
			mIVector *mappers = (mIVector*) mMalloc(sizeof(mIVector)); /* Will be freed by main */

			global->multi_mapper_count++;

			/* Make a list of all targets hit by this insert */
			mInitIVector(mappers, size);
			for (i=0; i<size; i++) {
				int fid = global->fmap[elem[i]->core.tid];
				if (!ub_target_hit[fid]) {
					mPushIVector(mappers, fid);
					ub_target_hit[fid] = 1;
				}
			}

			/* Share or add the value 1 (for d_insert_count) or 2 (for ui_insert_count) across all mappers depending on the share type */
			switch(share_type) {
				case MULTI_ADD_ALL:
					for (i=0; i<mappers->size; i++) {
						global->ui_insert_count[mappers->elem[i]] += 2;
					}
					break;

				case MULTI_SHARE_EQUAL:
				{
					double share = 1.0/mappers->size;
					for (i=0; i<mappers->size; i++) {
						global->d_insert_count[mappers->elem[i]] += share;
					}
					break;
				}

				case MULTI_SHARE_PROPORTIONAL:
					mPushVector(global->multi_mappers, mappers);
					break;

				default:
					mDie("Do not understand share_type=%d", share_type);
					break;
			}

			/* Reset target_hit flags to 0 */
			for (i=0; i<mappers->size; i++) ub_target_hit[mappers->elem[i]] = 0;

			/* Free mappers IVector */
			if (share_type != MULTI_SHARE_PROPORTIONAL) { /* For MULTI_SHARE_PROPORTIONAL, will be freed by main */
				mFreeIVector(mappers);
				mFree(mappers);
			}
		}
	}
}

/* Estimate abundance on a given bamfile, by iterating through the alignment groups. */

int mEstimateInsertCountOnFile(samfile_t *input, int share_type) {

	int       pool_limit   = 64;
	uint32_t  insert_count = 0;

	bam1_t   *current;
	mBamPool *pool      = (mBamPool*) mCalloc(1, sizeof(mBamPool));
	char     *prev_read = (char*) mCalloc(128, sizeof(char));

	/* init pool */

	mInitBamPool(pool, pool_limit);
	current = pool_current(pool);

	/* In all other apps we also differentiate between mate-pairs. */
	/* But here, mate-pairs are considered together as inserts.    */
	/* So no need to check for mate-pair differences. Only qnames. */
	prev_read[0] = '\0';
	while (samread(input, current) >= 0) {
		if (current->core.tid == -1) {
			continue;
		}
		if ( prev_read[0] != '\0' && strcmp(bam1_qname(current), prev_read) != 0) {
			mEstimateInsertCountOnPool(pool, share_type);
			mReOriginateBamPool(pool);
			current = pool_current(pool);
			insert_count++;
		}
		strncpy(prev_read, bam1_qname(current), 127);
		current = mAdvanceBamPool(pool);
	}
	mEstimateInsertCountOnPool(pool, share_type);
	insert_count++;

	mFreeBamPool(pool);
	mFree(pool);
	mFree(prev_read);
	return insert_count;
}

/* We will always reserve the first column in abundance for unmapped.
 * For normal calculations to work properly, array will be hacked 
 * so that 0th element is a real sequence by doing ++ on array */
mMatrix* mInsertCountToAbundanceMatrix(int row, const char* label, int share_type) {
	int        i;
	int        n_features      = global->n_features;
	char     **feature_name    = global->feature_name;
	double    *abundance       = (double*) mMalloc(n_features*sizeof(double));
	mVector   *multi_mappers   = global->multi_mappers;

	mMatrix  *m       = (mMatrix*) mMalloc(sizeof(mMatrix));
	mInitMatrix(m, 1, n_features+1);

	m->row_names[row] = (char*) mCalloc(strlen(label)+1, sizeof(char));
	strcpy(m->row_names[row], label);

	/* Set Unknown */

	m->col_names[0] = (char*) mMalloc((1+strlen("Unknown"))*sizeof(char));
	strcpy(m->col_names[0], "Unknown");

	/* Hack to hide Unknown (index 0) in the arrays */

	m->col_names++;
	m->elem[row]++;
	m->ncols--;

	/* End hack */

	/* Copy feature names into result matrix */
	for (i=0; i<n_features; i++) {
		m->col_names[i] = (char*) mMalloc((strlen(feature_name[i])+1)*sizeof(char));
		strcpy(m->col_names[i], feature_name[i]);
	}

	/* 1. Get abundance based on uniquely mapped reads -- U(i) for each refseq i */

	/* Get insert_counts */

	for (i=0; i<n_features; i++) {
		/* Since we used 2 for each insert, now let us divide by 2 */
		abundance[i] = 1.0 * global->ui_insert_count[i] / 2;
		/* also reset insert_count so that it can be reused for the next sample */
		global->ui_insert_count[i] = 0;
	}

	switch(share_type) {

		/* Every hit is counted. They were already in abundance. Just copy it */
		case MULTI_ADD_ALL:
			memcpy(m->elem[row], abundance, n_features*sizeof(double));
			break;

		/* Unique hits are already in abundance. Now add the equally-shared hits */
		case MULTI_SHARE_EQUAL:
		{
			/* Add the MULTI_SHARE_EQUAL double values in the d_insert_count array */
			for (i=0; i<n_features; i++) {
				/* Since we used only 1 for each insert in double mode, no need to divide by 2 */
				abundance[i] += global->d_insert_count[i];
				/* also reset insert_count so that it can be reused for the next sample */
				global->d_insert_count[i] = 0;
			}
			memcpy(m->elem[row], abundance, n_features*sizeof(double));
			break;
		}

		/* Unique hits are already in abundance. Now add the proportionally-shared hits */
		/********
		 * 2. Iteratively adjust U(i) by going through all multi-mappers and sharing them between hits 
		 ********/
		case MULTI_SHARE_PROPORTIONAL: /* Proportional sharing! */ 
		{
			int j, k;
			double *abundance_t_k         = (double*) mMalloc(n_features*sizeof(double));
			double *abundance_t_k_minus_1 = (double*) mMalloc(n_features*sizeof(double));

			/* 2.1. Initialize Abundance at t(k), denoted a(i, k), with U(i): *
			 *                     a(i, k) = U(i), for each i                 */

			memcpy(abundance_t_k, abundance, n_features*sizeof(double));

			/* 2.2. Iterate for max 20 times */

			fprintf(stderr, "# Start PropSharing:\n"); 
			for (k=1; k<20; k++) {
				double delta = 0;
				double *increment = (double*) mMalloc(n_features*sizeof(double));

				for (j=0; j<n_features; j++) increment[j] = 0.0f;

				/* 2.2.1. In a new iteration, a(i,k-1) <-- a(i,k) */
				memcpy(abundance_t_k_minus_1, abundance_t_k, n_features*sizeof(double));

				/* 2.2.2. Go through all multi-mappers */
				for (j=0; j<multi_mappers->size; j++) {
					mIVector* multi = (mIVector*) multi_mappers->elem[j];
					int* elem = multi->elem;

					/* 2.2.2.1. Get S = Sigma{ a(i,k) } for all the hits of this multimapper */
					double sum = 0;
					for (i=0; i<multi->size; i++) {
						sum += abundance_t_k[elem[i]];
					}

					/* 2.2.2.2. Share this multi-mapper proportionately according to a(i,k)/S
					 *          Let F(i) = Fraction(read belonging to refseq i)
					 *          It will depend on the abundance of i
					 *                      a(i)      
					 *          F(i) =  -------------
					 *                  Sigma{ a(i) }
					 */
					if (sum > 0) {
						for (i=0; i<multi->size; i++) {
							/* Division is to get the weighted fraction. */
							/* delta(i) += a_t(k) / Sigma(a_t(k)) / */
							increment[elem[i]] += (abundance_t_k[elem[i]] / sum);
						}
					}
				}

				/* Did this iteration do anything useful? */
				delta = 0;
				for (j=0; j<n_features; j++) {
					double diff;
					abundance_t_k[j] = abundance[j] + increment[j];
					/* In my observations, very low numbers are just movements 
					 * towards 0 over many iterations. Better help them converge. */
					if (abundance_t_k[j] < 1e-20) {
						abundance_t_k[j] = 0;
					}
					diff = abundance_t_k[j] - abundance_t_k_minus_1[j];
					delta += diff*diff;
				}
				delta /= n_features;
				fprintf(stderr, "#     PropSharing Iteration: %2d; DELTA^2=%g", k, delta); 
				mFree(increment);
				if (delta < 1e-10) {
					fprintf(stderr, ". CONVERGED!\n");
					break;
				} else {
					fprintf(stderr, "\n");
				}
			}
			fprintf(stderr, "# End   PropSharing!\n"); 
			memcpy(m->elem[row], abundance_t_k, n_features*sizeof(double));

			/* Purge multimapped reads that mapped to features that had no unique reads */
			for (j=0; j<multi_mappers->size; j++) {
				mIVector* multi = (mIVector*) multi_mappers->elem[j];
				int* elem = multi->elem;
				double sum = 0;
				for (i=0; i<multi->size; i++) {
					sum += abundance_t_k[elem[i]];
				}
				if (sum == 0) {
					global->purged_insert_count++;
				}
			}
			fprintf(stderr, "# Purged %d inserts that mapped to features without unique inserts.\n", global->purged_insert_count); 

			mFree(abundance_t_k);
			mFree(abundance_t_k_minus_1);
			break;
		}

		default:
			mDie("Do not understand share_type=%d", share_type);
			break;
	}

	/* Unhide the hidden row and set Unknown abundance to 0 */
	m->elem[row]--;
	m->col_names--;
	m->ncols++;
	m->elem[row][0] = 0.0;

	mFree(abundance);
	return(m);
}

int mCompareTemplate(const void *a, const void *b) {
	const int *ia = (const int *)a;
	const int *ib = (const int *)b;
	int ret = strcmp(global->header->target_name[*ia], global->header->target_name[*ib]);
	return ret;
}

#define subprogram "profile"

int msam_profile_main(int argc, char* argv[]) {

	/* common variables */

	const char      *infile;
	const char      *inmode;
	samfile_t       *input  = NULL;
	mMatrix         *abundance;
	int              share_type    = -1;
	int              unit_type     = -1;
	int              length_normalize = 1;
	int              total_inserts = -1;
	int              mapped_inserts = 0;
	int              this_sample = 0; /* index of this sample in the matrix is 0 */
	int              n_targets;
	FILE            *def_stream;
	char             line[LINE_MAX];

	/* argtable related */

	void           **argtable;
	struct arg_lit  *arg_samin;
	struct arg_file *arg_samfile;
	struct arg_lit  *arg_help;
	struct arg_end  *end;

	struct arg_str  *arg_out;
	struct arg_str  *arg_label;
	struct arg_str  *arg_genome;
	struct arg_int  *arg_total;
	struct arg_int  *arg_mincount;
	struct arg_str  *arg_unit;
	struct arg_lit  *arg_skiplen;
	struct arg_str  *arg_multi;
	int              set_argcount = 0;

	/* Program related */

	int i;
	zoeHash sequences = NULL;
	zoeHash genomes = NULL;
	zoeTVec keys;

	gzFile  output = NULL;

	arg_samin        = arg_lit0("S",  NULL,                     "input is SAM (default: false)");
	arg_samfile      = arg_filen(NULL, NULL, "<bamfile>", 1, 1, "input SAM/BAM file");
	arg_help         = arg_lit0(NULL, "help",                   "print this help and exit\n\n"
								    "Specific options:\n"
								    "-----------------\n");

	/* Specific args */
	arg_out             = arg_str1("o",  NULL,       "<file>", "name of output file (required)");
	arg_label           = arg_str1(NULL, "label",    NULL,     "label to use for the profile; typically the sample id (required)");
	arg_genome          = arg_str0(NULL, "genome",   NULL,     "tab-delimited genome definition file - 'genome-id<tab>seq-id' (default: none)");
	arg_mincount        = arg_int0(NULL, "mincount", NULL,     "minimum number of inserts mapped to a feature, below which the feature is counted as absent (default: 0)");
	arg_total           = arg_int0(NULL, "total",    NULL,     "number of high-quality inserts (mate-pairs/paired-ends) that were input to the aligner (default: 0)");
	arg_unit            = arg_str0(NULL, "unit",     NULL,     "unit of abundance to report {ab | rel | fpkm | tpm} (default: rel)");
	arg_skiplen         = arg_lit0(NULL, "nolen",              "do not normalize the abundance (only relevant for ab or rel) for sequence length (default: normalize)");
	arg_multi           = arg_str0(NULL, "multi",    NULL,     "how to deal with multi-mappers {all | equal | proportional} (default: proportional)\n"
                                                                 "\n"
                                                                 "Description\n"
                                                                 "-----------\n"
                                                                 "\n"
                                                                 "Produces an abundance profile of all reference sequences in a BAM file\n"
                                                                 "based on the number of read-pairs (inserts) mapping to each reference sequence.\n"
                                                                 "It can work with genome-scale reference sequences while mapping to a database \n"
                                                                 "of sequenced genomes, but can also work with gene-scale sequences such as in the\n"
                                                                 "Integrated Gene Catalog from human gut microbiome (Li et al, Nat biotech 2014).\n"
                                                                 "\n"
                                                                 "In the output file, each sequence in the BAM file gets a line with its abundance.\n"
                                                                 "They are presented in the order in which they appear in the BAM header. <label>\n"
                                                                 "is used as the first line, so that reading or 'joining' these files is easier.\n"
                                                                 "\n"
                                                                 "--total option:      In metagenomics, an unmapped read could still be a valid\n"
                                                                 "                     sequence, just missing in the database being mapped against.\n"
                                                                 "                     This is the purpose of the '--total' option to track the\n"
                                                                 "                     fraction of 'unknown' entities in the sample. If --total\n"
                                                                 "                     is ignored or specified as --total=0, then tracking the \n"
                                                                 "                     'unknown' fraction is disabled. However, if the total \n"
                                                                 "                     sequenced inserts were given, then there will be a new\n"
                                                                 "                     feature added to denote the 'unknown' fraction.\n"
                                                                 "Units of abundance:  Currently four different units are available.\n"
                                                                 "                          'ab': raw insert-count abundance\n"
                                                                 "                         'rel': relative abundance (default)\n"
                                                                 "                        'fpkm': fragments per kilobase of sequence per million reads\n"
                                                                 "                         'tpm': transcripts per million\n"
                                                                 "                     If number of reads input to the aligner is given via --total,\n"
                                                                 "                     fpkm and tpm will behave differently than in RNAseq data,\n"
                                                                 "                     as there is now a new entity called 'unknown'.\n"
                                                                 "Alignment filtering: 'profile' expects that every alignment listed is considered \n"
                                                                 "                     valid. For example, if one needs to filter alignments \n"
                                                                 "                     based on alignment length, read length, alignment percent\n"
                                                                 "                     identity, etc, this should have been done prior to \n"
                                                                 "                     '"subprogram"'. Please see 'filter' for such filtering.\n"
                                                                 "Multi-mapper reads:  Reads mapping to multiple references need to be considered\n"
                                                                 "                     carefully, as spurious mappings of promiscuous regions or\n"
                                                                 "                     short homology could lead to incorrect abundances of \n"
                                                                 "                     sequences. '"subprogram"' offers three options for this purpose.\n"
                                                                 "                     If an insert maps to N references at the same time:\n"
                                                                 "                   'all': each reference gets 1 insert added.\n"
                                                                 "                 'equal': each reference gets 1/N insert added.\n"
                                                                 "          'proportional': each reference gets a fraction proportional to its \n"
                                                                 "                          reference-sequence-length-normalized relative \n"
                                                                 "                          abundance estimated only based on uniquely\n"
                                                                 "                          mapped reads."
                                                                 );
	end    = arg_end(20); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */

	argtable = (void**) mCalloc(12, sizeof(void*));

	/* Common args */
	set_argcount = 0;
	argtable[set_argcount++] = arg_samin;
	argtable[set_argcount++] = arg_samfile;
	argtable[set_argcount++] = arg_help;

	/* Specific args */
	argtable[set_argcount++] = arg_out;
	argtable[set_argcount++] = arg_label;
	argtable[set_argcount++] = arg_genome;
	argtable[set_argcount++] = arg_total;
	argtable[set_argcount++] = arg_mincount;
	argtable[set_argcount++] = arg_unit;
	argtable[set_argcount++] = arg_skiplen;
	argtable[set_argcount++] = arg_multi;
	argtable[set_argcount++] = end;

	/* Set Global */
	mInitGlobal();
	global->multiple_input = 0;

	/* parse command line */

	{
		int nerrors;
		if (arg_nullcheck(argtable) != 0) {
			mDie("insufficient memory");
		}

		/* parse */
		nerrors = arg_parse(argc, argv, argtable);

		/* help needed? */
		if (arg_help->count > 0) {
			mPrintHelp(subprogram, argtable);
			exit(EXIT_SUCCESS);
		}

		/* parse error? */
		if (nerrors > 0) {
			arg_print_errors(stdout, end, PROGRAM);
			fprintf(stdout, "Use --help for usage instructions!\n");
			mQuit("");
		}

		if (arg_samfile->count > 1 && global->multiple_input == 0) {
			mMultipleFileError(subprogram, argtable);
		}
	}

	if (arg_label->count != 1 || arg_out->count != 1) {
		fprintf(stdout, "requires --label and -o\n");
		mPrintHelp(subprogram, argtable);
		mQuit("");
	}

	if (arg_total->count > 0) {
		total_inserts = arg_total->ival[0];
		if (total_inserts < 0) {
			fprintf(stdout, "--total must be a positive integer\n");
			mPrintHelp(subprogram, argtable);
			mQuit("");
		}
	}

	/* Set input/output modes */

	inmode = M_INPUT_MODE(arg_samin);

	infile = arg_samfile->filename[0];
	input = mOpenSamFile(infile, inmode, NULL);
	global->header = input->header;

	/* multi-mapper share type */

	share_type = -1;
	if (arg_multi->count > 0) {
		/* Check with #defines MULTI_ADD_ALL=1; MULTI_SHARE_EQUAL=2; MULTI_SHARE_PROPORTIONAL=3; */
		const char *types[4] = {"", "all", "equal", "proportional"};
		for (i=1; i<=3; i++) {
			/* Any prefix of the options will work. Make sure that the options don't share a prefix. */
			if (strncmp(arg_multi->sval[0], types[i], strlen(arg_multi->sval[0])) == 0) {
				share_type = i;
				break;
			}
		}
		if (share_type == -1) {
			mDie("Do not understand --multi=%s", arg_multi->sval[0]);
		}
	} else {
		share_type = MULTI_SHARE_PROPORTIONAL;
	}

	/* unit of measurement */

	unit_type = -1;
	if (arg_unit->count > 0) {
		/* Check with #defines UNIT_REL=1; UNIT_FPKM=2; UNIT_TPM=3; */
		const char *types[5] = {"", "relative", "fpkm", "tpm", "abundance"};
		for (i=1; i<=4; i++) {
			/* Any prefix of the options will work. Make sure that the options don't share a prefix. */
			if (strncmp(arg_unit->sval[0], types[i], strlen(arg_unit->sval[0])) == 0) {
				unit_type = i;
				break;
			}
		}
		if (unit_type == -1) {
			mDie("Do not understand --unit=%s", arg_unit->sval[0]);
		}
	} else {
		unit_type = UNIT_REL;
	}

	/* Should I normalize for length */

	length_normalize = 1;
	if (unit_type == UNIT_REL || unit_type == UNIT_ABN) {
		length_normalize = (arg_skiplen->count == 0);
	}
	
	n_targets = global->header->n_targets;
	global->fmap = (int*) mMalloc(n_targets*sizeof(int));
	for (i=0; i<n_targets; i++) {
		global->fmap[i] = -1;
	}

	/* Assign seq->genome maps for counting */
	if (arg_genome->count > 0) {
		char *key;
		sequences = zoeNewHash();
		for (i=0; i<n_targets; i++) {
			int *code = (int*) mMalloc(sizeof(int));
			*code = i;
			zoeSetHash(sequences, global->header->target_name[i], code);
		}

		/* Make a hash of genomes from the def file */
		def_stream = mSafeOpenFile(arg_genome->sval[0], "r", 0);
		genomes = zoeNewHash();
		key  = (char*) mMalloc(LINE_MAX*sizeof(char));
		while (fgets(line, LINE_MAX, def_stream) != NULL) {
			int  *code = (int*) mMalloc(sizeof(int));
			char *seqname  = (char*) mMalloc(LINE_MAX*sizeof(char));
			char  *genome_name  = (char*) mMalloc(LINE_MAX*sizeof(char));
			if (sscanf(line, "%s\t%s", genome_name, seqname) != 2) {
				mDie("GENOME DEFINITION LINE ERROR");
			}

			*code = 1;
			zoeSetHash(genomes, genome_name, code);
			mFree(seqname);
		}
		mSafeCloseFile(def_stream, 0);
		keys = zoeKeysOfHash(genomes);
		global->n_features = keys->size;
		for (i=0; i<keys->size; i++) {
			int  *code = (int*) mMalloc(sizeof(int));
			key = keys->elem[i];
			*code = i;
			zoeSetHash(genomes, key, code);
		}
		/* Now genome_name ==> index */

		/* Make a hash of sequence names from the list file */
		def_stream = mSafeOpenFile(arg_genome->sval[0], "r", 0);
		while (fgets(line, LINE_MAX, def_stream) != NULL) {
			int *genome_id;
			int *seq_id;
			char  *genome_name  = (char*) mMalloc(LINE_MAX*sizeof(char));
			char  *seqname = (char*) mMalloc(LINE_MAX*sizeof(char));
			if (sscanf(line, "%s\t%s", genome_name, seqname) != 2) {
				mDie("GENOME DEFINITION LINE ERROR");
			}

			genome_id = (int*) zoeGetHash(genomes, genome_name);
			seq_id    = (int*) zoeGetHash(sequences, seqname);
			if (genome_id == NULL) {
				mDie("Genome '%s' not found in BAM file", genome_name);
			}
			if (seq_id == NULL) {
				mDie("Sequence '%s' not found in BAM file", seqname);
			}
			global->fmap[*seq_id] = *genome_id;
		}
		mSafeCloseFile(def_stream, 0);

		/* Estimate feature lengths */
		global->feature_len = (uint32_t*) mMalloc(global->n_features*sizeof(uint32_t));
		for (i=0; i<global->n_features; i++) {
			global->feature_len[i] = 0;
		}
		for (i=0; i<n_targets; i++) {
			if (global->fmap[i] == -1) {
				mDie("Sequence '%s' not found in genome definition", global->header->target_name[i]);
			}
			global->feature_len[global->fmap[i]] += global->header->target_len[i];
		}

		/* Set feature names */
		keys = zoeKeysOfHash(genomes);
		global->feature_name = (char**) mMalloc(keys->size*sizeof(char*));
		for (i=0; i<keys->size; i++) {
			char *fname = (char*)keys->elem[i];
			global->feature_name[i] = (char*) mMalloc((1+strlen(fname))*sizeof(char));
			strcpy(global->feature_name[i], fname);
		}

	} else {
		for (i=0; i<n_targets; i++) {
			global->fmap[i] = i;
		}
		global->n_features   = n_targets;
		global->feature_name = global->header->target_name;
		global->feature_len  = global->header->target_len;
	}

	/* Calculate abundance */
	mInitInsertCounts(share_type);
	mapped_inserts = mEstimateInsertCountOnFile(input, share_type);
	abundance      = mInsertCountToAbundanceMatrix(this_sample, arg_label->sval[0], share_type); /* Make abundance table */
	if (arg_mincount->count > 0) {
		int    mincount = arg_mincount->ival[0];
		double purged_inserts = 0;
		/* Mask features with fewer than min_count reads */
		/* We do that by moving inserts from them to Unknown */
		for (i=1; i<abundance->ncols; i++) {
			if (abundance->elem[this_sample][i] < mincount) {
				purged_inserts += abundance->elem[this_sample][i];
				abundance->elem[this_sample][i] = 0;
			}
		}
		purged_inserts = round(purged_inserts); /* Round it */
		fprintf(stderr, "# Purged %d inserts from low-abundance features based on --mincount.\n", (int)purged_inserts);
		global->purged_insert_count += (int)purged_inserts;
	}

	/* Introduce unmapped if necessary */

	if (total_inserts > 0 && total_inserts < mapped_inserts) {
		fprintf(stderr, "# Ignoring 'unknown' fraction, as total inserts (%d) < mapped inserts (%d)!\n", total_inserts, mapped_inserts);
		total_inserts = -1;
	}

	/* Open output file */
	if (strcmp(arg_out->sval[0], "-") == 0) { /* If '-' was given as output, redirect to stdout */
		output = gzdopen(fileno(stdout), "wb");
	} else {
		output = gzopen(arg_out->sval[0], "wb");
	}

	/* Print command-line for book-keeping */
	mPrintCommandLineGzip(output, argc, argv);

	if (total_inserts > 0) {
		gzprintf(output, "#   Total inserts: %d\n", total_inserts);
		gzprintf(output, "#  Mapped inserts: %d (%.2f%%)\n", mapped_inserts - global->purged_insert_count, 100.0*(mapped_inserts - global->purged_insert_count)/total_inserts);
		if (arg_mincount->count > 0) {
			gzprintf(output, "# Uniquely mapped: NA (NA%%)\n");
			gzprintf(output, "# Multiple mapped: NA (NA%%)\n");
		} else {
			gzprintf(output, "# Uniquely mapped: %d (%.2f%%)\n", global->uniq_mapper_count, 100.0*global->uniq_mapper_count/total_inserts);
			gzprintf(output, "# Multiple mapped: %d (%.2f%%)\n", global->multi_mapper_count, 100.0*global->multi_mapper_count/total_inserts);
		}
	} else {
		gzprintf(output, "#   Total inserts: NA\n");
		gzprintf(output, "#  Mapped inserts: %d (NA%%)\n", mapped_inserts - global->purged_insert_count);
		gzprintf(output, "# Uniquely mapped: %d (NA%%)\n", global->uniq_mapper_count);
		gzprintf(output, "# Multiple mapped: %d (NA%%)\n", global->multi_mapper_count);
		gzprintf(output, "# Estimated seq. length for 'Unknown': NA\n");
	}

	/* Create new feature for 'unknown' if total_inserts is valid */
	if (total_inserts > 0) {
		int       i;

		uint32_t *feature_len = global->feature_len;

		/* We just add "purged inserts" back to Unknown */
		abundance->elem[this_sample][0] = total_inserts - mapped_inserts + global->purged_insert_count;

		if (length_normalize) {
			/* Remember: abundance has elem[0] as Unknown, which is missing in target_len[] */
			int      count = 0;
			uint64_t sum   = 0;
			uint32_t unknown_size;
			for (i=0; i<global->n_features; i++) {
				sum += feature_len[i];
				count++;
			}
			unknown_size = sum / count;
			gzprintf(output, "# Estimated seq. length for 'Unknown': %dbp\n", unknown_size);
			abundance->elem[this_sample][0] = 1.0 * abundance->elem[this_sample][0] / unknown_size ;
		} else {
			gzprintf(output, "# Estimated seq. length for 'Unknown': NA\n");
		}
	}

	/* Length normalize */
	if (length_normalize) {
		abundance->col_names++;
		abundance->elem[this_sample]++;
		abundance->ncols--;
		for (i=0; i<global->n_features; i++) {
			abundance->elem[this_sample][i] /= global->feature_len[i];
		}
		abundance->col_names--;
		abundance->elem[this_sample]--;
		abundance->ncols++;
	}

	/* Convert length-normalized abundance into the right unit */
	switch(unit_type) {
		case UNIT_FPKM:  /* Convert to FPKM */
			/* We have fragments divided by seq length */
			/* So multiply by 1000 to get per kb seq length */
			/* To get per million reads, divide by total/10^6 */
			/* Altogether, *10^9/total */
			if (total_inserts > 0) {
				mMultiplyMatrixByScalar(abundance, 1.0E9/total_inserts);
			} else {
				mMultiplyMatrixByScalar(abundance, 1.0E9/mapped_inserts);
			}
			break;

		case UNIT_TPM:  /* Normalize to relative abundance and them multiple by 1E06 */
			mRowNormalizeMatrix(abundance);
			mMultiplyMatrixByScalar(abundance, 1.0E6);
			break;

		case UNIT_REL:  /* Normalize abundance to relative abundance */
			mRowNormalizeMatrix(abundance);
			break;

		case UNIT_ABN:  /* Do nothing */
		default:
			break;
	}

	/* Write output */
	mWriteRMatrixTransposedGzip(output, abundance);

	/* Close output */
	gzclose(output);

	/* Wind-up operations */

	mFreeInsertCounts(share_type);
	mFreeMatrix(abundance);
	mFree(abundance);

	/* Free the hashes */
	if (arg_genome->count > 0) {
		keys = zoeKeysOfHash(genomes);
		for (i=0; i<keys->size; i++) {
			mFree(zoeGetHash(genomes, keys->elem[i]));
		}
		zoeDeleteTVec(keys);
		zoeDeleteHash(genomes);
		keys = zoeKeysOfHash(sequences);
		for (i=0; i<keys->size; i++) {
			mFree(zoeGetHash(sequences, keys->elem[i]));
		}
		zoeDeleteTVec(keys);
		zoeDeleteHash(sequences);
	}

	samclose(input);
	arg_freetable(argtable, set_argcount);
	mFree(argtable);

	mFreeGlobal();

	return 0;
}
