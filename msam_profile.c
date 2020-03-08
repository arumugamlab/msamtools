#include "msam.h"

#define ADD_ALL (1)
#define SHARE_EQUAL (2)
#define SHARE_PROPORTIONAL (3)

void mInitInsertCounts(int share_type);
int  mEstimateInsertCountOnFile(samfile_t *input, int share_type);
void mWriteCompressedSeqAbundance();
mMatrix* mInsertCountToAbundanceMatrix(int row, const char* rowname, int share_type);

void mInitInsertCounts(int share_type) {
	int n_targets  = global->header->n_targets;

	global->ui_insert_count = (uint32_t*) mCalloc(n_targets, sizeof(uint32_t));
	global->ub_target_hit   = (uint8_t*) mCalloc(n_targets, sizeof(uint8_t));

	if (share_type == SHARE_EQUAL) {
		int i;
		double *f = (double*) mMalloc(n_targets*sizeof(double));
		for (i=0; i<n_targets; i++) f[i] = 0.0f;
		global->d_insert_count = f;
	}

	if (share_type == SHARE_PROPORTIONAL) {
		global->multi_mappers = (mVector*) mMalloc(sizeof(mVector));
		mInitVector(global->multi_mappers, 65536);
	}
}

void mFreeInsertCounts(int share_type) {
	if (share_type == SHARE_PROPORTIONAL) {
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

	if (share_type == SHARE_EQUAL) {
		mFree(global->d_insert_count);
	}

	mFree(global->ui_insert_count);
	mFree(global->ub_target_hit);
}

/* This mode makes sure that the total insert counts are preserved in the end */
/* It shares a total count of 1 among the seen sequences */
void mEstimateInsertCountOnPool(mBamPool *pool, int share_type) {
	bam1_t **elem  = pool->elem;
	int      size  = pool->size;
	int      i;

/* Note the multiple return's in this function for avoiding too many conditionals.
 * Maintaining might be difficult, but I chose this over if-else-if-else-etc */

	switch(size) {
		case 1: /* With one hit, you got to be unique hit, so add 2 no matter what share_type */
			global->ui_insert_count[elem[0]->core.tid] += 2;
			return;

		case 2: /* could be 2 mate-pairs, or 1 mate mapped twice*/
		{
			int tid0  = elem[0]->core.tid;
			int tid1  = elem[1]->core.tid;

			if (tid0 == tid1) { /* Hits twice the same target, so add 2 independent of share_type and mate status */
				global->ui_insert_count[tid0] += 2;
				return;
			}

			/* If we are here, tid0 and tid1 are different, so multi-mappers */
	/* TO DO: Should I differentiate between mate-types? Probably not! */
			switch(share_type) {
				case ADD_ALL:     /* For 2 distinct hits, SHARE_EQUAL adds 1 to each independent of mate status */
					global->ui_insert_count[tid0]+=2;
					global->ui_insert_count[tid1]+=2;
					break;
				case SHARE_EQUAL:     /* For 2 distinct hits, SHARE_EQUAL adds 1 to each independent of mate status */
					global->ui_insert_count[tid0]++;
					global->ui_insert_count[tid1]++;
					break;
				case SHARE_PROPORTIONAL:
				{
					int fmate = (BAM_FREAD1 | BAM_FREAD2); /* Test which mates are these, so capture these two flags */
					int mate1 = elem[0]->core.flag & fmate;
					int mate2 = elem[1]->core.flag & fmate;
					if (mate1 != mate2) { /* Unique    : Mate-pairs, so this is still unique. Add 1 each for SHARE_PROPORTIONAL and SHARE_EQUAL*/
						global->ui_insert_count[tid0]++;
						global->ui_insert_count[tid1]++;
					} else {              /* Non-unique: Same read twice; works also for mate1=mate2=0 when mapped in unpaired mode */
						mIVector *mappers = (mIVector*) mMalloc(sizeof(mIVector)); /* Will be freed by main */
						mInitIVector(mappers, 2);
						mPushIVector(mappers, tid0);
						mPushIVector(mappers, tid1);
						mPushVector(global->multi_mappers, mappers);
					}
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
			mInitIVector(mappers, size);
			for (i=0; i<size; i++) {
				int tid = elem[i]->core.tid;
				if (!ub_target_hit[tid]) {
					mPushIVector(mappers, tid);
					ub_target_hit[tid] = 1;
				}
			}
			switch(share_type) {
				case ADD_ALL:
					for (i=0; i<mappers->size; i++) {
						global->ui_insert_count[mappers->elem[i]] += 2;
					}
					break;

				case SHARE_EQUAL:
				{
					double share = 1.0f/mappers->size;
					for (i=0; i<mappers->size; i++) {
						global->d_insert_count[mappers->elem[i]] += share;
					}
					break;
				}

				case SHARE_PROPORTIONAL:
					mPushVector(global->multi_mappers, mappers);
					break;

				default:
					mDie("Do not understand share_type=%d", share_type);
					break;
			}

			/* Reset target_hit flags to 0 */
			for (i=0; i<mappers->size; i++) ub_target_hit[mappers->elem[i]] = 0;

			/* Free mappers IVector */
			if (share_type != SHARE_PROPORTIONAL) { /* For SHARE_PROPORTIONAL, will be freed by main */
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
	bam1_t   *prev         = NULL;
	mBamPool *pool         = (mBamPool*) mCalloc(1, sizeof(mBamPool));

	/* init pool */

	mInitBamPool(pool, pool_limit);
	current = pool_current(pool);

	while (samread(input, current) >= 0) {
		if (prev != NULL && (bam1_templatecmp(current, prev) != 0)) {
			mEstimateInsertCountOnPool(pool, share_type);
			mReOriginateBamPool(pool);
			current = pool_current(pool);
			insert_count++;
		}
		prev    = current;
		current = mAdvanceBamPool(pool);
	}
	mEstimateInsertCountOnPool(pool, share_type);
	insert_count++;

	mFreeBamPool(pool);
	mFree(pool);
	/*fprintf(stderr, "Insert Count=%d\n", insert_count);*/
	return insert_count;
}

/* We will always reserve the first column in abundance for unmapped.
 * For normal calculations to work properly, array will be hacked 
 * so that 0th element is a real sequence by doing ++ on array */
mMatrix* mInsertCountToAbundanceMatrix(int row, const char* label, int share_type) {
	int        i;
	int        n_targets       = global->header->n_targets;
	uint32_t  *target_len      = global->header->target_len;
	uint32_t  *ui_insert_count = global->ui_insert_count;
	double    *abundance       = (double*) mMalloc(n_targets*sizeof(double));
	double    *normalize       = (double*) mMalloc(n_targets*sizeof(double));
	char     **target_name     = global->header->target_name;
	mVector   *multi_mappers   = global->multi_mappers;

	mMatrix  *m       = (mMatrix*) mMalloc(sizeof(mMatrix));
	mInitMatrix(m, 1, n_targets+1);

	m->row_names[row] = (char*) mCalloc(strlen(label)+1, sizeof(char));
	strcpy(m->row_names[row], label);

	/* Set Unknown */

	m->col_names[0] = (char*) mMalloc((1+strlen("Unknown"))*sizeof(char));
	strcpy(m->col_names[0], "Unknown");

	/* Hack to hide the actual first (index 0) element in the arrays */

	m->col_names++;
	m->elem[row]++;
	m->ncols--;

	/* End hack */

	for (i=0; i<n_targets; i++) {
		m->col_names[i] = (char*) mMalloc((strlen(target_name[i])+1)*sizeof(char));
		strcpy(m->col_names[i], target_name[i]);
	}

	/* normalization factor for length normalization */
	/* instead of dividing by length every time, we divide once and multiply it every time */

	for (i=0; i<n_targets; i++) {
		normalize[i] = 1.0f / target_len[i];
	}

	/* 1. Get abundance based on uniquely mapped reads -- U(i) for each refseq i */

	/* Get insert_counts and length-normalize */

	for (i=0; i<n_targets; i++) {
		/* Since we used 2 for each insert, now let us divide by 2 */
		abundance[i] = 0.5 * ui_insert_count[i] * normalize[i];
		/* also reset insert_count so that it can be reused */
		ui_insert_count[i] = 0;
	}

	switch(share_type) {

		case ADD_ALL:
			memcpy(m->elem[row], abundance, n_targets*sizeof(double));
			break;

		case SHARE_EQUAL:
		{
			double *d_insert_count = global->d_insert_count;
			/* Add the SHARE_EQUAL double values in the d_insert_count array */
			for (i=0; i<n_targets; i++) {
				/* Since we used only 1 for each insert in double mode, no need to divide by 2 */
				abundance[i] += d_insert_count[i] * normalize[i];
				/* also reset insert_count so that it can be reused */
				d_insert_count[i] = 0;
			}
			memcpy(m->elem[row], abundance, n_targets*sizeof(double));
			break;
		}

		/********
		 * 2. Iteratively adjust U(i) by going through all multi-mappers and sharing them between hits 
		 ********/
		case SHARE_PROPORTIONAL: /* Proportional sharing! */ 
		{
			int j, k;
			double *abundance_t_k         = (double*) mMalloc(n_targets*sizeof(double));
			double *abundance_t_k_minus_1 = (double*) mMalloc(n_targets*sizeof(double));

			/* 2.1. Initialize Abundance at t(k), denoted a(i, k), with U(i): *
			 *                     a(i, k) = U(i), for each i                 */

			memcpy(abundance_t_k, abundance, n_targets*sizeof(double));

			/* 2.2. Iterate for max 20 times */

			fprintf(stderr, "# Start PropSharing:\n"); 
			for (k=1; k<20; k++) {
				double delta = 0;
				double *increment = (double*) mMalloc(n_targets*sizeof(double));

				for (j=0; j<n_targets; j++) increment[j] = 0.0f;

				/* 2.2.1. In a new iteration, a(i,k-1) <-- a(i,k) */
				memcpy(abundance_t_k_minus_1, abundance_t_k, n_targets*sizeof(double));

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
					 *          It will depend on the abundance of i, and also its length
					 *                      a(i)            1
					 *          F(i) =  -------------  x  ------
					 *                  Sigma{ a(i) }     len(i)
					 */
					if (sum > 0) {
						for (i=0; i<multi->size; i++) {
							/* Division is to get the weighted fraction. Length-normalization is done by multiplying */
							/* delta(i) += a_t(k) / Sigma(a_t(k)) / */
							increment[elem[i]] += (normalize[elem[i]] * abundance_t_k[elem[i]] / sum);
						}
					} else {
						for (i=0; i<multi->size; i++) {
							increment[elem[i]] += (normalize[elem[i]] / multi->size);
						}
					}
				}

				/* Did this iteration do anything useful? */
				delta = 0;
				for (j=0; j<n_targets; j++) {
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
				delta /= n_targets;
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
			memcpy(m->elem[row], abundance_t_k, n_targets*sizeof(double));
			mFree(abundance_t_k);
			mFree(abundance_t_k_minus_1);
			break;
		}

		default:
			mDie("Do not understand share_type=%d", share_type);
			break;
	}

	mFree(abundance);
	mFree(normalize);
	return(m);
}

int mCompareTemplate(const void *a, const void *b) {
	const int *ia = (const int *)a;
	const int *ib = (const int *)b;
	int ret = strcmp(global->header->target_name[*ia], global->header->target_name[*ib]);
	return ret;
}

void mWriteCompressedSeqAbundance(mMatrix *m_abundance) {
/*
	mWriteUIntArrayBinary(stream, header->n_targets, ui_insert_count);
	int i;
	uint32_t size;
	FILE *stream = mSafeOpenFile("out.bin", "wb", 0);
	mWriteUIntArrayBinary(stream, header->n_targets, ui_insert_count);
	mSafeCloseFile(stream, 0);
	mSafeOpenFile("out.bin", "rb", 0);
	mFree(ui_insert_count);
	mReadUIntArrayBinary(stream, &size, &ui_insert_count);
*/
	int     sorted = 0;
	int     i;
	int     row_id = 0; /* If we ever bring back multi-bam processing, this would be useful */
	FILE   *compressor = fdopen(global->pipe_fd[1], "w");
	int     n_targets = global->header->n_targets;
	double *abundance = m_abundance->elem[row_id];
	int    *indices = mCalloc(n_targets,sizeof(int));

	for (i=0; i<n_targets; i++) indices[i] = i;
	if (sorted == 1) {
		qsort(indices, n_targets, sizeof(int), mCompareTemplate);
	}

	fprintf(compressor, "%s\n", m_abundance->row_names[row_id]);
	for (i=0; i<n_targets; i++) {
		fprintf(compressor, "%s\t", global->header->target_name[indices[i]]);
		fprintf(compressor, "%g", abundance[indices[i]]);
		fprintf(compressor, "\n");
	}
	fclose(compressor);
	mFree(indices);
}

#define subprogram "profile"

int msam_profile_main(int argc, char* argv[]) {

	/* common variables */

	char             outfile[256];
	const char      *infile;
	const char      *inmode;
	samfile_t       *input  = NULL;
	mMatrix         *abundance;
	int              share_type = -1;
	int              total_inserts = -1;
	int              mapped_inserts = 0;
	int              this_sample = 0; /* index of this sample in the matrix is 0 */

	/* argtable related */

	void           **argtable;
	struct arg_lit  *arg_samin;
	struct arg_file *arg_samfile;
	struct arg_lit  *arg_help;
	struct arg_end  *end;

	struct arg_str  *arg_out;
	struct arg_str  *arg_label;
	struct arg_str  *arg_multi;
	struct arg_int  *arg_total;
	struct arg_lit  *arg_gzip;
	int              set_argcount = 0;

	/* Program related */

	int              gzip = 0;
	FILE            *output = NULL;

	mInitGlobal();
	global->multiple_input = 0;

	arg_samin        = arg_lit0("S",  NULL,                     "input is SAM (default: false)");
	arg_samfile      = arg_filen(NULL, NULL, "<bamfile>", 1, 1, "input SAM/BAM file");
	arg_help         = arg_lit0(NULL, "help",                   "print this help and exit\n\n"
								    "Specific options:\n"
								    "-----------------\n");

	/* Specific args */
	arg_out             = arg_str1("o",  NULL,     "<file>", "name of output file (required)");
	arg_label           = arg_str1(NULL, "label",  NULL,     "label to use for the profile; typically the sample id (required)");
	arg_multi           = arg_str0(NULL, "multi",  NULL,     "how to deal with multi-mappers {all | equal | proportional} (default: proportional)");
	arg_total           = arg_int0(NULL, "total",  NULL,     "number of high-quality inserts (mate-pairs/paired-ends) that were input to the aligner (default: 0)");
	arg_gzip            = arg_lit0("z", "gzip",              "compress output file using gzip (default: false)\n"
                                                                 "\n"
                                                                 "Description\n"
                                                                 "-----------\n"
                                                                 "\n"
                                                                 "Produces a relative abundance profile of all reference sequences in a BAM file\n"
                                                                 "based on the number of read-pairs (inserts) mapping to each reference sequence.\n"
                                                                 "It can work with genome-scale reference sequences while mapping to a database \n"
                                                                 "of sequenced genomes, but can also work with gene-scale sequences such as in the\n"
                                                                 "Integrated Gene Catalog from human gut microbiome (Li et al, Nat biotech 2014).\n"
                                                                 "\n"
                                                                 "In the output file, each sequence in the BAM file gets a line with its abundance.\n"
                                                                 "They are presented in the order in which they appear in the BAM header. <label>\n"
                                                                 "is used as the first line, so that reading or 'joining' these files is easier.\n"
                                                                 "\n"
                                                                 "If using '-z', output file does NOT automatically get '.gz' extension. This is \n"
                                                                 "up to the user to specify the correct full output file name.\n"
                                                                 "\n"
                                                                 "Alignment filtering: It expects that every alignment listed is considered \n"
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
	end    = arg_end(9); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */

	argtable = (void**) mCalloc(9, sizeof(void*));

	/* Common args */
	set_argcount = 0;
	argtable[set_argcount++] = arg_samin;
	argtable[set_argcount++] = arg_samfile;
	argtable[set_argcount++] = arg_help;

	/* Specific args */
	argtable[set_argcount++] = arg_out;
	argtable[set_argcount++] = arg_label;
	argtable[set_argcount++] = arg_multi;
	argtable[set_argcount++] = arg_total;
	argtable[set_argcount++] = arg_gzip;
	argtable[set_argcount++] = end;

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
			mQuit("");
		}

		/* parse error? */
		if (nerrors > 0) {
			arg_print_errors(stderr, end, PROGRAM);
			fprintf(stderr, "Use --help for usage instructions!\n");
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

	/* General operations */

	infile = arg_samfile->filename[0];
	fprintf(stderr, "%s\n", infile);
	input = mOpenSamFile(infile, inmode, NULL);
	strcpy(outfile, arg_out->sval[0]);
	global->header = input->header;

	/* Specific operations */

	if (arg_gzip->count > 0) {
		gzip = 1;
	}

	share_type = -1;
	if (arg_multi->count > 0) {
		/* Check with #defines ADD_ALL=1; SHARE_EQUAL=2; SHARE_PROPORTIONAL=3; */
		const char *types[4] = {"", "all", "equal", "proportional"};
		int         i;
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
		share_type = SHARE_PROPORTIONAL;
	}

	/* Calculate abundance */
	mInitInsertCounts(share_type);
	mapped_inserts = mEstimateInsertCountOnFile(input, share_type);
	abundance      = mInsertCountToAbundanceMatrix(this_sample, arg_label->sval[0], share_type); /* Normalize for sequence length and make abundance */

	/* Introduce unmapped if necessary */

	if (total_inserts > 0 && total_inserts < mapped_inserts) {
		fprintf(stderr, "# Ignoring 'unknown' fraction, as total inserts (%d) < mapped inserts (%d)!\n", total_inserts, mapped_inserts);
		total_inserts = -1;
	}

	if (total_inserts > 0) {
		int       i;
		int       count = 0;
		uint32_t  sum   = 0;

		int       this_sample = 0; /* This sample is row 0 in matrix */
		double   *this_elem = abundance->elem[this_sample];

		uint32_t *target_len = global->header->target_len;
		int       n_targets  = global->header->n_targets;

		int       unmapped   = total_inserts - mapped_inserts;
		uint32_t  unknown_size;

		for (i=0; i<n_targets; i++) {
			if (this_elem[i] > 0) {
				sum += target_len[i];
				count++;
			}
		}

		/* Unhide the hidden row */
		abundance->elem[this_sample]--;
		abundance->col_names--;
		abundance->ncols++;

		unknown_size = sum / count;
		fprintf(stderr, "# Estimated seq. length for 'Unknown': %dbp\n", unknown_size);
		abundance->elem[this_sample][0] = 1.0f * unmapped / unknown_size ;
	}

	/* Normalize to relative abundance */

	mRowNormalizeMatrix(abundance);

	/* Write output */
	output = mInitOutputStream(outfile, gzip);
	mWriteRMatrixTransposed(output, abundance);
	mFreeOutputStream(output, gzip);

	/* Wind-up operations */

	mFreeInsertCounts(share_type);

	/* Unhide hidden element, if it is not hidden already */
	if (total_inserts <= 0) {
		abundance->elem[this_sample]--; 
		abundance->col_names--; 
		abundance->ncols++;
	}
	mFreeMatrix(abundance);
	mFree(abundance);

	samclose(input);
	arg_freetable(argtable, set_argcount);
	mFree(argtable);

	mFreeGlobal();

	return 0;
}
