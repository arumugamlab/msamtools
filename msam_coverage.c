#include "msam.h"

void mEstimateCoverageOnFile(samfile_t *input);
void mInitCoverage();
void mWriteCoverageToStream(FILE* outstream, int skip_uncovered, int wordsize);

/*
 * Note: samtools' bam_calend() function returns a 1-based coordinate instead of
 * a 0-based coordinate that the documentation promises. This is potentially 
 * because they want to easily estimate 
 *	len = end - beg;
 *
 * Therefore, the interpretation of bam_calend is not exactly as reported in 
 * 	http://samtools.sourceforge.net/samtools/bam/Functions/Functions.html#//apple_ref/c/func/bam_calend
 *
 * This is why the loop in mUpdateCoverageForAlignment() below goes
 *	(i=start; i<end; i++)
 */

void mInitCoverage() {
	int     n_targets  = global->header->n_targets;
	float **f_coverage = (float**) mMalloc(n_targets*sizeof(float*));
	int    *covered    = (int*) mCalloc(n_targets, sizeof(int));

	global->f_coverage = f_coverage;
	global->covered    = covered;
}

void mFreeCoverage() {
	int i;
	int      n_targets  = global->header->n_targets;
	int     *covered    = global->covered;
	float  **f_coverage = global->f_coverage;
	for (i=0; i<n_targets; i++) {
		if (covered[i]) {
			mFree(f_coverage[i]);
		}
	}
	mFree(f_coverage);
	mFree(covered);
}

void mUpdateCoverageForAlignment(bam1_t *bam, float coverage) {
	int i;
	int tid   = bam->core.tid;
	int start = bam->core.pos;
	int end   = bam_calend(&bam->core, bam1_cigar(bam));

	int    *covered    = global->covered;
	float **f_coverage = global->f_coverage;
	float  *this_coverage;

#ifdef DEBUG
	fprintf(stderr, "Tagging target %d (%s) with %f: from %d to %d\n", tid, header->target_name[tid], coverage, start, end);
	fprintf(stderr, "This tag was %d query bases long\n", bam_cigar2qlen(&bam->core, bam1_cigar(bam)));
#endif

	/* Why check before you write? Read takes less time than write, and write happens only n_targets times MAX */

	if (covered[tid] == 0) {
		int tlen        = global->header->target_len[tid];
		covered[tid]    = 1; 
		f_coverage[tid] = (float*) mCalloc(tlen, sizeof(float));
	}

	/* It's important to assign this_coverage after the previous conditional clause.
	 * It is possible that f_coverage[tid] is not created yet, in which case the previous clause will create it.
	 * So wait until it finishes allocating memory before using it here!
	 */

	this_coverage = f_coverage[tid];
	for (i=start; i<end; i++) {
		this_coverage[i] += coverage;
	}
}

void mEstimateCoverageOnPool(mBamPool *pool) {
	int i;
	bam1_t **elem = pool->elem;
	int size      = pool->size;
	int n         = size;

	/* This pool of n alignments for 1 read needs to be shared between all the n refsequences 1/n each */
	float cov = 1.0/n;
	for (i=0; i<size; i++) {
		mUpdateCoverageForAlignment(elem[i], cov);
	}
}

void mEstimateCoverageOnFile(samfile_t *input) {

	int       pool_limit = 64;
	bam1_t   *current;
	char     *prev_read = (char*) mCalloc(128, sizeof(char));
	mBamPool *pool      = (mBamPool*) mCalloc(1, sizeof(mBamPool));

	uint32_t  mutual_pairs = (BAM_FREAD1 | BAM_FREAD2);
	uint32_t  prev_flag    = 0;

	/* init pool */

	mInitBamPool(pool, pool_limit);
	current = pool_current(pool);

	prev_read[0] = '\0';
	while (samread(input, current) >= 0) {
		if ( (prev_read[0] != '\0') && 
		     ((strcmp(bam1_qname(current), prev_read) != 0) || 
		      (((current->core.flag | prev_flag) & mutual_pairs) == mutual_pairs)
		     )
		   ) {
			mEstimateCoverageOnPool(pool);
			mReOriginateBamPool(pool);
			current = pool_current(pool);
		}
		prev_flag = current->core.flag;
		strncpy(prev_read, bam1_qname(current), 127);
		current   = mAdvanceBamPool(pool);
	}
	mEstimateCoverageOnPool(pool);
	mFreeBamPool(pool);
	mFree(pool);
	mFree(prev_read);
}

/* Write compressed coverage data using pipe: child is the compressor */

void mWriteCoverageToStream(FILE* outstream, int skip_uncovered, int wordsize) {
	int      tid;
	int      n_targets  = global->header->n_targets;
	int     *covered    = global->covered;
	float  **f_coverage = global->f_coverage;
	char    *s1         = (char*) mMalloc(sizeof(char) * 20);
	char    *s2;

	for (tid=0; tid<n_targets; tid++) {
		int32_t  i;
		float    sum  = 0.0;
		int32_t  tlen = global->header->target_len[tid];
		int      integer_length;

		/* If this target is not covered, deal with it */

		if (!covered[tid]) {
			/*****
			 * global->f_coverage[tid] will only be initialized if covered[tid]==1.
			 * We will print 0's for the tlen if skip_uncovered==0, or do nothing if skip_uncovered==1.
			 */
			if (!skip_uncovered) {
				for (i=0; i<tlen-1; i++) {
					if ((i+1)%wordsize == 0) {
						fprintf(outstream, "0\n");
					} else {
						fprintf(outstream, "0 ");
					}
				}
				fprintf(outstream, "0\n");
			}
			continue;
		}

		/* Write to the given output stream */
		/* We use 2 decimal precision when the number is not an integer, and thus integer_length+2 for %g */

		fprintf(outstream, ">%s\n", global->header->target_name[tid]);
		for (i=0; i<tlen-1; i++) {
			float val = f_coverage[tid][i];
			sum += val;
			if (val > 0) {
				sprintf(s1, "%f", val);
				s2 = strchr(s1, '.');
				integer_length = s2 - s1;
				if ((i+1)%wordsize == 0) {
					fprintf(outstream, "%.*g\n", integer_length+2, val);
				} else {
					fprintf(outstream, "%.*g ",  integer_length+2, val);
				}
			} else {
				if ((i+1)%wordsize == 0) {
					fprintf(outstream, "0\n");
				} else {
					fprintf(outstream, "0 ");
				}
			}
		}
		fprintf(outstream, "%.4f\n", f_coverage[tid][tlen-1]);
		/*fprintf(stderr, "%s\t%f\n", global->header->target_name[tid], sum/tlen);*/
	}
	mFree(s1);
}

#define subprogram "coverage"

int msam_coverage_main(int argc, char* argv[]) {

	/* common variables */

	char             outfile[256];
	const char      *infile;
	char            *headerfile = NULL;

	samfile_t       *input  = NULL;

	const char      *inmode;

	/* argtable related */

	void           **argtable;
	struct arg_lit  *arg_samin;
	struct arg_file *arg_samfile;
	struct arg_lit  *arg_help;
	struct arg_end  *end;

	struct arg_lit  *arg_skip_uncovered;
	struct arg_lit  *arg_gzip;
/*
	struct arg_lit  *arg_share_multi_hit;
*/
	struct arg_int  *arg_wordsize;
	struct arg_str  *arg_out;
	int              set_argcount = 0;

	/* Program related */

	int              wordsize = 17;
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
	arg_out             = arg_str1("o",  NULL,      "<file>","name of output file (required)");
/*
	arg_share_multi_hit = arg_lit0("s", "sharemultihit",     "share the read that hits multiple locations equally across locations (default: true)");
*/
	arg_skip_uncovered  = arg_lit0("x", "skipuncovered",     "do not report coverage for sequences without aligned reads (default: false)");
	arg_wordsize        = arg_int0("w", "wordsize", NULL,    "number of words (coverage values) per line (default: 17)");
	arg_gzip            = arg_lit0("z", "gzip",              "compress output file using gzip (default: false)\n"
                                                                 "\n"
                                                                 "Description:\n"
                                                                 "------------\n"
                                                                 "\n"
                                                                 "Produces per-position sequence coverage information for all reference sequences\n"
                                                                 "in the BAM file. Output is similar to old-style quality files from the Sanger \n"
                                                                 "sequencing era, with a fasta-style header followed by lines of space-delimited \n"
                                                                 "numbers.\n"
                                                                 "\n"
                                                                 "For large datasets, option '-x' comes in handy when only a small fraction of \n"
                                                                 "reference sequences are covered.\n"
                                                                 "\n"
                                                                 "If using '-z', output file does NOT automatically get '.gz' extension. This is \n"
                                                                 "up to the user to specify the correct full output file name."
                                                                 );
	end    = arg_end(20); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */

	argtable = (void**) mCalloc(8, sizeof(void*));

	/* Common args */
	set_argcount = 0;
	argtable[set_argcount++] = arg_samin;
	argtable[set_argcount++] = arg_samfile;
	argtable[set_argcount++] = arg_help;

	/* Specific args */
	argtable[set_argcount++] = arg_out;
/* This is default, so no need to add this */
/*
	argtable[set_argcount++] = arg_share_multi_hit;
*/
	argtable[set_argcount++] = arg_skip_uncovered;
	argtable[set_argcount++] = arg_wordsize;
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

	if (arg_out->count != 1) {
		fprintf(stdout, "requires -o\n");
		mPrintHelp(subprogram, argtable);
		mQuit("");
	}

	/* Set output stream */
	/* I set this early enough, because gzip output uses a fork/pipe where a child exits. */
	/* Memory allocated before the fork is all reported as lost. */
	/* To minimize the reported loss, I now open the output stream as early as possible */

	if (arg_gzip->count > 0) {
		gzip = 1;
	}
	strcpy(outfile, arg_out->sval[0]);
	output = mInitOutputStream(outfile, gzip);

	/* Setup input operations */

	inmode = M_INPUT_MODE(arg_samin);
	infile = arg_samfile->filename[0];
	input = mOpenSamFile(infile, inmode, headerfile);
	global->header = input->header;

	if (arg_wordsize->count > 0) {
		wordsize = arg_wordsize->ival[0];
	}

	/* Specific operations */

	mInitCoverage();
	mEstimateCoverageOnFile(input);

	/* Write output */

	mWriteCoverageToStream(output, arg_skip_uncovered->count > 0, wordsize);
	mFreeOutputStream(output, gzip);

	/* Wind-up operations */

	mFreeCoverage();
	samclose(input);
	arg_freetable(argtable, 8);
	mFree(argtable);

	mFreeGlobal();

	return 0;
}
