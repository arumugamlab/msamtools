#include "msam.h"

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
	int          n_targets  = global->header->n_targets;
	coverage_t **coverage = (coverage_t**) mMalloc(n_targets*sizeof(coverage_t*));
	int         *covered    = (int*) mCalloc(n_targets, sizeof(int));

	global->coverage = coverage;
	global->covered    = covered;
}

void mFreeCoverage() {
	int i;
	int          n_targets = global->header->n_targets;
	int         *covered   = global->covered;
	coverage_t **coverage  = global->coverage;
	for (i=0; i<n_targets; i++) {
		if (covered[i]) {
			mFree(coverage[i]);
		}
	}
	mFree(coverage);
	mFree(covered);
}

void mUpdateCoverageForAlignment(bam1_t *bam, coverage_t in_coverage) {
	int i;
	int tid   = bam->core.tid;
	int start = bam->core.pos;
	int end   = bam_calend(&bam->core, bam1_cigar(bam));

	int         *covered  = global->covered;
	coverage_t **coverage = global->coverage;
	coverage_t  *this_coverage;

#ifdef DEBUG
	fprintf(stderr, "Tagging target %d (%s) with %f: from %d to %d\n", tid, header->target_name[tid], in_coverage, start, end);
	fprintf(stderr, "This tag was %d query bases long\n", bam_cigar2qlen(&bam->core, bam1_cigar(bam)));
#endif

	/* Why check before you write? Read takes less time than write, and write happens only n_targets times MAX */

	if (covered[tid] == 0) {
		int tlen      = global->header->target_len[tid];
		covered[tid]  = 1; 
		coverage[tid] = (coverage_t*) mCalloc(tlen, sizeof(coverage_t));
	}

	/* It's important to assign this_coverage after the previous conditional clause.
	 * It is possible that coverage[tid] is not created yet, in which case the previous clause will create it.
	 * So wait until it finishes allocating memory before using it here!
	 */

	this_coverage = coverage[tid];
	for (i=start; i<end; i++) {
		this_coverage[i] += in_coverage;
	}
}

void mEstimateCoverageOnPool(mBamPool *pool) {
	int i;
	bam1_t **elem = pool->elem;
	int size      = pool->size;

	/* This pool of n alignments for 1 read needs to be shared between all the n refsequences 1/n each */
	/* I believe we don't need to share. Repeated regions are covered multiple times, so let them be */
	/* int n         = size; */
	/* coverage_t cov = 1.0/n; */

	/* This pool of n alignments for 1 read get 1 each */
	coverage_t cov = 1;
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

void mWriteCoverageToStream(gzFile stream, int skip_uncovered, int wordsize) {
	int          tid;
	int          n_targets  = global->header->n_targets;
	int         *covered    = global->covered;
	coverage_t **coverage = global->coverage;

	for (tid=0; tid<n_targets; tid++) {
		int32_t    i;
		/* coverage_t sum  = 0; */
		int32_t    tlen = global->header->target_len[tid];

		/* If this target is not covered, deal with it */

		if (!covered[tid]) {
			/*****
			 * global->coverage[tid] will only be initialized if covered[tid]==1.
			 * We will print 0's for the tlen if skip_uncovered==0, or do nothing if skip_uncovered==1.
			 */
			if (!skip_uncovered) {
				gzprintf(stream, ">%s\n", global->header->target_name[tid]);
				for (i=0; i<tlen-1; i++) {
					if ((i+1)%wordsize == 0) {
						gzprintf(stream, "0\n");
					} else {
						gzprintf(stream, "0 ");
					}
				}
				gzprintf(stream, "0\n");
			}
			continue;
		}

		/* Write to the given output stream */

		gzprintf(stream, ">%s\n", global->header->target_name[tid]);
		for (i=0; i<tlen-1; i++) {
			coverage_t val = coverage[tid][i];
			/*sum += val; */
			if ((i+1)%wordsize == 0) {
				gzprintf(stream, COVERAGE_T_FORMAT"\n", val);
			} else {
				gzprintf(stream, COVERAGE_T_FORMAT" ", val);
			}
		}
		gzprintf(stream, COVERAGE_T_FORMAT"\n", coverage[tid][tlen-1]);
		/*fprintf(stderr, "%s\t%f\n", global->header->target_name[tid], sum/tlen);*/
	}
}

void mWriteCoverageSummaryToStream(gzFile stream, int skip_uncovered) {
	int          tid;
	int          n_targets  = global->header->n_targets;
	int         *covered    = global->covered;
	coverage_t **coverage = global->coverage;

	for (tid=0; tid<n_targets; tid++) {
		int32_t i;
		int32_t touched = 0;
		int64_t sum     = 0;
		int32_t tlen    = global->header->target_len[tid];

		/* If this target is not covered, deal with it */

		if (!covered[tid]) {
			if (!skip_uncovered) {
				gzprintf(stream, "%s\t%d\t%d\n", global->header->target_name[tid], 0, 0);
			}
			continue;
		}

		/* Write to the given output stream */

		for (i=0; i<tlen-1; i++) {
			coverage_t val = coverage[tid][i];
			touched += (val != 0);
			sum += val;
		}
		gzprintf(stream, "%s\t%.8f\t%.2f\n", global->header->target_name[tid], 1.0*touched/tlen, 1.0*sum/tlen);
	}
}

#define subprogram "coverage"

int msam_coverage_main(int argc, char* argv[]) {

	/* common variables */

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
	struct arg_lit  *arg_summary_only;
	struct arg_lit  *arg_gzip;
/*
	struct arg_lit  *arg_share_multi_hit;
*/
	struct arg_int  *arg_wordsize;
	struct arg_str  *arg_out;
	int              set_argcount = 0;

	/* Program related */

	int              wordsize = 17;
	gzFile           output = NULL;

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
	arg_summary_only    = arg_lit0(NULL, "summary",          "do not report per-position coverage but report fraction of sequence covered (default: false)");
	arg_skip_uncovered  = arg_lit0("x", "skipuncovered",     "do not report coverage for sequences without aligned reads (default: false)");
	arg_wordsize        = arg_int0("w", "wordsize", NULL,    "number of words (coverage values) per line (default: 17)");
	arg_gzip            = arg_lit0("z", "gzip",              "compress output file using gzip (default: true)\n"
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

	argtable = (void**) mCalloc(9, sizeof(void*));

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
	argtable[set_argcount++] = arg_summary_only;
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
			exit(EXIT_SUCCESS);
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

	/* Open output file */
	if (strcmp(arg_out->sval[0], "-") == 0) { /* If '-' was given as output, redirect to stdout */
		output = gzdopen(fileno(stdout), "wb");
	} else {
		output = gzopen(arg_out->sval[0], "wb");
	}


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

	if (arg_summary_only->count > 0) {
		mWriteCoverageSummaryToStream(output, arg_skip_uncovered->count > 0);
	} else {
		mWriteCoverageToStream(output, arg_skip_uncovered->count > 0, wordsize);
	}

	/* Close output */
	gzclose(output);


	/* Wind-up operations */

	mFreeCoverage();
	samclose(input);
	arg_freetable(argtable, set_argcount);
	mFree(argtable);

	mFreeGlobal();

	return 0;
}
