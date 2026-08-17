#include "msam.h"
#include "mBamVector.h"

void mFilterFileWrapper(mSamFile *input, mSamFile *output, int uniqbesthit_only, int besthit_only, int rescore, int invert, int keep_unmapped);
void mFilterFile(mSamFile *input, mSamFile *output, int (*filter)(mAlignmentSummary*), void (*writer)(mSamFile*, mBamPool*), int rescore, int invert, int keep_unmapped);
void mFilterFileLite(mSamFile *input, mSamFile *output, void (*writer)(mSamFile*, mBamPool*));
void mWriteBestHitBamPool(mSamFile *stream, mBamPool *pool);
void mWriteUniqueBestHitBamPool(mSamFile *stream, mBamPool *pool);

/* Filters return 1 when the condition fails. This makes it easy to call
 * each filter and if one of them fails, you already quit.
 */


/****
 * PPT is a special case, since it takes both positive and negative values.
 * It is implemented as follows using edit distance:
 *
 * The following trick is to make sure that
 *  (1) if ppt is like  950, then you need above 950
 *  (2) if ppt is like -950, then you need below 950
 *
	dist = 1000*(alignment->length-alignment->edit);
	if (PPT < 0) dist = -dist;
	if (dist < alignment->length*PPT) continue;

 * This has now been converted into three #defined statements below.
 * Beware that they use three global variables that are not #define arguments:
 *  MIN_LENGTH, MAX_CLIP and PPT.
 ****/

#define _FILTER_L(a)     (a->length < global->MIN_LENGTH)
#define _FILTER_Z(a)     (100*a->query_clip > global->MAX_CLIP*a->query_length)
#define _FILTER_P_POS(a) (1000*(a->length - a->edit) < a->length*global->PPT)
#define _FILTER_P_NEG(a) (1000*(a->edit - a->length) < a->length*global->PPT)
#define _FILTER_P(a)     ((global->PPT<0)?_FILTER_P_NEG(a):_FILTER_P_POS(a))

static int filter_l(mAlignmentSummary* a) {
	return _FILTER_L(a);
}

static int filter_p(mAlignmentSummary* a) {
	return _FILTER_P(a);
}

static int filter_z(mAlignmentSummary* a) {
	return _FILTER_Z(a);
}

static int filter_lp(mAlignmentSummary* a) {
	return _FILTER_L(a) || _FILTER_P(a);
}

static int filter_lz(mAlignmentSummary* a) {
	return _FILTER_L(a) || _FILTER_Z(a);
}

static int filter_pz(mAlignmentSummary* a) {
	return _FILTER_P(a) || _FILTER_Z(a);
}

static int filter_lpz(mAlignmentSummary* a) {
	return _FILTER_L(a) || _FILTER_P(a) || _FILTER_Z(a);
}

void mFilterFileWrapper(mSamFile *input, mSamFile *output, int uniqbesthit_only, int besthit_only, int rescore, int invert, int keep_unmapped) {

	/* filter function */

	/* -l, -p, -z get 001, 010, 100 flags, so combination are ordered as follows. */

	int filter_choice = 0;
	int (*filter)(mAlignmentSummary *alignment);
	int (*filters[8])(mAlignmentSummary *alignment) = {NULL, filter_l, filter_p, filter_lp, filter_z, filter_lz, filter_pz, filter_lpz};

	/* writer function */
	void (*writer)(mSamFile*, mBamPool*);

	/* Assign filter function */
	if (global->MIN_LENGTH > 0) filter_choice |= 1; /* 001 for l */
	if (global->PPT       != 0) filter_choice |= 2; /* 010 for p */
	if (global->MAX_CLIP < 100) filter_choice |= 4; /* 100 for z */
	if (filter_choice == 0 && besthit_only == 0 && uniqbesthit_only == 0) {
		mDie("'filter' command requires atleast one of --ppt, -l, -p, -z, --besthit or --uniqhit");
	}
	filter = filters[filter_choice];

	/* Assign writer function */
	if (uniqbesthit_only)
		writer = mWriteUniqueBestHitBamPool;
	else if (besthit_only)
		writer = mWriteBestHitBamPool;
	else
		writer = mWriteBamPool;

	if (filter_choice == 0) { /* No alignment filtering, just besthit or uniqhit */
		mFilterFileLite(input, output, writer);
	} else { /* Need alignment filter and then perhaps besthit/uniqhit */
		mFilterFile(input, output, filter, writer, rescore, invert, keep_unmapped);
	}
}

void mFilterFile(mSamFile *input, mSamFile *output, int (*filter)(mAlignmentSummary*), void (*writer)(mSamFile*, mBamPool*), int rescore, int invert, int keep_unmapped) {

	/* parameters for rescoring alignment score */
	int       hit=1, miss=-1;

	int       pool_limit = 64;
	bam1_t   *current;
	mBamPool *pool = (mBamPool*) mCalloc(1, sizeof(mBamPool));
	char      prev_read[BAM_MAX_QNAME_LEN + 1];

	/* features to filter on */
	uint8_t  *nm;
	mAlignmentSummary* alignment = (mAlignmentSummary*) mMalloc(sizeof(mAlignmentSummary));

	/* init pool */

	mInitBamPool(pool, pool_limit);
	current = pool_current(pool);
	prev_read[0] = '\0';
	while (mSamRead(input, current) >= 0) {
		if ( (prev_read[0] != '\0') &&
		     (strcmp(bam_get_qname(current), prev_read) != 0) ) {
			writer(output, pool);
			mReOriginateBamPool(pool);
			current = pool_current(pool);
		}
		/*****
		 * Ignore an unmapped read, unless your selection uses upper limits.
		 *
		 * Normally, filter uses lower-limits, e.g. minimum percent id, minimum alignment length, etc.
		 * In all these cases, there is no need to process "unmapped" reads, because they are below threshold.
		 *
		 * However, when we use -v, then we are converting limits to upper-limits.
		 * Here, unmapped reads technically satisfy the criteria, but usually they are useless in output files.
		 * Thus, the default behavior is to ignore unmapped reads and avoid processing them.
		 *
		 * There are some cases, e.g. writing out unaligned read names or fastq sequences, where they do become useful.
		 * In those cases, using --keep_unmapped will write out unmapped read entries from the BAM file.
		 *
		 * This is also not done during --besthit and --uniqhit modes, as they will be combined with lower-limits.
		 *****/
		if (current->core.flag & BAM_FUNMAP) {
			if (keep_unmapped != 0) {
				if (global->PPT >= 0 && invert == 1)
					current = mAdvanceBamPool(pool);
			}
			continue;
		}

		/***
		 * For filtering, we need the following: length, query_length, edit.
		 * For rescoring, we need the same as above.

			if (MD is present):
				calculate #hits and #misses from MD
			elsif (NM is present):
				#misses = NM; #hits = length - #misses;
		 */

		/* Get alignment stats */
		if (bam_aux_get(current, "MD")) {
			bam_get_summary(current, alignment);
		} else {
			nm = bam_aux_get(current, "NM");
			if (!nm) {
				mDie("Either NM or MD must be present in SAM/BAM input for 'filter' command. Type '%s filter -h' for details.", PROGRAM);
			}

			bam_cigar2details(&current->core, bam_get_cigar(current), &alignment->length, &alignment->query_length, &alignment->query_clip);
			alignment->edit  = bam_aux2i(nm);
		}

		/* Rescore if necessary */
		if (rescore) {
			int score = (alignment->length - alignment->edit)*hit + alignment->edit*miss;
			uint8_t *as = bam_aux_get(current, "AS");
			if (as) {
				bam_aux_del(current, as);
			}
			bam_aux_append(current, "AS", 'i', 4, (uint8_t*)&score);
		}

		strcpy(prev_read, bam_get_qname(current));

		/***
		 * Do I pass the filter? "filter(alignment) == 1" means failed.
		 * under invert=0, filter=0,1 means advance=1,0 resp.
		 * under invert=1, filter=0,1 means advance=0,1 resp.
		 * The combination that covers this is: advance=(invert==filter)
		 */
		if (filter(alignment) == invert)
			current = mAdvanceBamPool(pool);
	}
	writer(output, pool);
	mFreeBamPool(pool);
	mFree(pool);
	mFree(alignment);
}

void mFilterFileLite(mSamFile *input, mSamFile *output, void (*writer)(mSamFile*, mBamPool*)) {

	/* parameters for rescoring alignment score */

	int       pool_limit = 64;
	bam1_t   *current;
	mBamPool *pool = (mBamPool*) mCalloc(1, sizeof(mBamPool));
	char      prev_read[BAM_MAX_QNAME_LEN + 1];

	/* features to filter on */

	/* init pool */

	mInitBamPool(pool, pool_limit);
	current = pool_current(pool);
	prev_read[0] = '\0';
	while (mSamRead(input, current) >= 0) {
		bam1_core_t  core = current->core;
		if ( (prev_read[0] != '\0') &&
		     (strcmp(bam_get_qname(current), prev_read) != 0) ) {
			writer(output, pool);
			mReOriginateBamPool(pool);
			current = pool_current(pool);
		}
		strcpy(prev_read, bam_get_qname(current));

		/* Ignore an unmapped read */

		if (core.flag & BAM_FUNMAP)
			continue;

		current = mAdvanceBamPool(pool);
	}
	writer(output, pool);
	mFreeBamPool(pool);
	mFree(pool);
}

static int mBamMatchesMate(const bam1_t *b, uint32_t mate_flag) {
	return (b->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == mate_flag;
}

static int mBamPoolIsPaired(mBamPool *pool) {
	int i;

	for (i=0; i<pool->size; i++) {
		if (pool->elem[i]->core.flag & (BAM_FREAD1 | BAM_FREAD2))
			return 1;
	}
	return 0;
}

static void mWriteBestHitBamPoolByMate(mSamFile *stream, mBamPool *pool, uint32_t mate_flag, int unique_only) {
	int        i;
	int        best_count = 0;
	int32_t    best_score = INT32_MIN;
	bam1_t   **elem       = pool->elem;

	for (i=0; i<pool->size; i++) {
		uint8_t *as;
		int32_t  score;

		if (!mBamMatchesMate(elem[i], mate_flag))
			continue;

		as = bam_aux_get(elem[i], "AS");
		if (!as)
			mDie("Required field AS not found in SAM/BAM input. Type '%s -h' for details.", PROGRAM);

		score = bam_aux2i(as);
		if (score > best_score) {
			best_score = score;
			best_count = 1;
		} else if (score == best_score) {
			best_count++;
		}
	}

	if (best_count == 0 || (unique_only && best_count != 1))
		return;

	for (i=0; i<pool->size; i++) {
		uint8_t *as;

		if (!mBamMatchesMate(elem[i], mate_flag))
			continue;

		as = bam_aux_get(elem[i], "AS");
		if (bam_aux2i(as) == best_score)
			mSamWrite(stream, elem[i]);
	}
}

void mWriteBestHitBamPool(mSamFile *stream, mBamPool *pool) {
	if (mBamPoolIsPaired(pool)) {
		mWriteBestHitBamPoolByMate(stream, pool, BAM_FREAD1, 0);
		mWriteBestHitBamPoolByMate(stream, pool, BAM_FREAD2, 0);
	} else {
		mWriteBestHitBamPoolByMate(stream, pool, 0, 0);
	}
}

void mWriteUniqueBestHitBamPool(mSamFile *stream, mBamPool *pool) {
	if (mBamPoolIsPaired(pool)) {
		mWriteBestHitBamPoolByMate(stream, pool, BAM_FREAD1, 1);
		mWriteBestHitBamPoolByMate(stream, pool, BAM_FREAD2, 1);
	} else {
		mWriteBestHitBamPoolByMate(stream, pool, 0, 1);
	}
}

#define subprogram "filter"

int msam_filter_main(int argc, char* argv[]) {

	/* common variables */

	const char      *infile = NULL;

	mSamFile        *input  = NULL;
	mSamFile        *output = NULL;
	sam_hdr_t       *output_header = NULL;
	mQNameCheckResult qname_check;

	const char      *inmode;
	char             outmode[6];

	/* argtable related */
	void           **argtable;
	struct arg_lit  *arg_samin;
	struct arg_lit  *arg_bamout;
	struct arg_lit  *arg_uncompressed;
	struct arg_lit  *arg_write_header;
	struct arg_file *arg_samfile;
	struct arg_lit  *arg_help;
	struct arg_end  *end;
	struct arg_lit  *arg_rescore;
	struct arg_int  *arg_minlength;
	struct arg_int  *arg_minpercentid;
	struct arg_int  *arg_minppt;
	struct arg_int  *arg_minqfrac;
	struct arg_lit  *arg_invertfilter;
	struct arg_lit  *arg_keepunmapped;
	struct arg_lit  *arg_besthitonly;
	struct arg_lit  *arg_uniqbesthitonly;

	int              set_argcount = 0;

	mInitGlobal();
	global->multiple_input = 0;
	arg_bamout       = arg_lit0("b", NULL, "output BAM (default: false)");
	arg_uncompressed = arg_lit0("u", NULL, "uncompressed BAM output (force -b) (default: false)");
	arg_write_header = arg_lit0("h", NULL, "print header for the SAM output (default: false)");
	arg_samin        = arg_lit0("S",  NULL,                     "input is SAM (default: false)");
	arg_samfile      = arg_filen(NULL, NULL, "<bamfile>", 1, 1, "input SAM/BAM file");
	arg_help         = arg_lit0(NULL, "help",                   "print this help and exit\n\n"
								    "Specific options:\n"
								    "-----------------\n");
	/* Specific args */
	arg_minlength       = arg_int0("l",  NULL,     NULL,           "min. length of alignment (default: 0)");
	arg_minpercentid    = arg_int0("p",  NULL,     NULL,           "min. sequence identity of alignment, in percentage, integer between 0 and 100; requires NM field to be present (default: 0)");
	arg_minppt          = arg_int0(NULL, "ppt",    NULL,           "min/max sequence identity of alignment, in parts per thousand, integer between -1000 and 1000; requires NM field to be present (default: 0)\n"
									"                            NOTE:\n"
									"                            -----\n"
									"                                  When using --ppt, +ve values mean minimum ppt and -ve values mean maximum ppt.\n"
									"                                  E.g., '--ppt 950' will report alignments with ppt>950,\n"
									"                                  and '--ppt -950' will report alignments with ppt<=950.");
	arg_minqfrac        = arg_int0("z",  NULL,     NULL,           "min. percent of the query that must be aligned, between 0 and 100 (default: 0)");
	arg_keepunmapped    = arg_lit0("k",  "keep_unmapped",          "report unmapped reads, when filtering using upper-limit thresholds (default: false)");
	arg_invertfilter    = arg_lit0("v",  "invert",                 "invert the effect of the filter (default: false)\n"
									"                            CAUTION:\n"
									"                            --------\n"
									"                                  When using --invert or -v, be precise in what needs to be inverted.\n"
									"                                  Adding -v gives you the complement of what you get without -v.\n"
									"                                  Sometimes, this might be counter-intuitive.\n"
									"                                  E.g., '-l 65 -p 95' will report alignments that are (>65bp AND >95%).\n"
									"                                        '-l 65 -p 95 -v' will not report (<65bp AND <95%), as one might think.\n"
									"                                        '-l 65 -p 95 -v' will report NOT(>65bp AND >95%) which is (<65bp OR <95%).\n"
									"                                        Notice the 'OR' in the final logical operation. This means that\n"
									"                                        an alignment that fails one condition will still be reported if it\n"
									"                                        satisfies the other condition.\n"
									"                                        If you only wanted alignments that are below 95%, then specify '-p 95 -v'");
	arg_rescore         = arg_lit0(NULL, "rescore",                "rescore alignments using MD or NM fields, in that order (default: false)\n\n"
									"Special filters:\n"
									"----------------\n\n"
									"The following special filters cannot be combined with -v, but require:\n"
									"  (1) the alignments to be sorted by name,\n"
									"  (2) AS field (alignment score) to be present.\n"
									"You can usually achieve sorting by:\n"
									"  samtools sort -n -T tmp.sort input.bam  | "PROGRAM" -m filter --besthit -\n"
									"If AS is missing, you can rescore alignments by:\n"
									"  samtools sort -n -T tmp.sort input.bam | "PROGRAM" -m filter --rescore --besthit -\n");
	arg_besthitonly     = arg_lit0(NULL, "besthit",                "keep all highest scoring hit(s) per read (default: false)");
	arg_uniqbesthitonly = arg_lit0(NULL, "uniqhit",                "keep only one highest scoring hit per read, only if it is unique (default: false)");
	end                 = arg_end(16); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */
	argtable = (void**) mCalloc(16, sizeof(void*));

	/* Common args */
	set_argcount = 0;
	argtable[set_argcount++] = arg_bamout;
	argtable[set_argcount++] = arg_uncompressed;
	argtable[set_argcount++] = arg_write_header;
	argtable[set_argcount++] = arg_samin;
	argtable[set_argcount++] = arg_samfile;
	argtable[set_argcount++] = arg_help;

	/* Specific args */
	argtable[set_argcount++] = arg_minlength;
	argtable[set_argcount++] = arg_minpercentid;
	argtable[set_argcount++] = arg_minppt;
	argtable[set_argcount++] = arg_minqfrac;
	argtable[set_argcount++] = arg_keepunmapped;
	argtable[set_argcount++] = arg_invertfilter;
	argtable[set_argcount++] = arg_rescore;
	argtable[set_argcount++] = arg_besthitonly;
	argtable[set_argcount++] = arg_uniqbesthitonly;
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
		if (arg_help->count > 0 || argc < 2) {
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
	if (arg_invertfilter->count > 0 && (arg_besthitonly->count > 0 || arg_uniqbesthitonly->count > 0)) {
		fprintf(stdout, "--invert cannot be combined with --besthit or --uniqhit\n");
		mPrintHelp(subprogram, argtable);
		mQuit("");
	} else if (arg_besthitonly->count > 0 && arg_uniqbesthitonly->count > 0) {
		fprintf(stdout, "--besthit cannot be combined with --uniqhit\n");
		mPrintHelp(subprogram, argtable);
		mQuit("");
	} else if (arg_minlength->count       == 0 &&
	           arg_minpercentid->count    == 0 &&
	           arg_minppt->count          == 0 &&
	           arg_uniqbesthitonly->count == 0 &&
	           arg_besthitonly->count     == 0 &&
	           arg_minqfrac->count        == 0) {
		fprintf(stdout, "--mode filter needs -l, -p, --ppt, -z, --besthit or --uniqhit\n");
		mPrintHelp(subprogram, argtable);
		mQuit("");
	} else {
		int32_t percent_id  = 0;
		if (arg_minpercentid->count > 0) {
			percent_id  = arg_minpercentid->ival[0];
		}
		global->PPT = 10*percent_id;
		if (arg_minppt->count > 0) {
			global->PPT = arg_minppt->ival[0];
		}

		if (global->PPT < -1000 || global->PPT > 1000) {
			fprintf(stdout, "-p or --ppt must be in the range [-100,100] or [-1000,1000], respectively\n");
			mPrintHelp(subprogram, argtable);
			mQuit("");
		}

		global->MAX_CLIP     = 100;
		if (arg_minqfrac->count > 0) {
			global->MAX_CLIP = 100 - arg_minqfrac->ival[0];
		}
		if (global->MAX_CLIP < 0 || global->MAX_CLIP > 100) {
			fprintf(stdout, "-z must be in the range [-100,100]\n");
			mPrintHelp(subprogram, argtable);
			mQuit("");
		}

		global->MIN_LENGTH   = 0;
		if (arg_minlength->count > 0) {
			global->MIN_LENGTH   = (int32_t) arg_minlength->ival[0];
		}
	}

	/* Set input mode */
	inmode = M_INPUT_MODE(arg_samin);

	/* Set output mode */
	strcpy(outmode, "w");
	if (arg_uncompressed->count > 0)
		strcat(outmode, "bu");
	else if (arg_bamout->count > 0)
		strcat(outmode, "b");
	else if (arg_write_header->count > 0)
		strcat(outmode, "h");

	/* General operations */

	infile = arg_samfile->filename[0];
	input = mOpenSamInput(infile, inmode);
	global->header = input->header;

	if (arg_besthitonly->count > 0 || arg_uniqbesthitonly->count > 0) {
		qname_check = mSamCheckQNameGrouping(input);
	} else {
		qname_check = mQNameCheckNotRequired();
	}

	/*
	 * Annotate only the output-header copy.  The input header remains
	 * untouched.
	 */
	output_header = sam_hdr_dup(global->header);
	if (output_header == NULL)
		mDie("Cannot duplicate SAM header for output provenance");
	mAddSamProvenance(output_header, argc, argv, &qname_check);
	output = mOpenSamOutput("-", outmode, output_header);
	sam_hdr_destroy(output_header);

	/* Specific operations */

	mFilterFileWrapper(input, output, arg_uniqbesthitonly->count > 0, arg_besthitonly->count > 0, arg_rescore->count > 0, arg_invertfilter->count > 0, arg_keepunmapped->count > 0);

	/* Wind-up operations */

	mSamClose(input);
	mSamClose(output);
	arg_freetable(argtable, set_argcount);
	mFree(argtable);
	mFreeGlobal();
	return 0;
}
