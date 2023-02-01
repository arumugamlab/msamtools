#include "msam.h"

#define MAP_STATS (0)
#define UNMAP_STATS (1)
#define EDIT_STATS (2)
#define SCOR_STATS (3)

void mSummarizeAlignmentsHQ(samfile_t *input, FILE *output) {

	int     i;
	int     min_edit = 4096;
	char   *prev_read = (char*) mCalloc(128, sizeof(char));
	bam1_t *b    = bam_init1();
	long   *dist = (long*) mCalloc(500, sizeof(long));
	mAlignmentSummary* alignment = (mAlignmentSummary*) mMalloc(sizeof(mAlignmentSummary));
	char  **target_name = global->header->target_name;

	prev_read[0] = '\0';
	while (samread(input, b) >= 0) {

		int edit;
		bam1_core_t *core = &b->core;
		int tid           = core->tid;
		int start         = core->pos;
		int end           = bam_calend(core, bam1_cigar(b));

		if (prev_read[0] != '\0' && strcmp(bam1_qname(b), prev_read) != 0) {
			if (min_edit != 4096) {
				dist[min_edit]++;
				min_edit = 4096;
			}
		}
		bam_get_extended_summary(b, alignment);
		if (bam_highquality_differences_only(b, 20)) {
			edit = alignment->mismatch + alignment->gapopen + alignment->gapextend;
			fprintf(output, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n", bam1_qname(b), alignment->query_length, target_name[tid], start, end, alignment->match, edit);
			if (edit < min_edit) min_edit = edit;
		}
		strncpy(prev_read, bam1_qname(b), 127);
	}
	bam_destroy1(b);
	if (min_edit != 4096) dist[min_edit]++;
	for (i=0; i<500; i++) {
		if (dist[i] > 0) {
			fprintf(stderr, "%d\t%ld\n", i, dist[i]);
		}
	}
	mFree(prev_read);
}

void mSummarizeAlignments_v1(samfile_t *input, FILE *output) {

	int     i;
	int     edit = 0;
	char   *prev_read = (char*) mCalloc(128, sizeof(char));
	bam1_t *b    = bam_init1();
	long   *dist = (long*) mCalloc(500, sizeof(long));
	mAlignmentSummary* alignment = (mAlignmentSummary*) mMalloc(sizeof(mAlignmentSummary));
	char  **target_name = global->header->target_name;
	uint32_t *target_len = global->header->target_len;

	prev_read[0] = '\0';
	while (samread(input, b) >= 0) {

		bam1_core_t *core = &b->core;
		int tid           = core->tid;
		int start         = core->pos + 1;
		int end           = bam_calend(core, bam1_cigar(b));

		if (prev_read[0] != '\0' && strcmp(bam1_qname(b), prev_read) != 0) {
/* wrong uninitialized reference. but not planning to fix this bug as this function is deprecated */
			dist[edit]++;
		}

		/* Skip reads mapping to the edges of ref */
		if ((start < 150) || (target_len[tid] - end < 150)) continue;

		/* Skip non-primary multi-mappers */
		if (core->flag & BAM_FSECONDARY) continue;

		bam_get_extended_summary(b, alignment);
		edit = alignment->edit;
		fprintf(output, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n", bam1_qname(b), alignment->query_length, target_name[tid], start, end, alignment->match, edit);
		strncpy(prev_read, bam1_qname(b), 127);
	}
	bam_destroy1(b);
	dist[edit]++;
	for (i=0; i<500; i++) {
		if (dist[i] > 0) {
			fprintf(stderr, "%d\t%ld\n", i, dist[i]);
		}
	}
	mFree(prev_read);
}

void mSummarizeAlignments_v2(samfile_t *input, FILE *output) {

	int     edit;
	bam1_t *b    = bam_init1();
	mAlignmentSummary* alignment = (mAlignmentSummary*) mMalloc(sizeof(mAlignmentSummary));
	char  **target_name = global->header->target_name;
	uint32_t *target_len = global->header->target_len;

	while (samread(input, b) >= 0) {

		bam1_core_t *core = &b->core;
		int tid           = core->tid;
		int start         = core->pos + 1;
		int end           = bam_calend(core, bam1_cigar(b));

		/* Skip reads mapping to the edges of ref */
		if ((start < 150) || (target_len[tid] - end < 150)) continue;

		/* Skip non-primary multi-mappers */
		if (core->flag & BAM_FSECONDARY) continue;

		bam_get_extended_summary(b, alignment);
		edit = alignment->edit;
		fprintf(output, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n", bam1_qname(b), alignment->query_length, target_name[tid], start, end, alignment->match, edit);
	}
	bam_destroy1(b);
}

/***************************************************************
 * Summarize alignment statistics for PRIMARY ALIGNMENTS ONLY! *
 * For statistics, counting multiple mappers more than once    *
 * will skew the numbers. E.g., number of mapped bases will be *
 * incorrect if we count the same read twice. Same for number  *
 * of mismatches. Therefore, I will only count primary reads.  *
 *  Also, I count softclip as mismatch. So beware!             *
 ***************************************************************/

int mCountInserts(samfile_t *input) {

	bam1_t *b         = bam_init1();
	char   *prev_read = (char*) mCalloc(128, sizeof(char));
	int     count     = 0;

	while (samread(input, b) >= 0) {

		bam1_core_t *core = &b->core;

		/* Skip unmapped */
		if (core->flag & BAM_FUNMAP) continue;

		if (strcmp(bam1_qname(b), prev_read) != 0) {
			count++;
		}

		strncpy(prev_read, bam1_qname(b), 127);
	}
	bam_destroy1(b);
	mFree(prev_read);
	return count;
}

void mSummarizeAlignments(samfile_t *input, FILE *output, uint32_t edge_len) {

	bam1_t *b    = bam_init1();
	mAlignmentSummary* alignment = (mAlignmentSummary*) mMalloc(sizeof(mAlignmentSummary));
	char  **target_name = global->header->target_name;
	uint32_t *target_len = global->header->target_len;

	while (samread(input, b) >= 0) {

		bam1_core_t *core = &b->core;
		int tid           = core->tid;
		uint32_t start    = core->pos + 1; /* I know I am assigning int to uint, but bam.c does it too! */
		uint32_t end      = bam_calend(core, bam1_cigar(b));

		int32_t  glocal_len;

		/* Skip unmapped */
		if (core->flag & BAM_FUNMAP) continue;

		/* Skip non-primary multi-mappers */
		if (core->flag & BAM_FSECONDARY) continue;

		/* Skip reads mapping to the edges of ref */
		if ((start < edge_len) || (target_len[tid] - end < edge_len)) continue;

		bam_get_extended_summary(b, alignment);
		glocal_len = alignment->length + alignment->query_clip;
		fprintf(output, "%s\t%d\t%s\t%d\t%d\t%.1f\n", bam1_qname(b), alignment->query_length, target_name[tid], glocal_len, alignment->match, 100.0 - 100.0*alignment->edit/glocal_len);

	}
	bam_destroy1(b);
	mFree(alignment);
}

void mSummarizeAlignmentsStats(int stats_type, samfile_t *input, FILE *output, uint32_t edge_len) {

	int     i;
	int     idx;
	bam1_t *b    = bam_init1();
	long   *dist = (long*) mCalloc(M_BAM_MAX_READ_LENGTH+1, sizeof(long));
	mAlignmentSummary* alignment = (mAlignmentSummary*) mMalloc(sizeof(mAlignmentSummary));
	uint32_t *target_len = global->header->target_len;
	uint64_t mapped   = 0;
	uint64_t unmapped = 0;
	uint64_t edits    = 0;
	uint64_t score    = 0;
	uint32_t stats[5];

	while (samread(input, b) >= 0) {

		bam1_core_t *core = &b->core;
		int      tid      = core->tid;
		int32_t  start    = core->pos + 1;
		uint32_t end      = bam_calend(core, bam1_cigar(b));

		/* Skip unmapped */
		if (core->flag & BAM_FUNMAP) continue;

		/* Skip non-primary multi-mappers */
		if (core->flag & BAM_FSECONDARY) continue;

		/* Skip reads mapping to the edges of ref */
		if ((start < (int32_t) edge_len) || (target_len[tid] - end < edge_len)) continue;

		bam_get_extended_summary(b, alignment);

		stats[MAP_STATS]   = alignment->match;
		stats[UNMAP_STATS] = alignment->query_length - alignment->match;
		stats[EDIT_STATS]  = alignment->edit;
		stats[SCOR_STATS]  = alignment->match - alignment->edit;

		  mapped += stats[MAP_STATS];
		unmapped += stats[UNMAP_STATS];
		   edits += stats[EDIT_STATS];
		   score += stats[SCOR_STATS];

		idx = stats[stats_type];
		if (idx > M_BAM_MAX_READ_LENGTH) idx = M_BAM_MAX_READ_LENGTH;
		if (idx < 0)                     idx = 0;
		dist[idx]++;
	}
	bam_destroy1(b);
	for (i=0; i<M_BAM_MAX_READ_LENGTH; i++) {
		if (dist[i] > 0) {
			fprintf(output, "%d\t%ld\n", i, dist[i]);
		}
	}
	if (dist[M_BAM_MAX_READ_LENGTH] > 0) {
		fprintf(output, "%d+\t%ld\n", M_BAM_MAX_READ_LENGTH, dist[i]);
	}

	/*fprintf(stderr, "  Mapped: %10ld\nUnmapped: %10ld\n   Edits: %10ld\n", mapped, unmapped, edits);*/
	mFree(dist);
	mFree(alignment);
}

#define subprogram "summary"

int msam_summary_main(int argc, char* argv[]) {

	/* common variables */

	const char      *infile;
	char            *headerfile = NULL;

	samfile_t       *input  = NULL;

	const char      *inmode;

	/* For this file */

	uint32_t         edge;

	/* argtable related */

	void           **argtable;
	struct arg_lit  *arg_samin;
	struct arg_file *arg_samfile;
	struct arg_lit  *arg_help;
	struct arg_end  *end;

	struct arg_int  *arg_edge;
	struct arg_str  *arg_stats;
	struct arg_lit  *arg_count;
	int              set_argcount = 0;

	mInitGlobal();
	global->multiple_input = 0;

	arg_samin        = arg_lit0("S",  NULL,                     "input is SAM (default: false)");
	arg_samfile      = arg_filen(NULL, NULL, "<bamfile>", 1, 1, "input SAM/BAM file");
	arg_help         = arg_lit0(NULL, "help",                   "print this help and exit\n\n"
								    "Specific options:\n"
								    "-----------------\n");


	/* Specific args */
	arg_edge         = arg_int0("e",  "edge",  "<num>", "ignore alignment if reads map to <num> bases at the edge of target sequence (default: 0)");
	arg_count        = arg_lit0("c",  "count",          "count number of unique inserts in BAM file (default: false)");
	arg_stats        = arg_str0(NULL, "stats", "<string>", "{mapped|unmapped|edit|score} only report readcount distribution for specified stats, not read-level stats (default: none)\n\n"
                                                             ""
                                                             "Description\n"
                                                             "-----------\n"
                                                             "Prints summary of alignments in the given BAM/SAM file. By default, it prints\n"
                                                             "a summary line per alignment entry in the file. The summary is a tab-delimited\n"
                                                             "line with the following fields:\n"
                                                             "\tqname,aligned_qlen,target_name,glocal_align_len,matches,percent_identity\n"
                                                             "glocal_align_len includes the unaligned qlen mimicing a global alignment \n"
                                                             "in the query and local alignment in target, thus glocal.\n\n"
                                                             "With --stats option, summary is consolidated as distribution of read counts\n"
                                                             "for a given measure. \n"
                                                             "   --stats=mapped   - distribution for number of mapped query bases\n"
                                                             "   --stats=unmapped - distribution for number of unmapped query bases\n"
                                                             "   --stats=edit     - distribution for edit distances\n"
                                                             "   --stats=score    - distribution for score=match-edit"
                                                             );
	end = arg_end(7); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */

	argtable = (void**) mCalloc(7, sizeof(void*));

	/* Common args */
	set_argcount = 0;
	argtable[set_argcount++] = arg_samin;
	argtable[set_argcount++] = arg_samfile;
	argtable[set_argcount++] = arg_help;

	/* Specific args */
	argtable[set_argcount++] = arg_edge;
	argtable[set_argcount++] = arg_count;
	argtable[set_argcount++] = arg_stats;
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

	/* Set input/output modes */

	inmode = M_INPUT_MODE(arg_samin);

	/* General operations */

	infile = arg_samfile->filename[0];
	input = mOpenSamFile(infile, inmode, headerfile);
	global->header = input->header;

	/* Specific operations */

	/* Is there an edge to avoid? */

	edge = 0;
	if (arg_edge->count > 0)
		edge = (uint32_t) arg_edge->ival[0];

	if (arg_stats->count > 0) {
		int i;
		int stat_mode = -1;
		/* IMPORTANT: Make sure that the following list matched #define 's on top */
		const char* modes[] = {"mapped", "unmapped", "edit", "score"};
		int n_modes = sizeof(modes)/sizeof(const char*);
		for (i=0; i<n_modes; i++) {
			if (strcmp(arg_stats->sval[0], modes[i]) == 0) {
				stat_mode = i;
				break;
			}
		}
		if (stat_mode == -1) {
			mDie("Do not understand %s as mode", arg_stats->sval[0]);
		}
		mSummarizeAlignmentsStats(stat_mode, input, stdout, edge);
	} else if (arg_count->count > 0) {
		fprintf(stdout, "%d\n", mCountInserts(input));
	} else {
		mSummarizeAlignments(input, stdout, edge);
	}

	/* Wind-up operations */

	samclose(input);
	arg_freetable(argtable, set_argcount);
	mFree(argtable);

	mFreeGlobal();

	return 0;
}
