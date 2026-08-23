#include "msam.h"
#include "msam_version.h"

#define QNAME_GROUP_CHECK_RECORDS 10000
#define COORD_ORDER_CHECK_RECORDS 100000
#define COORD_ORDER_MIN_RECORDS   10000

typedef struct {
	size_t first_record;
	size_t last_record;
} mClosedQNameGroup;

/*
 * Handle msam_global global variable
 */
void mInitGlobal() {
	global = (msam_global*) mMalloc(sizeof(msam_global));
	global->MIN_LENGTH = 0;
	global->PPT = 0;
	global->MAX_CLIP = 0;
	global->multiple_input = 0;

	/* Init all pointers to NULL */
	global->header = NULL;
	global->fmap = NULL;
	global->coverage = NULL;
	global->seq_touched = NULL;
	global->ui_insert_count = NULL;
	global->d_insert_count = NULL;
	global->multi_mappers = NULL;
	global->multi_mapper_count = 0;
	global->uniq_mapper_count = 0;
	global->purged_insert_count = 0;
	global->ub_target_hit = NULL;
}

void mFreeGlobal() {
	if (global->fmap != NULL) {
		mFree(global->fmap);
	}
	mFree(global);
	global = NULL;
}

/*
 * Useful functions for printing information
 */

void mPrintHelp (const char *subprogram, void **argtable) {
	fprintf(stdout, "Usage:\n------\n\n%s %s", PROGRAM, subprogram);
	arg_print_syntax(stdout, argtable, "\n");
	fprintf(stdout, "\nGeneral options:\n"
				"----------------\n\n"
				"These options specify the input/output formats of BAM/SAM files\n"
				"(same meaning as in 'samtools view'):\n");
	arg_print_glossary(stdout, argtable, "  %-25s %s\n");
}

static char *mBuildCommandLine(int argc, char *argv[]) {
	char **full_argv;
	char  *command;
	int    i;

	full_argv = (char**) mMalloc((argc + 1) * sizeof(char*));
	full_argv[0] = (char*) PROGRAM;
	for (i=0; i<argc; i++) {
		full_argv[i+1] = argv[i];
	}

	command = stringify_argv(argc + 1, full_argv);
	mFree(full_argv);
	if (command == NULL)
		mDie("Cannot construct command line for provenance");

	return command;
}

mQNameCheckResult mQNameCheckNotRequired(void) {
	mQNameCheckResult result;

	result.status = MSAM_QNAME_NOT_REQUIRED;
	result.qname_records_checked = 0;
	result.input_records_checked = 0;
	result.mapped_records_checked = 0;
	return result;
}

static void mFormatQNameCheck(const mQNameCheckResult *result, char *buffer, size_t buffer_size) {
	switch (result->status) {
		case MSAM_QNAME_NOT_REQUIRED:
			snprintf(
				buffer, buffer_size,
				"QNAME grouping check: not required for this operation"
			);
			break;

		case MSAM_QNAME_HEADER_CONFIRMED:
			snprintf(
				buffer, buffer_size,
				"QNAME grouping check: confirmed by input header SO:queryname"
			);
			break;

		case MSAM_QNAME_SAMPLE_OK:
			if (result->input_records_checked < QNAME_GROUP_CHECK_RECORDS) {
				snprintf(
					buffer, buffer_size,
					"QNAME grouping check: no QNAME grouping violation detected "
					"in all %zu records",
					result->qname_records_checked
				);
			} else {
				snprintf(
					buffer, buffer_size,
					"QNAME grouping check: no QNAME grouping violation detected "
					"in first %zu records",
					result->qname_records_checked
				);
			}
			break;

		case MSAM_QNAME_SAMPLE_WARNING:
			snprintf(
				buffer, buffer_size,
				"QNAME grouping check: WARNING - no QNAME grouping violation "
				"detected in first %zu records; %zu mapped records among the "
				"first %zu input records were consistent with coordinate ordering",
				result->qname_records_checked,
				result->mapped_records_checked,
				result->input_records_checked
			);
			break;

		default:
			mDie("Unknown QNAME grouping check status");
	}
}

void mPrintProfileProvenanceGzip(gzFile output, int argc, char *argv[], const mQNameCheckResult *qname_check) {
	char *command = mBuildCommandLine(argc, argv);
	char  qname_message[1024];

	mFormatQNameCheck(qname_check, qname_message, sizeof(qname_message));

	gzprintf(output, "# msamtools version: %s\n", PACKAGE_VERSION);
	gzprintf(output, "# msamtools git commit: %s\n", MSAM_GIT_COMMIT);
	gzprintf(output, "# Command line: %s\n", command);
	gzprintf(output, "# %s\n", qname_message);

	free(command);
}

void mAddSamProvenance(sam_hdr_t *header, int argc, char *argv[], const mQNameCheckResult *qname_check) {
	char *command = mBuildCommandLine(argc, argv);
	char  description[1252];
	char  qname_message[1024];

	mFormatQNameCheck(qname_check, qname_message, sizeof(qname_message));
	snprintf(
		description,
		sizeof(description),
		"git=%s; %s",
		MSAM_GIT_COMMIT, qname_message
	);

	/*
	 * HTSlib's sam_hdr_add_pg() generates a unique ID if "msamtools"
	 * already exists and maintains existing @PG chains.
	 */
	if (sam_hdr_add_pg(
		header,
		PROGRAM,
		"PN", PROGRAM,
		"VN", PACKAGE_VERSION,
		"CL", command,
		"DS", description,
		NULL
	) < 0) {
		free(command);
		mDie("Cannot add msamtools @PG record to SAM/BAM header");
	}

	free(command);
}

void mMultipleFileError(const char *subprogram, void **argtable) {
	fprintf(stderr, "Multiple input files not supported in %s.\n", subprogram);
	fprintf(stderr, "Use 'samtools merge' to combine BAM/SAM files.\n");
	mPrintHelp(subprogram, argtable);
	mQuit("");
}

/*
 * Handle sam file object
 */
mSamFile* mOpenSamInput(const char *filename, const char *inmode) {
	mSamFile *input;

	input = malloc(sizeof(mSamFile));
	if (input == NULL)
		mDie("Out of memory");

	input->file = sam_open(filename, inmode);
	if (input->file == NULL)
		mDie("Cannot open %s for reading", filename);

	input->header = sam_hdr_read(input->file);
	if (input->header == NULL)
		mDie("Cannot read header from %s", filename);

	input->replay_buffer = NULL;
	input->replay_count = 0;
	input->replay_index = 0;

	return input;
}

mSamFile* mOpenSamOutput(const char *filename, const char *outmode,
                         const sam_hdr_t *header) {
	mSamFile *output;

	output = malloc(sizeof(mSamFile));
	if (output == NULL)
		mDie("Out of memory");

	output->file = sam_open(filename, outmode);
	if (output->file == NULL)
		mDie("Cannot open %s for writing", filename);

	output->header = sam_hdr_dup(header);
	if (output->header == NULL)
		mDie("Cannot duplicate SAM header");

	output->replay_buffer = NULL;
	output->replay_count = 0;
	output->replay_index = 0;

	if (strchr(outmode, 'b') != NULL || strchr(outmode, 'h') != NULL) {
		if (sam_hdr_write(output->file, output->header) < 0)
			mDie("Cannot write SAM header");
	}

	return output;
}

int mSamRead(mSamFile *input, bam1_t *b) {
	if (input->replay_index < input->replay_count) {
		bam1_t *buffered = input->replay_buffer[input->replay_index];

		if (bam_copy1(b, buffered) == NULL)
			mDie("Out of memory while replaying buffered SAM/BAM record");

		bam_destroy1(buffered);
		input->replay_buffer[input->replay_index] = NULL;
		input->replay_index++;

		if (input->replay_index == input->replay_count) {
			free(input->replay_buffer);
			input->replay_buffer = NULL;
			input->replay_count = 0;
			input->replay_index = 0;
		}

		return 0;
	}

	return sam_read1(input->file, input->header, b);
}

int mSamWrite(mSamFile *output, bam1_t *b) {
	return sam_write1(output->file, output->header, b);
}

int mSamClose(mSamFile *stream) {
	int ret;
	size_t i;

	if (stream == NULL)
		return 0;

	if (stream->replay_buffer != NULL) {
		for (i=stream->replay_index; i<stream->replay_count; i++) {
			if (stream->replay_buffer[i] != NULL)
				bam_destroy1(stream->replay_buffer[i]);
		}
		free(stream->replay_buffer);
	}

	ret = sam_close(stream->file);
	sam_hdr_destroy(stream->header);
	free(stream);
	return ret;
}

mQNameCheckResult mSamCheckQNameGrouping(mSamFile *input) {
	mQNameCheckResult result;
	kstring_t        sort_order = {0, 0, NULL};
	bam1_t          *record;
	zoeHash          closed_qnames;
    zoeVec           allocated_groups;
	char            *current_qname = NULL;
	size_t           current_group_first = 0;
	size_t           record_number = 0;
	size_t           i;
	int              read_status;
	int              coordinate_ordered = 1;
	int              coordinate_relevant = 0;
	int              have_previous_mapped = 0;
	int32_t          previous_tid = -1;
	hts_pos_t        previous_pos = -1;

	result.status = MSAM_QNAME_SAMPLE_OK;
	result.qname_records_checked = 0;
	result.input_records_checked = 0;
	result.mapped_records_checked = 0;

	/* Explicit header declarations are authoritative. */
	read_status = sam_hdr_find_tag_hd(input->header, "SO", &sort_order);
	if (read_status >= 0 && sort_order.s != NULL) {
		if (strcmp(sort_order.s, "queryname") == 0) {
			result.status = MSAM_QNAME_HEADER_CONFIRMED;
			free(sort_order.s);
			return result;
		}
		if (strcmp(sort_order.s, "coordinate") == 0) {
			free(sort_order.s);
			mDie(
				"Input SAM/BAM declares 'SO:coordinate', but this operation "
				"requires records to be grouped by QNAME.\n"
				"             Please name-sort the input, for example with "
				"'samtools sort -n input.bam -o input.name_sorted.bam'."
			);
		}
	}
	free(sort_order.s);

	record = bam_init1();
	if (record == NULL)
		mDie("Out of memory while checking QNAME grouping");

	input->replay_buffer = (bam1_t**) calloc(
		COORD_ORDER_CHECK_RECORDS, sizeof(bam1_t*)
	);
	if (input->replay_buffer == NULL) {
		bam_destroy1(record);
		mDie("Out of memory while buffering SAM/BAM records");
	}
	input->replay_count = 0;
	input->replay_index = 0;

	closed_qnames    = zoeNewHash();
	allocated_groups = zoeNewVec();

	for (i=0; i<COORD_ORDER_CHECK_RECORDS; i++) {
		const char *qname;

		read_status = sam_read1(input->file, input->header, record);
		if (read_status < 0)
			break;

		record_number++;
		result.input_records_checked++;

		input->replay_buffer[input->replay_count] = bam_dup1(record);
		if (input->replay_buffer[input->replay_count] == NULL)
			mDie("Out of memory while buffering SAM/BAM record");
		input->replay_count++;

		qname = bam_get_qname(record);

		/*
		 * A QNAME that reappears after its contiguous group has closed is
		 * definitive evidence that the input is unsuitable.
		 */
		if (record_number <= QNAME_GROUP_CHECK_RECORDS) {
			result.qname_records_checked++;

			if (current_qname == NULL) {
				current_qname = (char*) mMalloc(strlen(qname) + 1);
				strcpy(current_qname, qname);
				current_group_first = record_number;
			} else if (strcmp(qname, current_qname) != 0) {
				mClosedQNameGroup *closed;
				mClosedQNameGroup *reopened;

				closed = (mClosedQNameGroup*) mMalloc(
					sizeof(mClosedQNameGroup)
				);
				closed->first_record = current_group_first;
				closed->last_record = record_number - 1;

				zoeSetHash(closed_qnames, current_qname, closed);
				zoePushVec(allocated_groups, closed);

				reopened = (mClosedQNameGroup*) zoeGetHash(
					closed_qnames, qname
				);
				if (reopened != NULL) {
					size_t intervening = (
						record_number - reopened->last_record - 1
					);
					mDie(
						"SAM/BAM file is not grouped by QNAME. "
						"Read '%s' reappears at record %zu after its previous "
						"group ended at record %zu (%zu intervening records). "
						"Please name-sort the input, for example with "
						"'samtools sort -n input.bam -o input.name_sorted.bam'.",
						qname,
						record_number,
						reopened->last_record,
						intervening
					);
				}

				mFree(current_qname);
				current_qname = (char*) mMalloc(strlen(qname) + 1);
				strcpy(current_qname, qname);
				current_group_first = record_number;
			}
		}

		/*
		 * Coordinate ordering is only a heuristic.  Inspect mapped records
		 * and warn, rather than fail, if the sampled prefix is perfectly
		 * coordinate ordered and the records indicate that templates may
		 * have multiple SAM/BAM records.
		 */
		if (!(record->core.flag & BAM_FUNMAP) && record->core.tid >= 0) {
			result.mapped_records_checked++;

			if (have_previous_mapped) {
				if (
					record->core.tid < previous_tid ||
					(
						record->core.tid == previous_tid &&
						record->core.pos < previous_pos
					)
				) {
					coordinate_ordered = 0;
				}
			}

			previous_tid = record->core.tid;
			previous_pos = record->core.pos;
			have_previous_mapped = 1;
		}

		if (
			record->core.flag &
			(BAM_FPAIRED | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)
		) {
			coordinate_relevant = 1;
		}
	}

	if (current_qname != NULL)
		mFree(current_qname);

	for (i=0; i<(size_t)allocated_groups->size; i++) {
		mFree(allocated_groups->elem[i]);
	}
	zoeDeleteVec(allocated_groups);
	zoeDeleteHash(closed_qnames);
	bam_destroy1(record);

	if (input->replay_count == 0) {
		free(input->replay_buffer);
		input->replay_buffer = NULL;
	}

	if (
		coordinate_ordered &&
		coordinate_relevant &&
		result.mapped_records_checked >= COORD_ORDER_MIN_RECORDS
	) {
		char warning[1024];

		result.status = MSAM_QNAME_SAMPLE_WARNING;
		mFormatQNameCheck(&result, warning, sizeof(warning));
		fprintf(stderr, "WARNING: %s\n", warning);
	}

	return result;
}
