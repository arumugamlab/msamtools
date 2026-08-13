#ifndef M_SAM_H
#define M_SAM_H

#include <stdlib.h>
#include <limits.h>
#include <regex.h>
#include <zlib.h>
#include "argtable2.h"
#include <htslib/sam.h>
#include "mMatrix.h"
#include "zoeTools.h"

/* I use integers for marking per-position coverage */
typedef int coverage_t;
#define COVERAGE_T_FORMAT "%d"

/* global variables grouped into an object */

struct msam_global {
	int            multiple_input;

	/* features to count vs sequences in bam file */
	sam_hdr_t     *header;
	int           *fmap; /* feature map: seq-->feature*/
	int            n_features;
	uint32_t      *feature_len;
	char         **feature_name;
	

	/* COVERAGE */
	coverage_t   **coverage;        /* read coverage per target_seq position */
	int           *covered;          /* has this target been covered? */
	int          **seq_touched;      /* part of the sequence that is covered */

	/* PROFILE */
	uint32_t       uniq_mapper_count;
	uint32_t       multi_mapper_count;
	uint32_t       purged_insert_count;
	uint32_t      *ui_insert_count;
	double        *d_insert_count;
	mVector       *multi_mappers;    /* A vector of integer vectors */
	uint8_t       *ub_target_hit;    /* Hit flag for current multimapper read */

	/* Convenience in filtering */

	int32_t          MIN_LENGTH;
	int32_t          PPT;
	int32_t          MAX_CLIP;
};
typedef struct msam_global msam_global;

extern msam_global *global;

/* New wrapper struct to replace old samfile_t */

typedef struct {
	samFile    *file;
	sam_hdr_t  *header;
	bam1_t    **replay_buffer;
	size_t      replay_count;
	size_t      replay_index;
} mSamFile;

/* objects to track whether sam file is qname-grouped */

typedef enum {
	MSAM_QNAME_NOT_REQUIRED = 0,
	MSAM_QNAME_HEADER_CONFIRMED,
	MSAM_QNAME_SAMPLE_OK,
	MSAM_QNAME_SAMPLE_WARNING
} mQNameCheckStatus;

typedef struct {
	mQNameCheckStatus status;
	size_t            qname_records_checked;
	size_t            input_records_checked;
	size_t            mapped_records_checked;
} mQNameCheckResult;

/* Helper functions for global variables struct */

void mInitGlobal();
void mFreeGlobal();

/* Wrapper functions for mSamFile */

mSamFile* mOpenSamInput(const char *filename, const char *inmode);
mSamFile* mOpenSamOutput(const char *filename, const char *outmode, const sam_hdr_t *header);
int       mSamRead(mSamFile *input, bam1_t *b);
int       mSamWrite(mSamFile *output, bam1_t *b);
int       mSamClose(mSamFile *stream);

/* Helper functions for qname-grouping test */

mQNameCheckResult mQNameCheckNotRequired(void);
mQNameCheckResult mSamCheckQNameGrouping(mSamFile *input);

/* Other helper functions */

void mAddSamProvenance(sam_hdr_t *header, int argc, char *argv[], const mQNameCheckResult *qname_check);
void mPrintProfileProvenanceGzip(gzFile output, int argc, char *argv[], const mQNameCheckResult *qname_check);
void mMultipleFileError(const char *subprogram, void **argtable);
void mPrintHelp (const char *subprogram, void **argtable);

#define M_INPUT_MODE(arg_samin) ((arg_samin->count == 0)?"rb":"r")

/* Main functions for the subprograms */

int msam_filter_main(int argc, char* argv[]);
int msam_profile_main(int argc, char* argv[]);
int msam_summary_main(int argc, char* argv[]);
int msam_coverage_main(int argc, char* argv[]);

#endif
