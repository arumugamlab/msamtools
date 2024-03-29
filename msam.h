#ifndef M_SAM_H
#define M_SAM_H

#include <stdlib.h>
#include <limits.h>
#include <regex.h>
#include <zlib.h>
#include "argtable2.h"
#include "bam.h"
#include "sam.h"
#include "mBamVector.h"
#include "mMatrix.h"
#include "zoeTools.h"

/* I use integers for marking per-position coverage */
typedef int coverage_t;
#define COVERAGE_T_FORMAT "%d"

struct msam_global {
	int            multiple_input;

	/* features to count vs sequences in bam file */
	bam_header_t  *header;
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

/* global variables grouped into an object */

extern msam_global *global;

/* Helper functions */

samfile_t* mOpenSamFile(const char *filename, const char *inmode, char *headerfile);
void       mMultipleFileError(const char *subprogram, void **argtable);
int        mHeaderCheck(int count, const char* filenames[], const char *inmode);

zoeHash mReadIntegerHash(const char *filename);
void    mEmptyZoeHash(zoeHash hash);

void mInitGlobal();
void mFreeGlobal();

FILE* mInitOutputStream(const char *filename, int gzip);
void  mFreeOutputStream(FILE *output, int gzip);

#define M_INPUT_MODE(arg_samin) ((arg_samin->count == 0)?"rb":"r")

void mPrintHelp (const char *subprogram, void **argtable);
void mPrintCommandLine(FILE *output, int argc, char *argv[]);
void mPrintCommandLineGzip(gzFile output, int argc, char *argv[]);

/* Main functions for the subprograms */

int msam_filter_main(int argc, char* argv[]);
int msam_profile_main(int argc, char* argv[]);
int msam_summary_main(int argc, char* argv[]);
int msam_coverage_main(int argc, char* argv[]);

#endif
