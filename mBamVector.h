#include <inttypes.h>
#include "mCommon.h"
#include "bam.h"
#include "sam.h"
#include "htslib/kstring.h"

/*
 *
 * bam1_t related
 *
 */

#define bam1_qfrag(b) (((char*)b->data)[strlen((char*)b->data)-1] - '1')
#define target_taxlen(x) (strchr(x, '.') - x)

/****
 *  2019.10.26: bam1_templatecmp was written for mate pairs with r001/1 & r001/2
 *              as names, and these names used to appear in the SAM files.
 *              SAM specification clearly defines QNAME to be template name.
 *              Therefore, bam1_templatecmp has become unnecessary.
 *              I leave the function for backward compatibility, but it is the 
 *              same as comparing full QNAMEs.
 */
#define bam1_templatecmp(a,b) (strcmp(bam1_qname(a), bam1_qname(b)))

#define M_BAM_FOUTIE   2048
#define M_BAM_FUNIDIR  4096

#define M_BAM_MAX_READ_LENGTH (4096)

/* get the length of the alignment (i.e., excluding the hard/soft masks in the ends) */

int32_t bam_cigar2alnlen(const bam1_core_t *c, const uint32_t *cigar);

/* get the details of the alignment */

void bam_cigar2details(const bam1_core_t *c, const uint32_t *cigar, int32_t *alen, int32_t *qlen, int32_t *qclip);

int32_t bam_highquality_differences_only(const bam1_t *b, uint8_t minqual);

struct mAlignmentSummary {
	int32_t match;
	int32_t mismatch;
	int32_t gapopen;
	int32_t gapextend;
	int32_t length;
	int32_t query_length;
	int32_t query_clip;
	int32_t edit;
};
typedef struct mAlignmentSummary mAlignmentSummary;

void bam_get_summary(const bam1_t *b, mAlignmentSummary *summary);
void bam_get_extended_summary(const bam1_t *b, mAlignmentSummary *summary);

/*
 *
 * mBamVector
 *
 */

struct mBamVector {
	int size;
	int limit;
	bam1_t **elem;
};
typedef struct mBamVector mBamVector;

void mInitBamVector(mBamVector *vec, int limit);
void mPushBamVector(mBamVector *vec, bam1_t *item);
void mExpandBamVector(mBamVector *vec);
void mFreeBamVector(mBamVector *vec);
void mEmptyBamVector(mBamVector *vec);
void mSortBamVector(mBamVector *vec, int(*compar)(const void *, const void *));
void mWriteBamVector(samfile_t *stream, mBamVector *bamvector);

/*
 *
 * mBamPool
 *
 */

struct mBamPool {
	int limit;     /* Max elements handled at the moment */
	int size;   /* Current size of the pool */
	bam1_t **elem;
};
typedef struct mBamPool mBamPool;

#define pool_current(a) (a->elem[a->size])

void mInitBamPool(mBamPool *pool, int limit);
void mExpandBamPool(mBamPool *pool);
void mFreeBamPool(mBamPool *pool);
bam1_t* mAdvanceBamPool(mBamPool *pool);
void mReOriginateBamPool(mBamPool *pool);
void mWriteBamPool(samfile_t *stream, mBamPool *pool);
