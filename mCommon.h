#ifndef M_COMMON_H
#define M_COMMON_H
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdint.h>

void mQuit(const char *format, ...);
void mDie(const char *format, ...);
void mWarn(const char *format, ...);
void PROGRESS_PRINT(const char *format, ...);

/***************************************************** 
 * Constants
 *****************************************************/

/***************************************************** 
 * File handling
 *****************************************************/

FILE* mSafeOpenFile(const char *filename, const char *mode, int gzip);
void  mSafeCloseFile(FILE *stream, int gzip);
void mWriteUIntArrayBinary(FILE *stream, uint32_t size, uint32_t *array);
void mReadUIntArrayBinary(FILE *stream, uint32_t *size, uint32_t **array);

/***************************************************** 
 * Text handling
 *****************************************************/

char* mCanonicalGroupName(const char *text, int length);
void  mGetFirstWord(char *source, char *word);
int   mGetIlluminaTemplate(char *source, char *template);

/***************************************************** 
 * Memory Allocation
 *****************************************************/

void* mMalloc(size_t size);
void* mCalloc(size_t nelem, size_t size);
void* mRealloc(void *p, size_t size);
void  mFree(void *p);

/***************************************************** 
 * Mathematical Functions
 *****************************************************/

#define min(a,b)         ((a<b)?a:b)

typedef double num_t;
#define FMAX (DBL_MAX)

int    mIntPow(int x, int y);
num_t  mDoublePow(int x, int y);
long   mLongPow(int x, int y);
num_t  mLog(num_t  x);
num_t  mScore(num_t  x);
num_t mMean(int count, num_t *numbers);
num_t mStdev(int count, num_t *numbers, num_t *out_mean);
void mMinMax(int count, num_t *numbers, num_t *min, num_t *max);

/***************************************************** 
 * Useful data structures
 *****************************************************/

/*
Memory allocation to all the elements should be
done by the calling function. The vector just stores
it in void* form. You are responsible for freeing it
later.
*/

struct mVector {
	int    size;
	int    limit;
	void **elem;
	void  *last;
};
typedef struct mVector mVector;

void mInitVector(mVector *vec, int limit);
void mPushVector(mVector *vec, void *item);
void mFreeVector(mVector *vec);
void mForceFreeVector(mVector *vec);
void mSortVector(mVector *vec, size_t struct_size, int(*compar)(const void *, const void *));

struct mIVector {
	int  size;
	int  limit;
	int *elem;
	int  last;
};
typedef struct mIVector mIVector;

void mInitIVector(mIVector *vec, int limit);
void mPushIVector(mIVector *vec, int item);
void mFreeIVector(mIVector *vec);
void mWriteIVector(FILE *fp, mIVector *vec, const char *delim);
void mFillIVector(mIVector *vec, int x);
int  mBinarySearchIVector(mIVector *vec, int item);

struct mLVector {
	int  size;
	int  limit;
	long *elem;
	long  last;
};
typedef struct mLVector mLVector;

void mInitLVector(mLVector *vec, int limit);
void mPushLVector(mLVector *vec, long item);
void mFreeLVector(mLVector *vec);
void mWriteLVector(FILE *fp, mLVector *vec, const char *delim);
void mFillLVector(mLVector *vec, long x);
int  mBinarySearchLVector(mLVector *vec, long item);

struct mFVector {
	int    size;
	int    limit;
	num_t *elem;
	num_t  last;
};
typedef struct mFVector mFVector;

void mInitFVector(mFVector *vec, int limit);
void mPushFVector(mFVector *vec, num_t item);
void mFreeFVector(mFVector *vec);
void mWriteFVector(FILE *fp, mFVector *vec, const char *delim);
void mFillFVector(mFVector *vec, num_t x);

void mDivideFVectorByScalar(mFVector *a, num_t n);
mFVector* mAddFVectors(mFVector *a, mFVector *b);
num_t mFKullbackLeiblerDivergence(mFVector *a, mFVector *b, int skip);
num_t mFJensenShannonDistance(mFVector *a, mFVector *b, int skip);
num_t mFJensenShannonDivergence(mFVector *a, mFVector *b, int skip);

void mCountFastaFile(char *filename, long *entries, long *bases);

int mFloatPtrCmp(const void *a, const void *b);
int mIntPtrCmp(const void *a, const void *b);
int mDoublePtrCmp(const void *a, const void *b);

#endif
