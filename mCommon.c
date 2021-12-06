#include "mCommon.h"

void mQuit(const char *format, ...) {
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit(EXIT_SUCCESS);
}

void mWarn(const char *format, ...) {
	va_list args;
	fflush(stdout);
	fprintf(stderr, "Warning: ");
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
}

void mDie(const char *format, ...) {
	va_list args;
	fflush(stdout);
	fprintf(stderr, "Fatal Error: ");
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

void PROGRESS_PRINT(const char *format, ...) {
	va_list args;
	fflush(stdout);
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fflush(stderr);
}

/***************************************************** 
 * File handling
 *****************************************************/

FILE* mSafeOpenFile(const char *filename, const char *mode, int gzip) {
	FILE *stream;
    if (gzip) {
        char command[512];
        if (strcmp(filename, "-") == 0) {
            if (strcmp(mode, "r") == 0) {
                sprintf(command, "gzip --stdout --decompress");
            } else if (strcmp(mode, "w") == 0) {
                sprintf(command, "gzip --stdout");
            }
            if ((stream = popen(command, mode)) == NULL) {
                mDie("Cannot open gzipped %s for mode %s", filename, mode);
            }
        } else {
            if (strcmp(mode, "r") == 0) {
                sprintf(command, "gzip --stdout --decompress %s", filename);
            } else if (strcmp(mode, "w") == 0) {
                sprintf(command, "gzip --stdout > %s", filename);
            }
            if ((stream = popen(command, mode)) == NULL) {
                mDie("Cannot open gzipped file %s for mode %s", filename, mode);
            }
        }
    } else {
        if (strcmp(filename, "-") == 0) {
            if (strcmp(mode, "r") == 0) {
                stream = stdin;
            } else if (strcmp(mode, "w") == 0) {
                stream = stdout;
            } else {
                mDie("Cannot open %s for mode %s", filename, mode);
            }
        } else {
            if ((stream = fopen(filename, mode)) == NULL) {
                mDie("Cannot open %s for mode %s", filename, mode);
            }
        }
    }
	return stream;
}

void mSafeCloseFile(FILE *stream, int gzip) {
	if (gzip) pclose(stream);
	else fclose(stream);
}

void mWriteUIntArrayBinary(FILE *stream, uint32_t size, uint32_t *array) {
	size_t bit32 = sizeof(uint32_t);
	fwrite(&size, bit32, 1, stream);
	fwrite(array, bit32, size, stream);
}

void mReadUIntArrayBinary(FILE *stream, uint32_t *size, uint32_t **array) {
	size_t bit32 = sizeof(uint32_t);
	size_t asize;
	if (fread(size, bit32, 1, stream) != bit32) {
		mDie("Could not read the header bytes in mUIntArrayBinary file");
	}
	asize = (*size)*bit32;
	*array = (uint32_t*) mMalloc(asize);
	if (fread(*array, bit32, *size, stream) != asize) {
		mDie("Could not read %d bytes in mUIntArrayBinary file", asize);
	}
}

/***************************************************** 
 * Text handling
 *****************************************************/

char* mCanonicalGroupName(const char *text, int length) {
	int i;
	char *ret   = (char*) mMalloc((length+1)*sizeof(char));
	char *slash = index(text, '/');
	if (slash == NULL) {
		strncpy(ret, text, length);
	} else {
		slash++; /* skip the slash */
		for (i=0; i<length/2; i++) {
			ret[i] = text[i];
		}
		for (i=length/2; i<length; i++) {
			ret[i] = slash[i-length/2];
		}
	}
	ret[length] = '\0';
	return ret;
}

void mGetFirstWord(char *source, char *word) {
	int i;
	int length = strlen(source);
	for (i=0; i<length; i++) {
		if (isspace(source[i])) {
			break;
		}
	}

	/* When there is no space, i=length, so it copies the whole string */

	strncpy(word, source, i);
	word[i] = '\0';
}

/* Get the mate-paired template name for Illumina platform.
 * Make sure that template has been malloc'ed with enough length.
 * No memory allocation here. It will just copy the template part
 * of the name into target.
 */

int mGetIlluminaTemplate(char *source, char *template) {
	int length = strlen(source);
	if (source[length-2] == '/') {
		strncpy(template, source, length-2);
		template[length-2] = '\0';
		return 1;
	} else {
		return 0;
	}
}

/***************************************************** 
 * Memory Allocation
 *****************************************************/

void* mMalloc(size_t size) {
	void *p;
	if ((p=malloc(size)) != NULL) {
		return p;
	}
	mQuit("malloc(%d) failed", size);
	return NULL;
}

void* mCalloc(size_t n, size_t size) {
	void *p;
	if ((p=calloc(n, size)) != NULL) {
		return p;
	}
	mQuit("calloc(%d, %d) failed", n, size);
	return NULL;
}

void* mRealloc(void *p, size_t size) {
	void *back = p;
	if ((p=realloc(p, size)) != NULL) {
		return p;
	}
	mFree(back);
	mQuit("realloc(%d) failed", size);
	return NULL;
}

void mFree(void *p) {
	free(p);
	p = NULL;
}

/***************************************************** 
 * Mathematical Functions
 *****************************************************/

int mIntPow(int x, int y) {
	return (int) pow((double)x, (double)y);
}

double mDoublePow(int x, int y) {
	return pow((double)x, (double)y);
}

long mLongPow(int x, int y) {
	return (long) pow((double)x, (double)y);
}

num_t mLog(num_t x) {
	if (x == 0) return -FMAX;
	else return log(x);
}

static const num_t  mLog2 = 0.301029995664;
num_t  mScore(num_t  x) {
	if (x < 0) mDie("Negative number found in mScore");

	if (x == 0) return -FMAX;
	else        return (num_t) 10*log(x)/mLog2;
}

/***************************************************** 
 * Useful data structures
 *****************************************************/

/***************************************************** 
 * Generic vector
 *****************************************************/

/*
Memory allocation to all the elements should be
done by the calling function. The vector just stores
it in void* form. You are responsible for freeing it
later.
*/

void mInitVector(mVector *vec, int limit) {
	vec->limit = (limit>0)?limit:1;
	/* Don't want to handle empty vectors. This is faster. Memory doesnt matter */
	vec->elem = (void**) mMalloc(limit*sizeof(void*));
	vec->size = 0;
	vec->last = NULL;
}

void mPushVector(mVector *vec, void *item) {
	if (vec->size == vec->limit) {
		vec->limit *= 2;
		vec->elem = (void**) mRealloc(vec->elem, vec->limit*sizeof(void*));
	}
	vec->elem[vec->size] = item;
	vec->last            = vec->elem[vec->size];
	vec->size++;
}

void mFreeVector(mVector *vec) {
	mFree(vec->elem);
}

/* Forcefully frees the pointers in the vector.
 * Use it with caution */
void mForceFreeVector(mVector *vec) {
	int i;
	for (i=0; i<vec->size; i++) {
		mFree(vec->elem[i]);
	}
	mFree(vec->elem);
}

/**********************
 ****** WARNING *******
 **********************
 * Don't mess with this unless you know what this does.
 * This relies on the fact that char is defined as always
 * having size of one byte. If that changes, this will fail
 */
void mSortVector(mVector *vec, size_t struct_size, int(*compar)(const void *, const void *)) {
	int i;
	int original_size = vec->size; /* vec->size cannot be used until we exit the function */
	char* elements = (char*) mMalloc(original_size*struct_size);
	for (i = 0; i < original_size; i++) {
		memcpy(elements+i*struct_size, vec->elem[i], struct_size);
		mFree(vec->elem[i]);
	}
	qsort(elements, original_size, struct_size, compar);

	mFreeVector(vec);
	mInitVector(vec, 2);
	for (i = 0; i < original_size; i++) {
		void *element = mMalloc(struct_size);
		memcpy(element, elements+i*struct_size, struct_size);
		mPushVector(vec, element);
	}
	mFree(elements);
}
/**********************
 ****** END WARNING ***
 **********************/

/***************************************************** 
 * mIVector - vector of integers
 *****************************************************/

void mInitIVector(mIVector *vec, int limit) {
	vec->limit = (limit>0)?limit:1;
	/* Don't want to handle empty vectors. This is faster. Memory doesnt matter */
	vec->elem = (int*) mMalloc(limit*sizeof(int));
	vec->size = 0;
}

void mPushIVector(mIVector *vec, int item) {
	if (vec->size == vec->limit) {
		vec->limit *= 2;
		vec->elem = (int*) mRealloc(vec->elem, vec->limit*sizeof(int));
	}
	vec->elem[vec->size] = item;
	vec->last       = vec->elem[vec->size];
	vec->size++;
}

void mFreeIVector(mIVector *vec) {
	mFree(vec->elem);
}

void mWriteIVector(FILE *fp, mIVector *vec, const char *delim) {
	int i;
	fprintf(fp, "<IVector size=\"%d\">\n", vec->size);
	for (i=0; i<vec->size-1; i++) {
		fprintf(fp, "%d%s", vec->elem[i], delim);
	}
	fprintf(fp, "%d\n", vec->elem[i]);
	fprintf(fp, "</IVector>\n");
}

void mFillIVector(mIVector *vec, int x) {
	int i;
	for (i=0; i<vec->size; i++)
		vec->elem[i] = x;
}

int  mBinarySearchIVector(mIVector *vec, int item) {
	int left = 0;
	int right = vec->size-1;
	int p, middle;
	do {
		p = (right-left)/2;
		middle = left+p;
		if (item < vec->elem[middle]) {
			right = middle;
		} else {
			if (item < vec->elem[middle+1]) {
				return middle;
			} else {
				left = middle;
			}
		}
	} while (p > 0);
	return -1;
}

/***************************************************** 
 * mLVector - vector of longs
 *****************************************************/

void mInitLVector(mLVector *vec, int limit) {
	vec->limit = (limit>0)?limit:1;
	/* Don't want to handle empty vectors. This is faster. Memory doesnt matter */
	vec->elem = (long*) mMalloc(limit*sizeof(long));
	vec->size = 0;
}

void mPushLVector(mLVector *vec, long item) {
	if (vec->size == vec->limit) {
		vec->limit *= 2;
		vec->elem = (long*) mRealloc(vec->elem, vec->limit*sizeof(long));
	}
	vec->elem[vec->size] = item;
	vec->last       = vec->elem[vec->size];
	vec->size++;
}

void mFreeLVector(mLVector *vec) {
	mFree(vec->elem);
}

void mWriteLVector(FILE *fp, mLVector *vec, const char *delim) {
	long i;
	fprintf(fp, "<LVector size=\"%d\">\n", vec->size);
	for (i=0; i<vec->size-1; i++) {
		fprintf(fp, "%ld%s", vec->elem[i], delim);
	}
	fprintf(fp, "%ld\n", vec->elem[i]);
	fprintf(fp, "</LVector>\n");
}

void mFillLVector(mLVector *vec, long x) {
	long i;
	for (i=0; i<vec->size; i++)
		vec->elem[i] = x;
}

int  mBinarySearchLVector(mLVector *vec, long item) {
	int left = 0;
	int right = vec->size-1;
	int p, middle;
	do {
		p = (right-left)/2;
		middle = left+p;
		if (item < vec->elem[middle]) {
			right = middle;
		} else {
			if (item < vec->elem[middle+1]) {
				return middle;
			} else {
				left = middle;
			}
		}
	} while (p > 0);
	return -1;
}

/***************************************************** 
 * mFVector - vector of floats
 *****************************************************/

void mInitFVector(mFVector *vec, int limit) {
	vec->limit = (limit>0)?limit:1;
	/* Don't want to handle empty vectors. This is faster. Memory doesnt matter */
	vec->elem = (num_t*) mMalloc(limit*sizeof(num_t));
	vec->size = 0;
}

void mPushFVector(mFVector *vec, num_t item) {
	if (vec->size == vec->limit) {
		vec->limit *= 2;
		vec->elem = (num_t*) mRealloc(vec->elem, vec->limit*sizeof(num_t));
	}
	vec->elem[vec->size] = item;
	vec->last       = vec->elem[vec->size];
	vec->size++;
}

void mFreeFVector(mFVector *vec) {
	mFree(vec->elem);
}

void mWriteFVector(FILE *fp, mFVector *vec, const char *delim) {
	int i;
	fprintf(fp, "<FVector size=\"%d\">\n", vec->size);
	for (i=0; i<vec->size-1; i++) {
		fprintf(fp, "%-2.10f%s", vec->elem[i], delim);
	}
	fprintf(fp, "%-2.10f\n", vec->elem[i]);
	fprintf(fp, "</FVector>\n");
}

void mFillFVector(mFVector *vec, num_t x) {
	int i;
	for (i=0; i<vec->size; i++)
		vec->elem[i] = x;
}

mFVector* mAddFVectors(mFVector *a, mFVector *b) {
	mFVector* sum;
	int i;
	if (a->size != b->size) {
		mDie("Incompatible vectors in addition: %d vs %d", a->size, b->size);
	}
	sum = (mFVector*) mMalloc(sizeof(mFVector));
	mInitFVector(sum, a->size);
	for (i=0; i<a->size; i++) {
		mPushFVector(sum, a->elem[i] + b->elem[i]);
	}
	return sum;
}

void mDivideFVectorByScalar(mFVector *a, num_t n) {
	int i;
	for (i=0; i<a->size; i++) 
		a->elem[i] /= n;
}

num_t mFKullbackLeiblerDivergence(mFVector *a, mFVector *b, int skip) {
	int i;
	num_t sum = 0;
	if (a->size != b->size) {
		mDie("Cannot calculate Kullback-Leibler divergence for incompatible vectors: %d vs %d", a->size, b->size);
	}
	for (i=skip; i<a->size; i++) {
		sum += a->elem[i]*log(a->elem[i]/b->elem[i]);
	}
	return sum;
}

num_t mFJensenShannonDivergence(mFVector *a, mFVector *b, int skip) {
	num_t divergence;
	mFVector* mean;

	mean = mAddFVectors(a, b);
	mDivideFVectorByScalar(mean, 2);
	divergence = (mFKullbackLeiblerDivergence(a, mean, skip) + mFKullbackLeiblerDivergence(b, mean, skip))/2;
	mFreeFVector(mean);
	mFree(mean);
	return divergence;
}

num_t mFJensenShannonDistance(mFVector *a, mFVector *b, int skip) {
	return sqrt(mFJensenShannonDivergence(a, b, skip));
}

/***************************************************** 
 * mDVector - vector of doubles
 *****************************************************/


/*
void mInitDVector(mDVector *vec, int limit) {
	vec->limit = (limit>0)?limit:1;
	vec->elem = (double*) mMalloc(limit*sizeof(double));
	vec->size = 0;
}

void mPushDVector(mDVector *vec, double item) {
	if (vec->size == vec->limit) {
		vec->limit *= 2;
		vec->elem = (double*) mRealloc(vec->elem, vec->limit*sizeof(double));
	}
	vec->elem[vec->size] = item;
	vec->last       = vec->elem[vec->size];
	vec->size++;
}

void mFreeDVector(mDVector *vec) {
	mFree(vec->elem);
}

*/

/***************************************************** 
 * comparison functions for qsort
 *****************************************************/

int mIntPtrCmp(const void *a, const void *b) {
	return *(int*)a - *(int*)b;
}

int mDoublePtrCmp(const void *a, const void *b) {
	double diff =  *(double*)a - *(double*)b;
	if (diff > 0)      return 1;
	else if (diff < 0) return -1;
	else               return 0;
}

int mFloatPtrCmp(const void *a, const void *b) {
	float diff = *(float*)a - *(float*)b;
	if (diff > 0)      return 1;
	else if (diff < 0) return -1;
	else               return 0;
}

int mNumTPtrCmp(const void *a, const void *b) {
	num_t diff = *(num_t*)a - *(num_t*)b;
	if (diff > 0)      return 1;
	else if (diff < 0) return -1;
	else               return 0;
}

/***************************************************** 
 * Statistical Analysis on arrays
 *****************************************************/

num_t mMean(int count, num_t *numbers) {
	int i;
	num_t sum = 0;
	for (i=0; i<count; i++) sum+=numbers[i];
	return (num_t) 1.0*sum/count;
}

num_t mStdev(int count, num_t *numbers, num_t *out_mean) {
	int i;
	num_t mean = mMean(count, numbers);
	num_t sum  = 0;
	for (i=0; i<count; i++) {
		num_t diff = numbers[i]-mean;
		sum += (diff*diff);
	}
	if (count == 1) {
		sum = 0;
	} else {
		sum = 1.0 * sum / (count-1);
	}
	if (out_mean != NULL)
		*out_mean = mean;
	return (num_t) sqrt(sum);
}

void mMinMax(int count, num_t *numbers, num_t *min, num_t *max) {
	int i;
	*min = FLT_MAX;
	*max = -FLT_MAX;
	for (i=0; i<count; i++) {
		num_t f = numbers[i];
		if (f > *max)
			*max = f;
		if (f < *min)
			*min = f;
	}
}

/*
 * Binary read/write of arrays
 */

/*
void mWriteIntArrayBinary(FILE *stream, int32_t length, uint16_t *data) {
}
*/

/***************************************************** 
 * Useful data structures
 *****************************************************/

/* Generic sequence handler without any struct constraints */

void mCountFastaFile(char *filename, long *entries, long *bases) {
	FILE *fp;
	long  base_count=0, entry_count=0;
	long *count = (long*) calloc(256, sizeof(long));
	char  nuc_symbols[] = "ACGTURYMKSWBDHVN";
	int   c, i;
	
	if ((fp = fopen(filename, "r")) == NULL) {
		mDie("Could not open file %s\n", filename);
	}
	while ((c = fgetc(fp)) != EOF){
		/* Skip the FASTA header */
		while ('>' == (char) c) {
			entry_count++;
			while ( (c = fgetc(fp)) != EOF && c != '\n' ) {
			}
			if ( c == EOF || (c = fgetc(fp)) == EOF) break;
		}
		c = toupper(c);
		count[c]++;
	}
	fclose(fp);
	for (i=0; i<(int)strlen(nuc_symbols); i++) {
		base_count += count[(int)nuc_symbols[i]];
	}
	if (bases != NULL)   *bases   = base_count;
	if (entries != NULL) *entries = entry_count;
	return;
}
