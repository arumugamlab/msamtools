#ifndef M_MATRIX_H
#define M_MATRIX_H

#include "mDefinitions.h"
#include "mCommon.h"

struct mMatrix {
	int      nrows;
	int      ncols;
	double **elem;
	char   **row_names;
	char   **col_names;
};
typedef struct mMatrix mMatrix;

void mInitMatrix(mMatrix *m, int nrows, int ncols);
void mFreeMatrix(mMatrix *m);

void mMatrixFromArray(mMatrix *m, int nrows, int ncols, double **array);
void mReadMatrix(FILE *stream, mMatrix *m, int nrows, int ncols);
void mReadMatrixWithLabels(FILE *stream, mMatrix *m, int nrows, int ncols, int label_source, mVector *labels);
void mWriteMatrix(FILE *stream, mMatrix *m);
void mReadRMatrix(FILE *stream, mMatrix *m, int nrows, int ncols, int header, int row_names);
void mWriteRMatrix(FILE *stream, mMatrix *m);
void mWriteRMatrixTransposed(FILE *stream, mMatrix *m);

/* destructive - change the matrix that is passed */
void mAddToMatrix(mMatrix *a, mMatrix *b);
void mDivideMatrixByScalar(mMatrix *m, double n);
void mZerofyMatrix(mMatrix *m);
void mFillMatrix(mMatrix *m, double x);
void mSquareMatrixElements(mMatrix *a);
void mLogTransformMatrix(mMatrix *m);
void mColumnNormalizeMatrix(mMatrix *m);
void mRowNormalizeMatrix(mMatrix *m);

/* non-destructive - arguments are untouched and the result is returned */
mMatrix* mMatrixTranspose(mMatrix *a);
mMatrix* mAddMatrices(mMatrix *a, mMatrix *b);
mMatrix* mSubtractMatrices(mMatrix *a, mMatrix *b);
mMatrix* mMultiplyMatrices(mMatrix *a, mMatrix *b);
mMatrix* mMatrixElementsSquared(mMatrix *a);
mMatrix* mMatrixInverse(mMatrix *m);

/* mMatrixI */
struct mMatrixI {
	int      nrows;
	int      ncols;
	int    **elem;
	char   **row_names;
	char   **col_names;
};
typedef struct mMatrixI mMatrixI;

void mInitMatrixI(mMatrixI *m, int nrows, int ncols);
void mFreeMatrixI(mMatrixI *m);
void mReadRMatrixI(FILE *stream, mMatrixI *m, int nrows, int ncols, int header, int row_names);
void mWriteRMatrixI(FILE *stream, mMatrixI *m);
void mWriteRMatrixTransposedI(FILE *stream, mMatrixI *m);

double mSquaredEuclideanDistance(mMatrix *a, mMatrix *b);
double mEuclideanDistance(mMatrix *a, mMatrix *b);
double mBinaryEuclideanDistance(mMatrix *a, mMatrix *b);
double mJensenShannonDistance(mMatrix *a, mMatrix *b);
double mKullbackLeiblerDivergence(mMatrix *a, mMatrix *b);

mMatrix* mMahalanobisDistance(mMatrix *a, mMatrix *b, mMatrix *covar);

#endif
