#include "mMatrix.h"

/* BEGIN MATRIX */

/***********
 * memory management
 ***********/

void mInitMatrix(mMatrix *m, int nrows, int ncols) {
	int i;
	m->elem      = (double**) mMalloc(nrows*sizeof(double*));
	m->row_names = (char**) mMalloc(nrows*sizeof(char*));
	m->col_names = (char**) mMalloc(ncols*sizeof(char*));
	for (i=0; i<nrows; i++)
		m->elem[i] = (double*) mMalloc(ncols*sizeof(double));
	m->nrows = nrows;
	m->ncols = ncols;
}

void mFreeMatrix(mMatrix *m) {
	int i;
	for (i=0; i<m->nrows; i++) {
		mFree(m->elem[i]);
		mFree(m->row_names[i]);
	}
	for (i=0; i<m->ncols; i++) {
		mFree(m->col_names[i]);
	}
	mFree(m->row_names);
	mFree(m->col_names);
	mFree(m->elem);
}

/***********
 * Matrix operations
 ***********/

mMatrix* mMatrixTranspose(mMatrix *a) {
	int i, j;
	mMatrix *tra = (mMatrix*) mMalloc(sizeof(mMatrix));

	mInitMatrix(tra, a->ncols, a->nrows);
	for (i=0; i<a->nrows; i++) {
		for (j=0; j<a->ncols; j++) {
			tra->elem[j][i] = a->elem[i][j];
		}
	}
	return tra;
}

mMatrix* mMatrixElementsSquared(mMatrix *a) {
	int i, j;
	mMatrix *sqr = (mMatrix*) mMalloc(sizeof(mMatrix));

	mInitMatrix(sqr, a->nrows, a->ncols);
	for (i=0; i<a->nrows; i++) {
		for (j=0; j<a->ncols; j++) {
			sqr->elem[i][j] = a->elem[i][j]*a->elem[i][j];
		}
	}
	return sqr;
}

void mSquareMatrixElements(mMatrix *a) {
	int i, j;
	for (i=0; i<a->nrows; i++) {
		for (j=0; j<a->ncols; j++) {
			a->elem[i][j] = a->elem[i][j]*a->elem[i][j];
		}
	}
}

void mAddToMatrix(mMatrix *a, mMatrix *b) {
	int i, j;
	if (a->nrows != b->nrows || a->ncols != b->ncols) {
		mDie("Incompatible matrices in addition");
	}
	for (i=0; i<a->nrows; i++) {
		for (j=0; j<a->ncols; j++) {
			a->elem[i][j] += b->elem[i][j];
		}
	}
}

mMatrix* mAddMatrices(mMatrix *a, mMatrix *b) {
	mMatrix* sum;
	int i, j;
	if (a->nrows != b->nrows || a->ncols != b->ncols) {
		mDie("Incompatible matrices in addition");
	}
	sum = (mMatrix*) mMalloc(sizeof(mMatrix));
	mInitMatrix(sum, a->nrows, a->ncols);
	for (i=0; i<sum->nrows; i++) {
		for (j=0; j<sum->ncols; j++) {
			sum->elem[i][j] = a->elem[i][j] + b->elem[i][j];
		}
	}
	return sum;
}

mMatrix* mSubtractMatrices(mMatrix *a, mMatrix *b) {
	mMatrix* diff;
	int i, j;
	if (a->nrows != b->nrows || a->ncols != b->ncols) {
		mDie("Incompatible matrices in subtraction");
	}
	diff = (mMatrix*) mMalloc(sizeof(mMatrix));
	mInitMatrix(diff, a->nrows, a->ncols);
	for (i=0; i<diff->nrows; i++) {
		for (j=0; j<diff->ncols; j++) {
			diff->elem[i][j] = a->elem[i][j] - b->elem[i][j];
		}
	}
	return diff;
}

mMatrix* mMultiplyMatrices(mMatrix *a, mMatrix *b) {
	mMatrix* product;
	int i, j, k;
	if (a->ncols != b->nrows) {
		mDie("Incompatible matrices in multiplication");
	}
	product = (mMatrix*) mMalloc(sizeof(mMatrix));
	mInitMatrix(product, a->nrows, b->ncols);
	for (i=0; i<product->nrows; i++) {
		for (j=0; j<product->ncols; j++) {
			double sum = 0;
			for (k=0; k<a->ncols; k++) {
				sum += a->elem[i][k]*b->elem[k][j];
			}
			product->elem[i][j] = sum;
		}
	}
	return product;
}

void mMultiplyMatrixByScalar(mMatrix *m, double d) {
	int i, j;
	for (i=0; i<m->nrows; i++) {
		for (j=0; j<m->ncols; j++) {
			m->elem[i][j] *= d;
		}
	}
}

void mDivideMatrixByScalar(mMatrix *m, double d) {
	int i, j;
	for (i=0; i<m->nrows; i++) {
		for (j=0; j<m->ncols; j++) {
			m->elem[i][j] /= d;
		}
	}
}

void mColumnNormalizeMatrix(mMatrix *m) {
	int i, j;
	for (j=0; j<m->ncols; j++) {
		double sum = 0;
		for (i=0; i<m->nrows; i++) {
			sum += m->elem[i][j];
		}
		for (i=0; i<m->nrows; i++) {
			m->elem[i][j] /= sum;
		}
	}
}

void mRowNormalizeMatrix(mMatrix *m) {
	int i, j;
	for (i=0; i<m->nrows; i++) {
		double sum = 0;
		for (j=0; j<m->ncols; j++) {
			sum += m->elem[i][j];
		}
		for (j=0; j<m->ncols; j++) {
			m->elem[i][j] /= sum;
		}
	}
}

void mFillMatrix(mMatrix *m, double x) {
	int i, j;
	for (i=0; i<m->nrows; i++) {
		for (j=0; j<m->ncols; j++) {
			m->elem[i][j] = x;
		}
	}
}

void mZerofyMatrix(mMatrix *m) {
	mFillMatrix(m, 0);
}

void mLogTransformMatrix(mMatrix *m) {
	int i, j;
	for (i=0; i<m->nrows; i++) {
		for (j=0; j<m->ncols; j++) {
			m->elem[i][j] = log(m->elem[i][j]);
		}
	}
}

/***********
 * inter-matrix methods
 ***********/
double mSquaredEuclideanDistance(mMatrix *a, mMatrix *b) {
	int i;
	double sum = 0;
	if (a->nrows != b->nrows || a->ncols != b->ncols || a->ncols != 1) {
		mDie("Cannot calculate Euclidean distances for %d X %d matrices", a->nrows, a->ncols);
	}
	for (i=0; i<a->nrows; i++) {
		double dist = (a->elem[i][0] - b->elem[i][0]);
		sum += dist*dist;
	}
	return sum;
}

double mEuclideanDistance(mMatrix *a, mMatrix *b) {
	return sqrt(mSquaredEuclideanDistance(a, b));
}

double mBinaryEuclideanDistance(mMatrix *a, mMatrix *b) {
	int i;
	double sum = 0;
	if (a->nrows != b->nrows || a->ncols != b->ncols || a->ncols != 1) {
		mDie("Cannot calculate Euclidean distances for %d X %d matrices", a->nrows, a->ncols);
	}
	for (i=0; i<a->nrows; i++) {
		double dist = (MIN(1,a->elem[i][0]) - MIN(1,b->elem[i][0]));
		sum += dist*dist;
	}
	return sqrt(sum);
}

double mKullbackLeiblerDivergence(mMatrix *a, mMatrix *b) {
	int i;
	double sum = 0;
	if (a->nrows != b->nrows || a->ncols != b->ncols || a->ncols != 1) {
		mDie("Cannot calculate Kullback-Leibler divergence for %d X %d matrices", a->nrows, a->ncols);
	}
	for (i=0; i<a->nrows; i++) {
		sum += a->elem[i][0]*log(a->elem[i][0]/b->elem[i][0]);
	}
	return sum;
}

double mJensenShannonDivergence(mMatrix *a, mMatrix *b) {
	double divergence;
	mMatrix* mean;

	if (a->nrows != b->nrows || a->ncols != b->ncols || a->ncols != 1) {
		mDie("Cannot calculate Jensen-Shannon divergence for %d X %d matrices", a->nrows, a->ncols);
	}
	
	mean = mAddMatrices(a, b);
	mDivideMatrixByScalar(mean, 2);
	divergence = (mKullbackLeiblerDivergence(a, mean) + mKullbackLeiblerDivergence(b, mean))/2;
	mFreeMatrix(mean);
	mFree(mean);
	return divergence;
}

double mJensenShannonDistance(mMatrix *a, mMatrix *b) {
	return sqrt(mJensenShannonDivergence(a, b));
}

/***********
 * Matrix read/write methods
 ***********/
void mMatrixFromArray(mMatrix *m, int nrows, int ncols, double **array) {
	int i, j;
	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			m->elem[i][j] = array[i][j];
		}
	}
}

void mReadMatrix(FILE *stream, mMatrix *m, int nrows, int ncols) {
	int i, j;
	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			if (fscanf(stream, "%lf", m->elem[i]+j) == EOF) {
				perror("mReadMatrix encountered unexpected EOF");
			}
		}
	}
}

/* label_source&01 == 01 => each row has label for itself, label_source&10 == 10 => top row has labels for every column */
void mReadMatrixWithLabels(FILE *stream, mMatrix *m, int nrows, int ncols, int label_source, mVector *labels) {
	int i, j;
	char *label;
	char *buffer = (char*) mMalloc(64*sizeof(char));
	if (label_source == 2) {
		for (i=0; i<ncols; i++) {
			if (fscanf(stream, "%s", buffer) == EOF) {
				perror("mReadMatrixWithLabels encountered unexpected EOF while reading Labels");
			}
			label = (char*) mMalloc(strlen(buffer)*sizeof(char));
			mPushVector(labels, label);
		}
	}
	for (i=0; i<nrows; i++) {
		if (label_source == 1) {
			if (fscanf(stream, "%s", buffer) == EOF) {
				perror("mReadMatrixWithLabels encountered unexpected EOF while reading Labels");
			}
			label = (char*) mMalloc(strlen(buffer)*sizeof(char));
			mPushVector(labels, label);
		}
		for (j=0; j<ncols; j++) {
			if (fscanf(stream, "%lf", m->elem[i]+j) == EOF) {
				perror("mReadMatrixWithLabels encountered unexpected EOF while reading Data");
			}
		}
	}
	mFree(buffer);
}

void mWriteMatrix(FILE *stream, mMatrix *m) {
	int i, j;
	for (i=0; i<m->nrows; i++) {
		fprintf(stream, "| ");
		for (j=0; j<m->ncols; j++) {
			fprintf(stream, "%-8.8f ", m->elem[i][j]);
		}
		fprintf(stream, "|\n");
	}
}

void mReadRMatrix(FILE *stream, mMatrix *m, int nrows, int ncols, int header, int row_names) {
	int i, j;
	char *label;
	char *buffer = (char*) mMalloc(128*sizeof(char));
	if (header == 1) {
		for (i=0; i<ncols; i++) {
			if (fscanf(stream, "%s", buffer) == EOF) {
				perror("mReadRMatrix encountered unexpected EOF while reading Labels");
			}
			label = (char*) mMalloc((strlen(buffer)+1)*sizeof(char));
			strcpy(label, buffer);
			m->col_names[i] = label;
		}
	}
	for (i=0; i<nrows; i++) {
		if (row_names == 1) {
			if (fscanf(stream, "%s", buffer) == EOF) {
				perror("mReadRMatrix encountered unexpected EOF while reading Labels");
			}
			label = (char*) mMalloc((strlen(buffer)+1)*sizeof(char));
			strcpy(label, buffer);
			m->row_names[i] = label;
		}
		for (j=0; j<ncols; j++) {
			if (fscanf(stream, "%lf", m->elem[i]+j) == EOF) {
				perror("mReadRMatrix encountered unexpected EOF while reading Data");
			}
		}
	}
	mFree(buffer);
}

void mWriteRMatrix(FILE *stream, mMatrix *m) {
	int i, j;
	fprintf(stream, "%s", m->col_names[0]);
	for (j=1; j<m->ncols; j++) {
		fprintf(stream, "\t%s", m->col_names[j]);
	}
	fprintf(stream, "\n");
	for (i=0; i<m->nrows; i++) {
		fprintf(stream, "%s", m->row_names[i]);
		for (j=0; j<m->ncols; j++) {
			fprintf(stream, "\t%.8g", m->elem[i][j]);
		}
		fprintf(stream, "\n");
	}
}

void mWriteRMatrixTransposed(FILE *stream, mMatrix *m) {
	int i, j;
	fprintf(stream, "%s", m->row_names[0]);
	for (j=1; j<m->nrows; j++) {
		fprintf(stream, "\t%s", m->row_names[j]);
	}
	fprintf(stream, "\n");
	for (i=0; i<m->ncols; i++) {
		fprintf(stream, "%s", m->col_names[i]);
		for (j=0; j<m->nrows; j++) {
			fprintf(stream, "\t%.8g", m->elem[j][i]);
		}
		fprintf(stream, "\n");
	}
}

/* IMatrix */

/***********
 * memory management
 ***********/

void mInitMatrixI(mMatrixI *m, int nrows, int ncols) {
	int i;
	m->elem      = (int**) mMalloc(nrows*sizeof(int*));
	m->row_names = (char**) mMalloc(nrows*sizeof(char*));
	m->col_names = (char**) mMalloc(ncols*sizeof(char*));
	for (i=0; i<nrows; i++)
		m->elem[i] = (int*) mMalloc(ncols*sizeof(int));
	m->nrows = nrows;
	m->ncols = ncols;
}

void mFreeMatrixI(mMatrixI *m) {
	int i;
	for (i=0; i<m->nrows; i++) {
		mFree(m->elem[i]);
		mFree(m->row_names[i]);
	}
	mFree(m->row_names);
	for (i=0; i<m->ncols; i++) {
		mFree(m->col_names[i]);
	}
	mFree(m->col_names);
	mFree(m->elem);
}

void mReadRMatrixI(FILE *stream, mMatrixI *m, int nrows, int ncols, int header, int row_names) {
	int i, j;
	char *label;
	char *buffer = (char*) mMalloc(128*sizeof(char));
	if (header == 1) {
		for (i=0; i<ncols; i++) {
			if (fscanf(stream, "%s", buffer) == EOF) {
				perror("mReadRMatrixI encountered unexpected EOF while reading Labels");
			}
			label = (char*) mMalloc((strlen(buffer)+1)*sizeof(char));
			strcpy(label, buffer);
			m->col_names[i] = label;
		}
	}
	for (i=0; i<nrows; i++) {
		int *elem = m->elem[i];
		if (row_names == 1) {
			if (fscanf(stream, "%s", buffer) == EOF) {
				perror("mReadRMatrixI encountered unexpected EOF while reading Labels");
			}
			label = (char*) mMalloc((strlen(buffer)+1)*sizeof(char));
			strcpy(label, buffer);
			m->row_names[i] = label;
		}
		for (j=0; j<ncols; j++) {
			if (fscanf(stream, "%d", elem+j) == EOF) {
				perror("mReadRMatrixI encountered unexpected EOF while reading Data");
			}
		}
	}
	mFree(buffer);
}

void mWriteRMatrixI(FILE *stream, mMatrixI *m) {
	int i, j;
	fprintf(stream, "%s", m->col_names[0]);
	for (j=1; j<m->ncols; j++) {
		fprintf(stream, "\t%s", m->col_names[j]);
	}
	fprintf(stream, "\n");
	for (i=0; i<m->nrows; i++) {
		int *elem = m->elem[i];
		fprintf(stream, "%s", m->row_names[i]);
		for (j=0; j<m->ncols; j++) {
			fprintf(stream, "\t%d", elem[j]);
		}
		fprintf(stream, "\n");
	}
}

void mWriteRMatrixTransposedI(FILE *stream, mMatrixI *m) {
	int i, j;
	fprintf(stream, "%s", m->row_names[0]);
	for (j=1; j<m->nrows; j++) {
		fprintf(stream, "\t%s", m->row_names[j]);
	}
	fprintf(stream, "\n");
	for (i=0; i<m->ncols; i++) {
		fprintf(stream, "%s", m->col_names[i]);
		for (j=0; j<m->nrows; j++) {
			fprintf(stream, "\t%d", m->elem[j][i]);
		}
		fprintf(stream, "\n");
	}
}

/* END_MATRIX */
