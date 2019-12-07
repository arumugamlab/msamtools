#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mSequence.h"

void mInitSeq(mSeq *seq) {
	char *type_string = (char*) mCalloc(32, sizeof(char));
	switch (seq->type) {
		case FASTA_ALPH:
		case FASTA_DNA:
		case FASTA_PROT:
			strcpy(type_string, "fasta");
			seq->window = 60;
			break;
		case FASTQ:
			strcpy(type_string, "fastq");
			seq->window = 500;
			break;
		case FASTA_QUAL:
			strcpy(type_string, "qual");
			seq->window = 17;
			break;
		default:
			mDie("Fatal error: mInitSeq called without setting seqtype!");
			break;
	}
	seq->type_string = type_string;
	seq->def  = NULL;
	seq->alph = NULL;
	seq->qual = NULL;
	seq->code = NULL;
	seq->aux  = NULL;
	seq->aux_count = 0;
}


/******************************
 *  Seq freeing  functions    *
 ******************************/

void mFreeSeq(mSeq *seq) {
	mFree(seq->def);
	mFree(seq->type_string);
	mFree(seq->code);
	switch (seq->type) {
		case FASTA_ALPH:
		case FASTA_DNA:
		case FASTA_PROT:
			mFree(seq->alph);
			if (seq->aux != NULL) {
				mFree((float*)seq->aux[AUX_DNA_GC]);
				mFree((coor_t*)  seq->aux[AUX_DNA_AMBIG]);
				mFree(seq->aux);
			}
			break;
		case FASTA_QUAL:
			mFree(seq->qual);
			break;
		case FASTQ:
			mFree(seq->alph);
			mFree(seq->qual);
			break;
		default:
			break;
	}
}

void mFreeSeqVector(mVector *vec) {
	int i;
	for (i=0; i<vec->size; i++) {
		mFreeSeq((mSeq*)vec->elem[i]);
		mFree(vec->elem[i]);
	}
	mFreeVector(vec);
}

/******************************
 *  Seq reading  functions    *
 ******************************/

int mReadSeqHeader(FILE *stream, mSeq *seq, char header_char) {
	int     c;
	coor_t  count = 0;
	char   *def;
	coor_t  hdr_limit = 512; /* initial mem alloc */

	def  = (char*) mCalloc(hdr_limit, sizeof(char));

	while (isspace(c = fgetc(stream)));
	if (header_char == (char) c) { /* Header line found */
		while ( (c = fgetc(stream)) != EOF && c != '\n') { 
			/* Process header */
			def[count++] = (char) c;
			if (count == hdr_limit) {
				hdr_limit *= 2;
				def = (char*) mRealloc(def, hdr_limit*sizeof(char));
			}
		}
		def[count] = '\0';
	} else {
		mDie("Not a valid %s file!", seq->type_string);
	}

	seq->def = def;
	return 1;
}

/* returns END_OF_STREAM    - if there is no more fasta entry left,
 *         HAS_MORE_ENTRIES - if there is more fasta entries left
 */
int mReadSeqUntil(FILE *stream, mSeq *seq, char header_char) {
	coor_t  seq_limit = 512; /* initial mem alloc */
	int     status = END_OF_STREAM;
	int     c;
	coor_t  count = 0;

	char *alph = (char*) mCalloc(seq_limit, sizeof(char));

	for(;;){
		c = fgetc(stream);
		if (c == EOF || c == header_char) {
			/* next record */
			break;
		}
		if ((c > 64 && c < 91) || (c > 96 && c < 123) || c == '*' || c == '-' || c == '.') {
			if (count == seq_limit) {
				seq_limit *= 2;
				alph = (char*) mRealloc(alph, seq_limit*sizeof(char));
			}
			alph[count] = c;
			count++;
		} else if (isspace(c)) {
		} else {
			mDie("Invalid character %c (%d) found in %s", c, c, seq->def);
		}
	}
	if (c == header_char) {
		/* next record */
		ungetc(c, stream);
		status = HAS_MORE_ENTRIES;
	}
	seq->alph   = alph;
	seq->length = count;

	return status;
}

/* returns END_OF_STREAM    - if there is no more fasta entry left,
 *         HAS_MORE_ENTRIES - if there is more fasta entries left
 */
int mReadAlphUntil(FILE *stream, mSeq *seq, char header_char) {
	coor_t  seq_limit = 512; /* initial mem alloc */
	int     status = END_OF_STREAM;
	int     c;
	coor_t  count = 0;

	char *alph = (char*) mCalloc(seq_limit, sizeof(char));

	for(;;){
		c = fgetc(stream);
		if (c == EOF || c == header_char) {
			/* next record */
			break;
		}
		if ((c > 64 && c < 91) || (c > 96 && c < 123)) {
			if (count == seq_limit) {
				seq_limit *= 2;
				alph = (char*) mRealloc(alph, seq_limit*sizeof(char));
			}
			alph[count] = c;
			count++;
/*
		} else if (isspace(c)) {
		} else {
			mDie("Invalid character %c (%d) found in %s", c, c, seq->def);
*/
		}
	}
	if (c == header_char) {
		/* next record */
		ungetc(c, stream);
		status = HAS_MORE_ENTRIES;
	}
	seq->alph   = alph;
	seq->length = count;

	return status;
}

/* returns END_OF_STREAM    - if there is no more fasta entry left,
 *         HAS_MORE_ENTRIES - if there is more fasta entries left
 */
int mReadLineAsChars(FILE *stream, mSeq *seq, char header_char) {
	coor_t  seq_limit = 512; /* initial mem alloc */
	int     status = END_OF_STREAM;
	int     c;
	coor_t  count = 0;

	char *alph = (char*) mCalloc(seq_limit, sizeof(char));

	for(;;){
		c = fgetc(stream);
		if (c == EOF || c == header_char) {
			/* next record */
			break;
		}
		if ((c > 64 && c < 91) || (c > 96 && c < 123) || c == '*' || c == '-' || c == '.') {
			if (count == seq_limit) {
				seq_limit *= 2;
				alph = (char*) mRealloc(alph, seq_limit*sizeof(char));
			}
			alph[count] = c;
			count++;
		} else if (isspace(c)) {
		} else {
			mDie("Invalid character %c (%d) found in %s", c, c, seq->def);
		}
	}
	if (c == header_char) {
		/* next record */
		ungetc(c, stream);
		status = HAS_MORE_ENTRIES;
	}
	seq->alph    = alph;
	seq->length = count;

	return status;
}

/* returns END_OF_STREAM    - if there is no more fasta entry left,
 *         HAS_MORE_ENTRIES - if there is more fasta entries left
 */
int mReadIntegersUntil(FILE *stream, mSeq *seq, char header_char) {
	coor_t  seq_limit = 512; /* initial mem alloc */
	int     status = END_OF_STREAM;
	int     c;
	coor_t  count = 0;

	char *alph = (char*) mCalloc(seq_limit, sizeof(char));
	for(;;){
		int val;
		if (fscanf(stream, "%d", &val) != 1) {
			c = fgetc(stream);
			if (c == EOF || c == '>') {
				/* next record */
				break;
			}
		}
		alph[count] = val;
		if (count == seq_limit) {
			seq_limit *= 2;
			alph = (char*) mRealloc(alph, seq_limit*sizeof(char));
		}
		count++;
	}
	if (c == header_char) {
		/* next record */
		ungetc(c, stream);
		status = HAS_MORE_ENTRIES;
	}
	seq->alph   = alph;
	seq->length = count;

	return status;
}

int mReadFastaAlph(FILE *stream, mSeq *seq) {
	mReadSeqHeader(stream, seq, '>');
	return mReadAlphUntil(stream, seq, '>');
}

int mReadFastaSeq(FILE *stream, mSeq *seq) {
	mReadSeqHeader(stream, seq, '>');
	return mReadSeqUntil(stream, seq, '>');
}

int mReadFastaQual(FILE *stream, mSeq *seq) {
	mReadSeqHeader(stream, seq, '>');
	return mReadIntegersUntil(stream, seq, '>');
}

int mReadFastq(FILE *stream, mSeq *seq) {
	int c;
	int seq_limit = 1024;
	coor_t length;
	char *line = (char*) mMalloc(seq_limit*sizeof(char));
	char *alph;
	char *qual;

	/* Read alph followed by @header */

	mReadSeqHeader(stream, seq, '@');
	if (fgets(line, seq_limit, stream) == NULL)
		mDie("Fatal error: encountered EOF while reading sequence");
	length = strlen(line)-1; /* last character is new-line! */
	alph = (char*) mMalloc(length*sizeof(char));
	memcpy(alph, line, length);
	seq->length = length;

	/* Read +header */

	if (fgets(line, seq_limit, stream) == NULL)
		mDie("Fatal error: encountered EOF while reading sequence");

	/* Read qual */

	if (fgets(line, seq_limit, stream) == NULL)
		mDie("Fatal error: encountered EOF while reading sequence");
	length = strlen(line)-1; /* last character is new-line! */
	if (length != seq->length)
		mDie("Fatal error: read and qual lengths do not match for %s!", seq->def);
	qual = (char*) mMalloc(length*sizeof(char));
	memcpy(qual, line, length);

	seq->alph = alph;
	seq->qual = qual;

	mFree(line);

	/* What's the status here? */

	c = fgetc(stream);
	if (c == EOF)
		return END_OF_STREAM;
	else {
		ungetc(c, stream);
		return HAS_MORE_ENTRIES;
	}
}

int mReadSeq(FILE *stream, mSeq *seq) {
	int status = END_OF_STREAM;
	mInitSeq(seq);
	switch (seq->type) {
		case FASTA_DNA:
			status = mReadFastaSeq(stream, seq);
			/*mDnaCalculateStats(dna);*/
			break;
		case FASTA_PROT:
			status = mReadFastaSeq(stream, seq);
			break;
		case FASTA_QUAL:
			status = mReadFastaQual(stream, seq);
			break;
		case FASTA_ALPH:
			status = mReadFastaAlph(stream, seq);
			break;
		case FASTQ:
			status = mReadFastq(stream, seq);
			break;
		default:
			break;
	}
	return status;
}

int mReadMultipleSeq(FILE *stream, mVector *multi) {
	int status;
	mSeq*  seq   = (mSeq*)  mCalloc(1, sizeof(mSeq));
	while((status = mReadSeq(stream, seq))) {
		mPushVector(multi, seq);
		if (status == END_OF_STREAM) { /* Last entry */
			break;
		}
		seq = (mSeq*) mCalloc(1, sizeof(mSeq));
	}
	return multi->size;
}

/******************************
 *  Seq writing  functions    *
 ******************************/

/* This is defined in the header. Copied just for reference */

void mWriteHeader(FILE *stream, mSeq *seq, char header_char) { fprintf(stream, "%c%s\n", header_char, seq->def); }

void mWriteFixedLengthChar(FILE *stream, coor_t length, char *string, int columns) {
	coor_t        i;
	int           full_lines = length - length%columns;
	for (i=0; i < full_lines; i+=columns) {
		fprintf(stream, "%.*s", columns, string+i);
		fprintf(stream, "\n");
	}
	for (i=full_lines; i < length; i++) {
		fprintf(stream, "%c", string[i]);
	}
	if (length%columns != 0) fprintf(stream, "\n");
}

/* Use this when writing quality values read in as integers */

void mWriteFixedWordsInt(FILE *stream, coor_t length, char *string, int words) {
	coor_t        i;
	for (i=0; i < length; i++) {
		fprintf(stream, "%d", string[i]);
		if (i%words == words-1) {
			fprintf(stream, "\n");
		} else {
			fprintf(stream, " ");
		}
	}
	if (length%words != 0) fprintf(stream, "\n");
}

/* Use this when writing quality values read in as string part of fastq *
 * Converted from Sanger scale text to Sanger scale quality score       */

void mWriteFixedWordsChar(FILE *stream, coor_t length, char *string, int words) {
	coor_t        i;
	for (i=0; i < length; i++) {
		fprintf(stream, "%d", string[i]-33);
		if (i%words == words-1) {
			fprintf(stream, "\n");
		} else {
			fprintf(stream, " ");
		}
	}
	if (length%words != 0) fprintf(stream, "\n");
}

void mWriteSeqN(FILE *stream, mSeq *seq, int columns) {
	if (columns == 0) columns = seq->window;
	switch (seq->type) {
		case FASTA_ALPH:
		case FASTA_DNA:
		case FASTA_PROT:
			mWriteHeader(stream, seq, '>');
			mWriteFixedLengthChar(stream, seq->length, seq->alph, columns);
			break;
		case FASTA_QUAL:
			mWriteHeader(stream, seq, '>');
			mWriteFixedWordsInt(stream, seq->length, seq->alph, columns);
			break;
		case FASTQ:
			mWriteHeader(stream, seq, '@');
			fprintf(stream, "%.*s\n", (int) seq->length, seq->alph);
			fprintf(stream, "+\n");
			fprintf(stream, "%.*s\n", (int) seq->length, seq->qual);
			break;
	}
}

void mWriteSeq(FILE *stream, mSeq *seq) {
	int columns = seq->window;
	mWriteSeqN(stream, seq, columns);
}

void mWriteMultipleSeq(FILE *stream, mVector *multi) {
	int i;
	for (i=0; i<multi->size; i++) {
		mWriteSeq(stream, (mSeq*)multi->elem[i]);
	}
}

/******************************
 *  Other generic utilities   *
 ******************************/

void mProcessSeqDef(mSeq *seq) {
	char *def    = seq->def;
	int   length = strlen(def);
	int   i;
	for (i=0; i<length; i++) {
		if (isspace(def[i])) {
			def[i] = '\0';
		}
	}
}

mSeq* mSubSeq(mSeq *seq, coor_t offset, coor_t length) {
	mSeq *sub;

	if (offset >= seq->length) return NULL;

	sub         = (mSeq*) mCalloc(1, sizeof(mSeq));
	sub->alph   = seq->alph + (int)offset;
	if (seq->qual != NULL) {
		sub->qual   = seq->qual + (int)offset;
	} else {
		sub->qual = NULL;
	}
	sub->def    = seq->def;
	sub->type   = seq->type;
	sub->window = seq->window;
	sub->type_string = seq->type_string;
	sub->length = length;
	if (offset + length > seq->length) {
		sub->length = seq->length - offset;
	}

	/* Don't calculate GC percentage. A call to mDnaGetGC will calculate it. */

	return sub;
}

void mReverseSeq(mSeq *seq) {
	coor_t i;
	coor_t length = seq->length;
	char   *alph = seq->alph;
	char   *qual = seq->qual;
	short int   *code = seq->code;
	char    tmp;
	int     qtmp;
	short int ctmp;
	for (i=0; i<length/2; i++) {
		tmp = alph[length-i-1];
		alph[length-i-1] = alph[i];
		alph[i] = tmp;
		if (qual != NULL) {
			qtmp = qual[length-i-1];
			qual[length-i-1] = qual[i];
			qual[i] = qtmp;
		}
		if (code != NULL) {
			ctmp = code[length-i-1];
			code[length-i-1] = code[i];
			code[i] = ctmp;
		}
	}
}

void mReverseComplementSeq(mSeq *seq) {
	if (seq->type == FASTA_DNA || seq->type == FASTQ) {
		mReverseComplementDNA(seq);
	} else {
		mReverseSeq(seq);
	}
}

/******************************
 *  DNA specific functions    *
 ******************************/

void mDnaCalculateStats(mSeq *dna) {
	coor_t  i;
	coor_t  at=0, gc=0;
	coor_t  length    = dna->length;
	int     aux_count = 3;
	void  **aux       = (void**)  mCalloc(aux_count, sizeof(void*));
	float  *fgc       = (float*)  mCalloc(1, sizeof(float));
	coor_t *ambiguous = (coor_t*) mCalloc(1, sizeof(coor_t));
	char   *alph      = dna->alph;

	if (dna->code == NULL) {
		dna->code = (short int*) mCalloc(length, sizeof(short int));
	}

	for(i=0; i<length; i++){
		int  c = alph[i];
		int lc = c;
		if (lc<91) lc+=32;
		switch(lc) {
			case 'a': /* 00 */
				dna->code[i] = 0; at++;
				break;
			case 'c': /* 01 */
				dna->code[i] = 1; gc++;
				break;
			case 'g': /* 10 */
				dna->code[i] = 2; gc++;
				break;
			case 't': /* 11 */
				dna->code[i] = 3; at++;
				break;
			case 'n': /* 100 */
			default:
				dna->code[i] = 4;
				break;
		}

		/* May become 1010 */
		if (c == lc) dna->code[i] |= 8;
	}

	*fgc  = 1.0 * gc / (gc+at);
	*ambiguous = length-gc-at;
	aux[AUX_DNA_GC]    = fgc;
	aux[AUX_DNA_AMBIG] = ambiguous;
	dna->aux = aux;
}

void mComplementDNA(mSeq *dna) {
	coor_t i;
	coor_t     length = dna->length;
	char      *alph   = dna->alph;
	for (i=0; i<length; i++) {
		switch(alph[i]) {
			case 'A':
				alph[i] = 'T';
				break;
			case 'a':
				alph[i] = 't';
				break;
			case 'C':
				alph[i] = 'G';
				break;
			case 'c':
				alph[i] = 'g';
				break;
			case 'G':
				alph[i] = 'C';
				break;
			case 'g':
				alph[i] = 'c';
				break;
			case 'U':
			case 'T':
				alph[i] = 'A';
				break;
			case 'u':
			case 't':
				alph[i] = 'a';
				break;
		}
	}


	if (dna->aux != NULL) {

		/* fix GC */

		short int *code = dna->code;
		float *gc = (float*) dna->aux[AUX_DNA_GC];
		*gc = (1 - *gc);

		/* weird error: for 'N', this becomes -1 and not 4.
		 * so I added a 5 in front just in case */

		for (i=0; i<length; i++) {
			code[i] = (5 + 3 - code[i])%5; 
		}
	}
}

void mReverseComplementDNA(mSeq *dna) {
	mReverseSeq(dna);
	mComplementDNA(dna);
}

char* mDNA2Char(mSeq *dna) {
	coor_t length = dna->length;
	char *seq = (char*) mCalloc(length+1, sizeof(char));
	memcpy(seq, dna->alph, length);
	seq[length] = '\0';
	return seq;
}

float mDnaGetGC(mSeq *dna) {
	/* Calculate GC percentage, which also calculates gc */
	if (dna->aux == NULL) {
		mDnaCalculateStats(dna);
	}
	return *(float*)dna->aux[AUX_DNA_GC];
}

coor_t mDnaGetUnambiguousLength(mSeq *dna) {
	/* Calculate GC percentage, which also calculates gc */
	if (dna->aux == NULL) {
		mDnaCalculateStats(dna);
	}
	return dna->length - *(int*)dna->aux[AUX_DNA_AMBIG];
}
