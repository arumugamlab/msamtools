#ifndef _M_SEQ_H_
#define _M_SEQ_H_

#include "mCommon.h"
#include "mDefinitions.h"

#define FASTA_DNA  (1)
#define FASTA_QUAL (2)
#define FASTA_PROT (3)
#define FASTQ      (4)
#define FASTA_ALPH (5)
#define UNKNOWN    (9)

#define AUX_DNA_GC     (0)
#define AUX_DNA_AMBIG  (1)

#define END_OF_STREAM    (21)
#define HAS_MORE_ENTRIES (22)

#define SEQ_COLUMNS (80)
#define QUAL_WORDS  (17)

struct mSeq {
	coor_t      length;
	char       *def;
	short int   type;
	char       *type_string;
	char       *alph;         /* fasta: nuc/prot sequence as char array; fastq: nuc as char array */
	char       *qual;         /* fastq: qual as char array; fasta_qual: qual as int array */
	short int  *code;         /* */
	void      **aux;          /* */
	int         aux_count;
	int         window;       /* window size (in chars/words) for printing sequence/qual */
};
typedef struct mSeq mSeq;

void  mInitSeq(mSeq *seq);
void  mFreeSeq(mSeq *seq);
int   mReadSeq(FILE *stream, mSeq *seq);
void  mWriteSeq(FILE *stream, mSeq *seq);
void  mWriteSeqN(FILE *stream, mSeq *seq, int columns);
void  mWriteMultipleSeq(FILE *stream, mVector *multi);
void  mWriteHeader(FILE *stream, mSeq *seq, char header_char);
void  mReverseSeq(mSeq *seq);
void  mReverseComplementSeq(mSeq *seq);
mSeq* mSubSeq(mSeq *seq, coor_t offset, coor_t length);

/* DNA specific */

void   mDnaCalculateStats(mSeq *dna);
float  mDnaGetGC(mSeq *dna);
coor_t mDnaGetUnambiguousLength(mSeq *dna);
void   mReverseComplementDNA(mSeq *dna);

/* Qual specific */

#define mFastqInt2Char(x) (33+x)
#define mFastqChar2Int(y) (y-33)
void mWriteFixedWordsChar(FILE *stream, coor_t length, char *string, int words);
void mWriteFixedLengthChar(FILE *stream, coor_t length, char *string, int columns);

#endif
