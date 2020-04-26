#include <sys/stat.h>
#include <limits.h>
#include "mSequence.h"
#include <argtable2.h>

#define FILE_SIZE   (1)
#define SEQ_LENGTH  (2)
#define FASTQ2FASTA (3)

int main(int argc, char* argv[]) {
	mSeq   *seq;
	FILE   *stream;
	int     status;
	int     i;
	int     max;
	int     mode;
	int     type;
	char   *prefix;
	char   *infile;
	char    suffix[16];
	float   factor;

	int              argcount = 0;
	int              nerrors;
	void           **argtable;

	struct arg_str  *arg_type;
	struct arg_str  *arg_infile;
	struct arg_str  *arg_prefix;
	struct arg_str  *arg_mode;
	struct arg_lit  *arg_gzip;
	struct arg_int  *arg_maxlength;
	struct arg_int  *arg_size;
	struct arg_lit  *arg_need_qual;
	struct arg_lit  *help;
	struct arg_end  *end;

	arg_mode            = arg_str1("m", "mode",   "<string>", "mode to split fasta files (filesize|seqlength|fastq2fasta)\n\n"
									"\tseqlength   : Creates one file for each length of sequence.\n"
									"\t              Files are named <prefix>.<n>bp.fa or such.\n"
									"\t              All sequences longer than <max> are lumped into one\n"
									"\t              file so that you dont end up a few hundred files!\n"
									"\tfilesize    : Creates several approximately equal-sized files, \n"
									"\t              each containing at most <size> bytes.\n"
									"\t              Files are named <prefix>.<nnn>.fa or such.\n"
									"\tfastq2fasta : Reads in a fastq file and writes out fasta and optionally quality files.\n"
									"\t              Files are named <prefix>.fa and <prefix>.qual.\n"
									"\t              Use '-' as prefix to write fasta sequence to stdout.\n"
									"\t              Assumes that fastq quality scores are in Sanger scale.\n"
									"Input/output options:\n"
									"---------------------\n\n"
									"These options specify the input/output formats of fasta/fastq files:\n");
	arg_infile          = arg_str1("i", "input",  "<file>",   "input sequence file");
	arg_type            = arg_str1("t", "type",   "<string>", "input sequence type or format (fasta|qual|fastq)");
	arg_gzip            = arg_lit0("z", "gzip",               "input is compressed using gzip (default: false)");
	arg_prefix          = arg_str1("p", "prefix", "<string>", "prefix for output fasta files\n\n"
									"Mode-specific options:\n"
									"----------------------\n\n"
									"1. seqlength: Creates one file for each length of sequence.\n");
	arg_maxlength       = arg_int0(NULL,"max",    "<int>",    "sequences longer than this will be lumped together (default: 100)\n\n"
									"2. filesize : Creates several files, each containing at most <size>bp of sequence.\n");
	arg_size            = arg_int0("s", "size",   "<int>",    "number of basepairs (approx.) that will be in each file (default: 1000000)\n\n"
									"3. fastq2fasta : convert input fastq to fasta/qual files.\n");
	arg_need_qual       = arg_lit0(NULL,"qual",               "also write a quality file (default: false) "
									"General options:\n"
									"----------------\n");
	help                = arg_lit0("h", "help",               "print this help and exit");
	end                 = arg_end(10); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */

	argtable          = (void**) mCalloc(10, sizeof(void*));
	argtable[argcount++] = arg_mode;
	argtable[argcount++] = arg_infile;
	argtable[argcount++] = arg_type;
	argtable[argcount++] = arg_gzip;
	argtable[argcount++] = arg_prefix;
	argtable[argcount++] = arg_maxlength;
	argtable[argcount++] = arg_size;
	argtable[argcount++] = arg_need_qual;
	argtable[argcount++] = help;
	argtable[argcount++] = end;

	/* defaults */

	arg_maxlength->ival[0] = 100;
	arg_size->ival[0] = 1000000;

	/* parse command line */

	if (arg_nullcheck(argtable) != 0) {
		mDie("insufficient memory");
	}
	nerrors = arg_parse(argc, argv, argtable);

	if (help->count > 0) {
		fprintf(stdout, "splitSeq: Split a single sequence file into many files.\n\n");
		fprintf(stdout, "Usage: splitSeq");
		arg_print_syntax(stdout, argtable, "\n");
		fprintf(stdout, "\nChoosing the right mode:\n");
		fprintf(stdout,   "------------------------\n\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		mQuit("");
	}

	if (nerrors > 0) {
		arg_print_errors(stderr, end, "splitSeq");
		fprintf(stderr, "try using -h\n");
		mQuit("");
	}

	/* figure out which mode we run under */

	if (strcmp(arg_mode->sval[0], "filesize") == 0) {
		mode = FILE_SIZE;
	} else if (strcmp(arg_mode->sval[0], "seqlength") == 0) {
		mode = SEQ_LENGTH;
	} else if (strcmp(arg_mode->sval[0], "fastq2fasta") == 0) {
		mode = FASTQ2FASTA;
	} else {
		mode = -1;
		fprintf(stderr, "Unknown mode: %s", arg_mode->sval[0]);
		fprintf(stdout, "Usage: splitSeq");
		arg_print_syntax(stdout, argtable, "\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		mQuit("");
	}

	/* arg_prefix for the output files */

	prefix = (char*) mCalloc(256, sizeof(char));
	infile = (char*) mCalloc(256, sizeof(char));

	strcpy(prefix, arg_prefix->sval[0]);
	strcpy(infile, arg_infile->sval[0]);

	if (strcmp(arg_type->sval[0], "fasta") == 0) {
		type = FASTA_DNA;
		strcpy(suffix, "fa");
		factor = 1.01*(SEQ_COLUMNS+1)/SEQ_COLUMNS; /* chars for seq + newlines */
	} else if (strcmp(arg_type->sval[0], "qual") == 0) {
		type = FASTA_QUAL;
		strcpy(suffix, "qual");
		factor = 1.01*(3*QUAL_WORDS-1)/QUAL_WORDS; /* chars for qual + spaces + newlines */
	} else if (strcmp(arg_type->sval[0], "fastq") == 0) {
		type = FASTQ;
		strcpy(suffix, "fq");
		factor = 2.01; /* chars for seq + qual + newlines */
	} else {
		type = UNKNOWN;
		factor = 1;
		mDie("Unknown type %s", arg_type->sval[0]);
	}
	
	/* open input stream */

	stream = mSafeOpenFile(infile, "r", arg_gzip->count > 0);

	seq = (mSeq*) mCalloc(1, sizeof(mSeq));
	seq->type = type;

	if (mode == SEQ_LENGTH) {

		FILE  **out;

		/* max length to separate - length>=max will go to a file for max */

		max = arg_maxlength->ival[0];

		/* initialize all the output stream variables */

		out = (FILE**) mCalloc((max+1), sizeof(FILE*));
		for (i=0; i<=max; i++) out[i] = NULL;

		/* start reading fasta files */

		while ((status=mReadSeq(stream, seq))) {
			int length = seq->length;
			if (length > max) length = max;
			if (out[length] == NULL) {
				char name[256];
				sprintf(name, "%s.%dbp.%s", prefix, length, suffix);
				if ((out[length] = fopen(name, "w")) == NULL) {
					mDie("Cannot open output file %s for reading", name);
				}
			}
			mWriteSeq(out[length], seq);
			mFreeSeq(seq);
			if (status==END_OF_STREAM) break;
		}

		/* close all the output streams */

		for (i=0; i<=max; i++) {
			if (out[i] != NULL) {
				fclose(out[i]);
			}
		}
		mFree(out);
	} else if (mode == FILE_SIZE) {

		char name[256];
		struct stat statbuf;
		int current_size = 0;
		int counter = 0;
		off_t filesize;
		off_t size;
		int padding;
		FILE *out;

		out = (FILE*) mCalloc(1, sizeof(FILE));
		stat(infile, &statbuf);
		filesize = statbuf.st_size;
		size     = (off_t) arg_size->ival[0];
		padding  = 1+(int)log10(7.0*filesize/size); /* assuming gzip compresses 7-fold */

		/* start reading fasta files */

		sprintf(name, "%s%0*d.%s", prefix, padding, counter, suffix);
		if ((out = fopen(name, "w")) == NULL) {
			mDie("Cannot open output file %s for writing", name);
		}

		while ((status=mReadSeq(stream, seq))) {
			int length = seq->length;
			current_size += factor*length + strlen(seq->def) + 2; // factors + headers + newlines
			if (current_size > size) {
				fclose(out);
				counter++;
				sprintf(name, "%s%0*d.%s", prefix, padding, counter, suffix);
				if ((out = fopen(name, "w")) == NULL) {
					mDie("Cannot open output file %s for writing", name);
				}
				current_size = length;
			}
			mWriteSeq(out, seq);
			mFreeSeq(seq);
			if (status==END_OF_STREAM) break;
		}
	} else if (mode == FASTQ2FASTA) {
		char sname[256];
		char qname[256];
		FILE *sout;
		FILE *qout = NULL;

		if (arg_gzip->count > 0) strcpy(suffix, ".gz");
		else                     strcpy(suffix, "");

		if (strcmp(prefix, "-") == 0) {
			strcpy(sname, "-");
		} else {
			sprintf(sname, "%s.fa%s", prefix, suffix);
		}
		sout = mSafeOpenFile(sname, "w", arg_gzip->count > 0);
		if (arg_need_qual->count > 0) {
			if (strcmp(prefix, "-") == 0) {
				mDie("Error: Cannot use - as prefix with --qual", arg_type->sval[0]);
			}
			sprintf(qname, "%s.qual%s", prefix, suffix);
			qout = mSafeOpenFile(qname, "w", arg_gzip->count > 0);
		}

		while ((status=mReadSeq(stream, seq))) {

			/* Write fasta seq */
			seq->window = 500; /* line length equal to fastq line length */
			seq->type = FASTA_DNA;
			mWriteSeq(sout, seq);

			/* Write qual if needed */
			if (arg_need_qual->count > 0) {
				seq->window = QUAL_WORDS;
				seq->type = FASTA_QUAL;
				mWriteHeader(qout, seq, '>');
				mWriteFixedWordsChar(qout, seq->length, seq->qual, QUAL_WORDS);
			}

			seq->type = FASTQ;
			mFreeSeq(seq);
			if (status==END_OF_STREAM) break;
		}
		mSafeCloseFile(sout, arg_gzip->count > 0);
		if (arg_need_qual->count > 0) {
			mSafeCloseFile(qout, arg_gzip->count > 0);
		}
	}

	/* close input stream */

	if (arg_gzip->count > 0) {
		pclose(stream);
	} else {
		fclose(stream);
	}

	/* Free memory etc */

	mFree(seq);
	mFree(prefix);
	mFree(infile);

	exit(0);
}
