#include <limits.h>
#include "mSequence.h"
#include "zoeTools.h"
#include "argtable2.h"

int main(int argc, char* argv[]) {
	mSeq   *seq;
	FILE   *stream, *list_stream, *out;
	char    line[LINE_MAX];
	int     status;
	int     i;
	int     exclude;
	int     pair;
	int     type;
	zoeHash keep;
	zoeTVec keys;
	char   *key;

	int              argcount = 0;
	int              nerrors;
	void           **argtable;

	struct arg_str  *seq_type;
	struct arg_str  *in_file;
	struct arg_str  *out_file;
	struct arg_str  *list_file;
	struct arg_lit  *arg_gzip;
	struct arg_lit  *arg_exclude;
	struct arg_lit  *arg_pair;
	struct arg_lit  *help;
	struct arg_end  *end;

	seq_type            = arg_str1("t", "type",   "<string>",       "input sequence type or format (fasta|qual|fastq|alpha)");
	in_file             = arg_str1("i", "input",  "<file>",         "input fasta file");
	out_file            = arg_str1("o", "output", "<file>",         "output fasta file");
	list_file           = arg_str1("l", "list",   "<file>",         "file containing list of fasta identifiers");
	arg_exclude         = arg_lit0("v", "exclude",                  "exclude sequences in this list (default is false)");
	arg_pair            = arg_lit0("p", "paired",                   "get both reads from a pair corresponding to the entry; needs pairs to be marked with /1 and /2 (default is false)");
	arg_gzip            = arg_lit0("z", "gzip",                     "compressed input/output using gzip (default is false)");
	help                = arg_lit0("h", "help",                     "print this help and exit");
	end                 = arg_end(20); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */


	argtable          = (void**) mMalloc(9*sizeof(void*));
	argtable[argcount++] = seq_type;
	argtable[argcount++] = in_file;
	argtable[argcount++] = out_file;
	argtable[argcount++] = list_file;
	argtable[argcount++] = arg_exclude;
	argtable[argcount++] = arg_pair;
	argtable[argcount++] = arg_gzip;
	argtable[argcount++] = help;
	argtable[argcount++] = end;

	if (arg_nullcheck(argtable) != 0) {
		mDie("insufficient memory");
	}
	nerrors = arg_parse(argc, argv, argtable);

	if (help->count > 0) {
		fprintf(stdout, "filterSeq: get a subset of sequences from a single sequence file.\n\n");
		fprintf(stdout, "Usage: filterSeq");
		arg_print_syntax(stdout, argtable, "\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		mQuit("");
	}

	if (nerrors > 0) {
		arg_print_errors(stderr, end, "filterSeq");
		fprintf(stderr, "try using -h\n");
		mQuit("");
	}

	exclude = (arg_exclude->count > 0);
	pair    = (arg_pair->count > 0);

	if (strcmp(seq_type->sval[0], "fasta") == 0) {
		type = FASTA_DNA;
	} else if (strcmp(seq_type->sval[0], "qual") == 0) {
		type = FASTA_QUAL;
	} else if (strcmp(seq_type->sval[0], "alpha") == 0) {
		type = FASTA_ALPH;
	} else if (strcmp(seq_type->sval[0], "fastq") == 0) {
		type = FASTQ;
	} else {
		type = UNKNOWN;
		mDie("Unknown type %s", seq_type->sval[0]);
	}
	
	/* Make a hash of sequence names from the list file */
	list_stream = mSafeOpenFile(list_file->sval[0], "r", 0);
	keep = zoeNewHash();
	key  = (char*) mMalloc(LINE_MAX*sizeof(char));
	while (fgets(line, LINE_MAX, list_stream) != NULL) {
		char  def[LINE_MAX];
		int  *code;
		if (sscanf(line, "%s", def) != 1) {
			mDie("LIST LINE ERROR");
		}

		code = (int*) mMalloc(sizeof(int));
		*code = 1;

		/* If pair, then check if there is template. Otherwise, just use def */
		if (pair && mGetIlluminaTemplate(def, key)) {
			zoeSetHash(keep, key, code);
		} else {
			zoeSetHash(keep, def, code);
		}
	}
	mSafeCloseFile(list_stream, 0);

	/* Process the input files one by one */
	out = mSafeOpenFile(out_file->sval[0], "w", arg_gzip->count > 0);
	seq = (mSeq*) mMalloc(sizeof(mSeq));
	seq->type = type;
	for (i=0; i<in_file->count; i++) {
		stream = mSafeOpenFile(in_file->sval[i], "r", arg_gzip->count > 0);
		while ((status=mReadSeq(stream, seq))) {
			int  *ptr;
			if (pair && mGetIlluminaTemplate(seq->def, key)) {
				ptr = (int*)zoeGetHash(keep, key);
			} else {
				mGetFirstWord(seq->def, key);
				ptr = (int*)zoeGetHash(keep, key);
			}
			if ((ptr != NULL) != exclude) { 
				mWriteSeq(out, seq);
			} else {
			}
			mFreeSeq(seq);
			if (status==END_OF_STREAM) break;
		}
		mSafeCloseFile(stream, arg_gzip->count > 0);
	}
	mSafeCloseFile(out, arg_gzip->count > 0);

	/* Free memory etc */
	keys = zoeKeysOfHash(keep);
	for (i=0; i<keys->size; i++) {
		mFree(zoeGetHash(keep, keys->elem[i]));
	}
	zoeDeleteTVec(keys);
	zoeDeleteHash(keep);
	mFree(seq);
	mFree(key);

	exit(0);
}
