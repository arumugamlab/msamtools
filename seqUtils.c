#include <argtable2.h>
#include "mSequence.h"

#define REVCOMP  (0)
#define BREAK    (1)
#define SEPARATE (2)
#define GC       (3)
#define CLIP     (4)
#define TRIM     (5)
#define LENGTH   (6)
#define LENGTHNOGAP   (7)
#define LENGTHFILTER  (8)
#define INTERLEAVE    (9)
#define SUBSEQ       (10)
#define RENAME       (11)

int window = 0;
#define WRITESEQWIN(stream, seq) (mWriteSeqN(stream, seq, window))

int scaf_break(int seq_type, int count, const char* files[], int gzip, int threshold) {
	mSeq *seq = (mSeq*) mMalloc(sizeof(mSeq));
	mSeq *sub;
	char def[256];
	FILE *stream;
	int i;
	int status;

	seq->type = seq_type;
	for (i=0; i<count; i++) {
		stream    = mSafeOpenFile(files[i], "r", gzip);
		while ((status=mReadSeq(stream, seq))) {
			int current=0, length = (int)seq->length, counter=0;
			int gap_stretch=0;
			int seq_length=0, seq_start, seq_end=0;

			/* Skip N's at the beginning */
			while (current<length && seq->code[current] == 4) {current++;}
			seq_start = current;

			/* Now we are guaranteed to have non-N here */

			while (current<length) {
				while (current<length && seq->code[current] != 4) {current++;}
				seq_end = current-1;
				/* Now, [seq_start,seq_end] is a sequence stretch */

				/* Now skip N's if there are any */
				while (current<length && seq->code[current] == 4) {gap_stretch++;current++;}

				if (gap_stretch >= threshold || current == length) {
					seq_length = seq_end - seq_start + 1;
					sub = mSubSeq(seq, seq_start, seq_length);
					sprintf(def, "%s_part%d", sub->def, counter++);
					sub->def = def;
					WRITESEQWIN(stdout, sub);
					seq_start = current;
				}
				gap_stretch = 0;
			}
			seq_length = current - seq_start;
			sub = mSubSeq(seq, seq_start, seq_length);
			if (sub != NULL) {
				sprintf(def, "%s_part%d", sub->def, counter++);
				sub->def = def;
				WRITESEQWIN(stdout, sub);
			}
			mFreeSeq(seq);
			if (status==END_OF_STREAM) break;
		}
		mSafeCloseFile(stream, gzip);
	}
	return 0;
}

int separate(int seq_type, int count, const char* files[], int gzip, int size) {
	mSeq *seq = (mSeq*) mMalloc(sizeof(mSeq));
	FILE *stream;
	int i;
	int status;
	FILE *out1, *out2;

	out1 = stdout;
	out2 = stderr;

	seq->type = seq_type;
	for (i=0; i<count; i++) {
		stream = mSafeOpenFile(files[i], "r", gzip);
		while((status = mReadSeq(stream, seq))) {
			if (seq->length > size)
				WRITESEQWIN(out1, seq);
			else
				WRITESEQWIN(out2, seq);
			mFreeSeq(seq);
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
		}
		mSafeCloseFile(stream, gzip);
	}
	return 0;
}

#define LIMIT 5000
int length_dist(int seq_type, int file_count, const char* files[], int gzip) {
	mSeq *seq = (mSeq*) mMalloc(sizeof(mSeq));
	FILE *stream;
	int i;
	int counts[LIMIT];
	int count = 0;
	int status;
	float cumulative = 0.0;
	int sum = 0;

	for (i=0; i<LIMIT; i++) counts[i] = 0;

	seq->type = seq_type;
	for (i=0; i<file_count; i++) {
		stream = mSafeOpenFile(files[i], "r", gzip);
		while((status = mReadSeq(stream, seq))) {
			counts[seq->length]++;
			count++;
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
			mFreeSeq(seq);
		}
		mSafeCloseFile(stream, gzip);
	}
	mFree(seq);

	for (i=0; i<LIMIT; i++) {
		float this = 1.0*counts[i]/count;
		cumulative += this;
		sum += (counts[i]*i);
		if (counts[i] > 0) printf("%d\t%d\t%.6f\t%.6f\n", i, counts[i], this, cumulative);
	}
	printf("Total\t%d\n", sum);
	exit(0);
}

/* interleave between two files and combine into one */

int interleave_and_combine(int seq_type, const char* files[], int gzip, const char *outfile) {
	mSeq *seq1 = (mSeq*) mMalloc(sizeof(mSeq));
	mSeq *seq2 = (mSeq*) mMalloc(sizeof(mSeq));
	FILE *stream1, *stream2, *outstream;
	/*int   i;*/
	int   status1, status2;
	char *template1 = (char*) mMalloc(LINE_MAX*sizeof(char));
	char *template2 = (char*) mMalloc(LINE_MAX*sizeof(char));

	seq1->type = seq_type;
	seq2->type = seq_type;
	stream1 = mSafeOpenFile(files[0], "r", gzip);
	stream2 = mSafeOpenFile(files[1], "r", gzip);
	outstream = mSafeOpenFile(outfile, "w", gzip);
	while((status1 = mReadSeq(stream1, seq1)) && (status2 = mReadSeq(stream2, seq2))) {

		if (!mGetIlluminaTemplate(seq1->def, template1))
			mDie("%s is not Illumina formatted paired end file", files[0]);
		if (!mGetIlluminaTemplate(seq2->def, template2))
			mDie("%s is not Illumina formatted paired end file", files[1]);
		if (strcmp(template1, template2) != 0) {
			mDie("Order mismatch in paired fastq: %s vs %s: %s, %s\n", seq1->def, seq2->def, files[0], files[1]);
		}

		WRITESEQWIN(outstream, seq1);
		WRITESEQWIN(outstream, seq2);

		if (status1 == END_OF_STREAM || status2 == END_OF_STREAM) { /* Last entry */
			break;
		}
		mFreeSeq(seq1);
		mFreeSeq(seq2);
	}
	mSafeCloseFile(stream1, gzip);
	mSafeCloseFile(stream2, gzip);
	mSafeCloseFile(outstream, gzip);
	mFree(seq1);
	mFree(seq2);

	exit(0);
}

int mocat_summary(int seq_type, int file_count, const char* files[], int gzip) {
	mSeq *seq = (mSeq*) mMalloc(sizeof(mSeq));
	FILE *stream;
	int i;
	int status;
	int max   = 0;
	long bases = 0;
	long reads = 0;

	seq->type = seq_type;
	for (i=0; i<file_count; i++) {
		stream = mSafeOpenFile(files[i], "r", gzip);
		while((status = mReadSeq(stream, seq))) {
			int len = seq->length;
			if (len > max) max = len;
			bases += len;
			reads++;
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
			mFreeSeq(seq);
		}
		mSafeCloseFile(stream, gzip);
	}
	mFree(seq);

	printf("%ld\t%ld\t%d\n", reads, bases, max);
	exit(0);
}

int clip_file(int seq_type, const char *seqfile, int gzip, const char *lfile, const char *rfile) {
	mSeq *seq = (mSeq*) mMalloc(sizeof(mSeq));
	mSeq *sub = (mSeq*) mMalloc(sizeof(mSeq));
	FILE *seqstream = NULL;
	FILE *lstream   = NULL;
	FILE *rstream   = NULL;
	int status;
	char lline[LINE_MAX];
	char rline[LINE_MAX];
	char ldef[LINE_MAX];
	char rdef[LINE_MAX];

	seq->type = seq_type;

	/* Open all the files */

	seqstream = mSafeOpenFile(seqfile, "r", gzip);
	/* lists are never gzipped, so send gzip=0 */
	if (lfile != NULL)
		lstream = mSafeOpenFile(lfile, "r", 0);
	if (rfile != NULL)
		rstream = mSafeOpenFile(rfile, "r", 0);

	/* Process sequences */

	if (rfile == NULL) { /* only left clip */
		int lclip;
		char *sdef = (char*) mMalloc(LINE_MAX*sizeof(char));
		while((status = mReadSeq(seqstream, seq)) && (fgets(lline, LINE_MAX, lstream) != NULL)) {
			int   len  = seq->length;
			mGetFirstWord(seq->def, sdef);

			if (sscanf(lline, "%s %d", ldef, &lclip) != 2) {
				mDie("LIST LINE ERROR");
			}
			if (lclip != -1) {
				if (strcmp(ldef, sdef) != 0) {
					mDie("Wrong order in lclip file. Expected %s, found %s", seq->def, ldef);
				}
				sub = mSubSeq(seq, lclip, len-lclip);
				WRITESEQWIN(stdout, sub);
			}
			mFreeSeq(seq);
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
		}
		mFree(sdef);
	} else if (lfile == NULL) { /* only right clip */
		int rclip;
		char *sdef = (char*) mMalloc(LINE_MAX*sizeof(char));
		while((status = mReadSeq(seqstream, seq)) && (fgets(rline, LINE_MAX, rstream) != NULL)) {
			int   len  = seq->length;
			mGetFirstWord(seq->def, sdef);

			if (sscanf(rline, "%s %d", rdef, &rclip) != 2) {
				mDie("LIST LINE ERROR");
			}
			if (rclip != -1) {
				if (strcmp(rdef, sdef) != 0) {
					mDie("Wrong order in rclip file. Expected %s, found %s", seq->def, rdef);
				}
				sub = mSubSeq(seq, 0, len-rclip);
				WRITESEQWIN(stdout, sub);
			}
			mFreeSeq(seq);
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
		}
		mFree(sdef);
	} else { /* both clips */
		int lclip, rclip;
		char *sdef = (char*) mMalloc(LINE_MAX*sizeof(char));
		while((status = mReadSeq(seqstream, seq)) && (fgets(lline, LINE_MAX, lstream) != NULL) && (fgets(rline, LINE_MAX, rstream) != NULL)) {
			int   len  = seq->length;
			mGetFirstWord(seq->def, sdef);

			if (sscanf(lline, "%s %d", ldef, &lclip) != 2) {
				mDie("LIST LINE ERROR");
			}
			if (sscanf(rline, "%s %d", rdef, &rclip) != 2) {
				mDie("LIST LINE ERROR");
			}
			if (lclip != -1 && rclip != -1) {
				if (strcmp(ldef, sdef) != 0) {
					mDie("Wrong order in lclip file. Expected %s, found %s", seq->def, ldef);
				}
				if (strcmp(ldef, sdef) != 0) {
					mDie("Wrong order in rclip file. Expected %s, found %s", seq->def, rdef);
				}
				sub = mSubSeq(seq, lclip, len-rclip-lclip);
				WRITESEQWIN(stdout, sub);
			}
			mFreeSeq(seq);
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
		}
		mFree(sdef);
	}

	mSafeCloseFile(seqstream, gzip);
	if (lfile != NULL) mSafeCloseFile(lstream, 0);
	if (rfile != NULL) mSafeCloseFile(rstream, 0);

	return 1;
}

/* generic() is reserved for any operation that is atomically operating
 * on each sequence in the file independently.
 */
int generic(int seq_type, int mode, int count, const char* files[], int gzip, int arg1, int arg2) {
	mSeq *seq = (mSeq*) mMalloc(sizeof(mSeq));
	mSeq *sub = (mSeq*) mMalloc(sizeof(mSeq));
	FILE *stream;
	int i;
	int j;
	int status;
	char *new_header;

	seq->type = seq_type;

	for (i=0; i<count; i++) {
		stream = mSafeOpenFile(files[i], "r", gzip);
		j = 0;
		while((status = mReadSeq(stream, seq))) {
			j++;
			switch(mode) {
				case RENAME:
					new_header = (char*) mCalloc(512, sizeof(char));
					sprintf(new_header, "scaffold_%d", j);
					mFree(seq->def);
					seq->def = new_header;
					WRITESEQWIN(stdout, seq);
					break;
				case REVCOMP:
					mReverseComplementSeq(seq);
					WRITESEQWIN(stdout, seq);
					break;
				case TRIM:
					sub = mSubSeq(seq, 0, (coor_t) arg1);
					WRITESEQWIN(stdout, sub);
					break;
				case BREAK:
					{
						coor_t current = 0;
						coor_t max     = seq->length;
						int    cc      = 1;
						char   def[64];
						while (current < max) {
							sub = mSubSeq(seq, current, (coor_t) arg1);
							sprintf(def, "%s_p%d_%ld-%ld", seq->def, cc, current+1, current+sub->length);
							sub->def = def;
							WRITESEQWIN(stdout, sub);
							current += arg1;
							cc++;
						}
					}
					break;
				case LENGTH:
					fprintf(stdout, "%s\t%ld\n", seq->def, seq->length);
					break;
				case LENGTHNOGAP:
					mDnaCalculateStats(seq);
					fprintf(stdout, "%s\t%ld\n", seq->def, mDnaGetUnambiguousLength(seq)); 
					break;
				case GC:
					fprintf(stdout, "%s\t%.4f\n", seq->def, mDnaGetGC(seq));
					break;
				case CLIP:
					sub = mSubSeq(seq, arg1, seq->length - arg1 - arg2);
					WRITESEQWIN(stdout, sub);
				case SUBSEQ:
					sub = mSubSeq(seq, arg1, arg2-arg1+1);
					WRITESEQWIN(stdout, sub);
					break;
				case LENGTHFILTER:
					if (seq->length >= arg1 && seq->length <= arg2)
						WRITESEQWIN(stdout, seq);
					break;
			}
			mFreeSeq(seq);
			if (status == END_OF_STREAM) { /* Last entry */
				break;
			}
		}
		mSafeCloseFile(stream, gzip);
	}
	mFree(sub);
	mFree(seq);
	return 0;
}

int print_help(void **argtable) {
	fprintf(stdout, "seqUtils: A suite of tools for sequence manipulation or summary.\n\n");
	fprintf(stdout, "Usage: seqUtils");
	arg_print_syntax(stdout, argtable, "\n");
	arg_print_glossary(stdout, argtable, "  %-25s %s\n");
	fprintf(stdout, "\nProcessing modes:\n");
	fprintf(stdout, "%10s    %s\n", "", "");
	fprintf(stdout, "%10s    %s\n", "revcomp",   "reverse complement each sequence");
	fprintf(stdout, "%10s    %s\n", "trim",      "trim each sequence to the first <size> bases");
	fprintf(stdout, "%10s    %s\n", "break",     "break sequences into fragments of size <size>");
	fprintf(stdout, "%10s    %s\n", "scafbreak", "break scaffolds at stretches of N's longer than <size>");
	fprintf(stdout, "%10s    %s\n", "separate",  "separate sequences into two files");
	fprintf(stdout, "%10s    %s\n", " ",         "  sequences longer than <size> go to STDOUT");
	fprintf(stdout, "%10s    %s\n", " ",         "  sequences shorter than or equal to <size> go to STDERR");
	fprintf(stdout, "%10s    %s\n", "subseq",    "get the subsequence starting at <left> and finishing at <right>, both 0-based coordinates");
	fprintf(stdout, "%10s    %s\n", "clip",      "clip first or last <n> bases from every read, combined with --left or --right");
	fprintf(stdout, "%10s    %s\n", "clip454",   "clip first 4 bases from the first cycle in every read");
	fprintf(stdout, "%10s    %s\n", "clipfile",  "clip bases from the file based on a clip file via following options");
	fprintf(stdout, "%10s    %s\n", "  ",        "  -lfile: clip info for clipping the 5' end");
	fprintf(stdout, "%10s    %s\n", " ",         "  -rfile: clip info for clipping the 3' end");
	fprintf(stdout, "%10s    %s\n", " ",         "  one or both of the -lfile and -rfile options must be used");
	fprintf(stdout, "%10s    %s\n", " ",         "  clip information should be on the first word of sequence header");
	fprintf(stdout, "%10s    %s\n", " ",         "  e.g., 'seq1 15' is valid but 'seq1 my seq 15' is not");
	fprintf(stdout, "%10s    %s\n", "gc",        "report GC content of each sequence");
	fprintf(stdout, "%10s    %s\n", "length",    "report length of each sequence in all files");
	fprintf(stdout, "%10s    %s\n", "lengthnogap", "report unambiguous length of each sequence in all files");
	fprintf(stdout, "%10s    %s\n", "dist",      "report length distribution of sequences in all files");
	fprintf(stdout, "%10s    %s\n", "mocat_stats", "report fastq stats for mocat");
	fprintf(stdout, "%10s    %s\n", "interleave","interleave between two files and combine into one file");
	fprintf(stdout, "%10s    %s\n", "rename",    "interleave between two files and combine into one file");
	mQuit("");
	return 0;
}

int main(int argc, char* argv[]) {
	struct arg_str  *arg_mode;
	struct arg_str  *arg_type;
	struct arg_str  *arg_out;
	struct arg_file *arg_files;
	struct arg_int  *arg_size;
	struct arg_int  *arg_min;
	struct arg_int  *arg_max;
	struct arg_int  *arg_left;
	struct arg_int  *arg_right;
	struct arg_int  *arg_window;
	struct arg_file *arg_lfile;
	struct arg_file *arg_rfile;
	struct arg_lit  *arg_gzip;
	struct arg_end  *end;
	int              nerrors;
	int              threshold;
	int              type;
	int              gzip = 0;
	int		 argcount = 0;
	struct arg_lit  *help;
	void           **argtable;

	arg_mode    = arg_str1("m",  "mode", "<string>", "revcomp|trim|subseq|clip|clip454|clipfile|break|scafbreak|separate|gc|length|dist|mocat_stats|interleave|rename");
	arg_type    = arg_str1("t",  "type", "<string>", "input sequence type or format (fasta|qual|fastq|alpha)");
	arg_gzip    = arg_lit0("z",  "gzip",               "compressed input/output using gzip (default is false)");
	arg_out     = arg_str0("o", "output", "<file>",   "output file");
	arg_size    = arg_int0("s",  "size", "<int>",    "size of output fragments for \"split\" (or) threshold for separating in \"separate\"");
	arg_min     = arg_int0(NULL, "min",  "<int>",    "min length of fragment to keep in \"lengthfilter\"");
	arg_max     = arg_int0(NULL, "max",  "<int>",    "max length of fragment to keep in \"lengthfilter\"");
	arg_left    = arg_int0(NULL, "left",  "<int>",   "5' clip length for \"-m clip\"");
	arg_right   = arg_int0(NULL, "right", "<int>",   "3' clip length for \"-m clip\"");
	arg_window  = arg_int0("w",  "window","<int>",   "number of chars/words per line (default: 80 for fasta, 500 for fastq, 17 for qual)");
	arg_lfile   = arg_file0(NULL, "lfile", "FILE",   "File with individual 5' clip length for \"-m clipfile\"");
	arg_rfile   = arg_file0(NULL, "rfile", "FILE",   "File with individual 3' clip length for \"-m clipfile\"");
	arg_files   = arg_filen(NULL, NULL,  "FILE", 0, argc+2, NULL);
	help        = arg_lit0("h",   "help",            "print this help and exit");
	end         = arg_end(20); /* this needs to be even, otherwise each element in end->parent[] crosses an 8-byte boundary */

	argtable    = (void**) mMalloc(16*sizeof(void*));
	argtable[argcount++] = arg_mode;
	argtable[argcount++] = arg_type;
	argtable[argcount++] = arg_out;
	argtable[argcount++] = arg_gzip;
	argtable[argcount++] = arg_size;
	argtable[argcount++] = arg_min;
	argtable[argcount++] = arg_max;
	argtable[argcount++] = arg_left;
	argtable[argcount++] = arg_right;
	argtable[argcount++] = arg_window;
	argtable[argcount++] = arg_lfile;
	argtable[argcount++] = arg_rfile;
	argtable[argcount++] = arg_files;
	argtable[argcount++] = help;
	argtable[argcount++] = end;

	arg_size->ival[0] = -1;

	if (arg_nullcheck(argtable) != 0) {
		mDie("insufficient memory");
	}
	nerrors = arg_parse(argc, argv, argtable);

	if (help->count > 0) {
		print_help(argtable);
		mQuit("");
	}

	if (nerrors > 0) {
		arg_print_errors(stderr, end, "seqUtils");
		fprintf(stderr, "try using -h\n");
		mQuit("");
	}

	threshold = arg_size->ival[0];

	if (strcmp(arg_type->sval[0], "fasta") == 0) {
		type = FASTA_DNA;
	} else if (strcmp(arg_type->sval[0], "qual") == 0) {
		type = FASTA_QUAL;
	} else if (strcmp(arg_type->sval[0], "alpha") == 0) {
		type = FASTA_ALPH;
	} else if (strcmp(arg_type->sval[0], "fastq") == 0) {
		type = FASTQ;
	} else {
		type = UNKNOWN;
		mDie("Unknown type %s", arg_type->sval[0]);
	}

	if (arg_gzip->count > 0) 
		gzip = 1;

	if (arg_window->count > 0) 
		window = arg_window->ival[0];

	if (strcmp(arg_mode->sval[0], "revcomp") == 0) {
		generic(type, REVCOMP, arg_files->count, arg_files->filename, gzip, 0, 0);
	} else if (strcmp(arg_mode->sval[0], "rename") == 0) {
		generic(type, RENAME, arg_files->count, arg_files->filename, gzip, 0, 0);
	} else if (strcmp(arg_mode->sval[0], "trim") == 0) {
		generic(type, TRIM, arg_files->count, arg_files->filename, gzip, threshold, 0);
	} else if (strcmp(arg_mode->sval[0], "gc") == 0) {
		generic(type, GC, arg_files->count, arg_files->filename, gzip, 0, 0);
	} else if (strcmp(arg_mode->sval[0], "subseq") == 0) {
		generic(type, SUBSEQ, arg_files->count, arg_files->filename, gzip, arg_left->ival[0], arg_right->ival[0]);
	} else if (strcmp(arg_mode->sval[0], "clip") == 0) {
		generic(type, CLIP, arg_files->count, arg_files->filename, gzip, arg_left->ival[0], arg_right->ival[0]);
	} else if (strcmp(arg_mode->sval[0], "length") == 0) {
		generic(type, LENGTH, arg_files->count, arg_files->filename, gzip, 0, 0);
	} else if (strcmp(arg_mode->sval[0], "lengthnogap") == 0) {
		generic(type, LENGTHNOGAP, arg_files->count, arg_files->filename, gzip, 0, 0);
	} else if (strcmp(arg_mode->sval[0], "lengthfilter") == 0) {
		generic(type, LENGTHFILTER, arg_files->count, arg_files->filename, gzip, arg_min->ival[0], arg_max->ival[0]);
	} else if (strcmp(arg_mode->sval[0], "break") == 0) {
		if (threshold == -1) {
			print_help(argtable);
			exit(-1);
		}
		generic(type, BREAK, arg_files->count, arg_files->filename, gzip, threshold, 0);
	} else if (strcmp(arg_mode->sval[0], "clipfile") == 0) {
		const char *lfile = NULL;
		const char *rfile = NULL;
		if (arg_files->count > 1) {
			mDie("mode=clipfile can only take one sequence file");
		}
		if (arg_lfile->count > 0) lfile=arg_lfile->filename[0];
		if (arg_rfile->count > 0) rfile=arg_rfile->filename[0];
		if (lfile == NULL && rfile == NULL) {
			mDie("mode=clipfile needs at least one clip file");
		}
		clip_file(type, arg_files->filename[0], gzip, lfile, rfile);
	} else if (strcmp(arg_mode->sval[0], "scafbreak") == 0) {
		if (threshold == -1) {
			print_help(argtable);
			exit(-1);
		}
		scaf_break(type, arg_files->count, arg_files->filename, gzip, threshold);
	} else if (strcmp(arg_mode->sval[0], "separate") == 0) {
		if (threshold == -1) {
			print_help(argtable);
			exit(-1);
		}
		separate(type, arg_files->count, arg_files->filename, gzip, threshold);
	} else if (strcmp(arg_mode->sval[0], "interleave") == 0) {
		if (arg_files->count != 2) {
			mDie("mode=interleave needs two files");
		}
		if (arg_out->count != 1) {
			mDie("mode=interleave needs one argument for --output");
		}
		interleave_and_combine(type, arg_files->filename, gzip, arg_out->sval[0]);
	} else if (strcmp(arg_mode->sval[0], "dist") == 0) {
		length_dist(type, arg_files->count, arg_files->filename, gzip);
	} else if (strcmp(arg_mode->sval[0], "mocat_stats") == 0) {
		mocat_summary(type, arg_files->count, arg_files->filename, gzip);
	}
	arg_freetable(argtable, 3);
	mFree(argtable);
	return 0;
}
