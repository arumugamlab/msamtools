#include "msam.h"

void mInitGlobal() {
	global = (msam_global*) mMalloc(sizeof(msam_global));
	global->MIN_LENGTH = 0;
	global->PPT = 0;
	global->MAX_CLIP = 0;
	global->multiple_input = 0;

	/* Init all pointers to NULL */
	global->header = NULL;
	global->f_coverage = NULL;
	global->seq_touched = NULL;
	global->ui_insert_count = NULL;
	global->d_insert_count = NULL;
	global->multi_mappers = NULL;
	global->ub_target_hit = NULL;
	global->pipe_fd = NULL;
}

void mFreeGlobal() {
	mFree(global);
	global = NULL;
}

void mPrintHelp (const char *subprogram, void **argtable) {
/*
	fprintf(stdout, "\n%s: Tools you wish samtools had included!\n", PROGRAM);
	fprintf(stdout,   "           Part of the awesome maui library (%s)\n", PACKAGE_NAME);
	fprintf(stdout,   "           Smash and libmaui: Making your life easy, step by step (TM)\n");
	if (strlen(BUILD) > 0)
		fprintf(stdout,   "Version  : %s %s build %s\n\n", PACKAGE_NAME, PACKAGE_VERSION, BUILD);
	else
		fprintf(stdout,   "Version  : %s %s\n\n", PACKAGE_NAME, PACKAGE_VERSION);
*/
	fprintf(stdout, "Usage:\n------\n\n%s %s", PROGRAM, subprogram);
	arg_print_syntax(stdout, argtable, "\n");
	fprintf(stdout,       "\nGeneral options:\n"
				"----------------\n\n"
				"These options specify the input/output formats of BAM/SAM files \n"
				"(same meaning as in 'samtools view'):\n");
	arg_print_glossary(stdout, argtable, "  %-25s %s\n");
}

void mPrintCommandLine(FILE *output, int argc, char *argv[]) {
	int i;
	if (strlen(BUILD) > 0)
		fprintf(output, "# %s version %s build %s\n", PACKAGE_NAME, PACKAGE_VERSION, BUILD);
	else
		fprintf(output, "# %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
	fprintf(output, "# Command: msamtools");
	for (i=0; i<argc; i++) fprintf(output, " %s", argv[i]);
	fprintf(output, "\n");
}

void mMultipleFileError(const char *subprogram, void **argtable) {
/*
	fprintf(stderr, "-t required when using multiple input files!\n");
*/
	fprintf(stderr, "Multiple input files not supported in %s.\n", subprogram);
	fprintf(stderr, "Use 'samtools merge' to combine BAM/SAM files.\n");
	mPrintHelp(subprogram, argtable);
	mQuit("");
}

samfile_t* mOpenSamFile(const char *filename, const char *inmode, char *headerfile) {
	samfile_t *input;

	/* open input file */
	/* since you pass headerfile, if there is a headerfile, it will be used */

	input = samopen(filename, inmode, headerfile);
	if (input == NULL) 
		mDie("Cannot open %s for reading", filename);

	return input;
}

int mHeaderCheck(int count, const char* filenames[], const char *inmode) {
	int i, j;
	samfile_t    *sam1    = mOpenSamFile(filenames[0], inmode, NULL);
	bam_header_t *header1 = sam1->header;
	if (header1 == 0)
		mDie("No header found in %s. Please fix your bamfile!", filenames[0]);
	for (i=1; i<count; i++) {
		samfile_t *sam2 = mOpenSamFile(filenames[i], inmode, NULL);
		bam_header_t *header2 = sam2->header;
		if (header2 == 0)
			mDie("No header found in %s. Please fix your bamfile!", filenames[i]);
		if (header1->n_targets != header2->n_targets) {
			fprintf(stderr, "%s vs %s:\n", filenames[0], filenames[i]);
			fprintf(stderr, "Sequence count in BAM header differs:\n");
			fprintf(stderr, "%d vs %d\n", header1->n_targets, header2->n_targets);
			mDie("@SQ count mismatch");
		}
		for (j=1; j<header1->n_targets; j++) {
			if (strcmp(header1->target_name[j], header2->target_name[j]) != 0) {
				fprintf(stderr, "%s vs %s:\n", filenames[0], filenames[i]);
				fprintf(stderr, "Sequence name for target %d in BAM header differs:\n", j);
				fprintf(stderr, "%s vs %s\n", header1->target_name[j], header2->target_name[j]);
				mDie("@SQ name mismatch");
			}
		}
		samclose(sam2);
	}
	samclose(sam1);
	return 1;
}

zoeHash mReadIntegerHash(const char *filename) {
	char     line[2048];
	char     key[128];
	int     *value;
	FILE    *stream;
	zoeHash  hash;

	stream = fopen(filename, "r");
	if (stream == NULL)
		mDie("Cannot open file: %s\n", filename);

	value = (int*) mCalloc(1, sizeof(int));;
	hash = zoeNewHash();
	while (fgets(line, sizeof(line), stream) != NULL) {
		if (sscanf(line, "%s\t%d", key, value) == 2) {
			zoeSetHash(hash, key, value);
		}
		value = (int*) mCalloc(1, sizeof(int));
	}
	fclose(stream);
	mFree(value);
	return hash;
}

void mEmptyZoeHash(zoeHash hash) {
	int    i;
	zoeVec vals;
	if (hash == NULL) 
		return;
	vals = zoeValsOfHash(hash); 
	for (i=0; i<vals->size; i++) {
		mFree(vals->elem[i]);
	}
	zoeDeleteVec(vals);
}

/***
 * Prepares a compressor that waits for you to write to it
 * and writes the compressed output to the given filename
 */

FILE* mInitCompressedOutputStream(const char *filename) {
	pid_t  child;
	int   *pipe_fd = (int*) mMalloc(2*sizeof(int));
	FILE  *compressor = NULL;

	if (pipe(pipe_fd) == -1) {
		perror("Error creating compressor using pipe()");
		exit(EXIT_FAILURE);
	}

	if ((child = fork()) == -1) {
		perror("Error creating compressor child using fork()");
		exit(EXIT_FAILURE);
	}

	if (child == 0) { /* Child will read in 'reader' and exit when it gets EOF */
		/* Compress the pipe input stream into outstrem */
		FILE *outstream = fopen(filename, "w");
		FILE *reader    = fdopen(pipe_fd[0], "r");

		close(pipe_fd[1]);
		def(reader, outstream);
		fclose(reader);
		fclose(outstream);

		mFree(pipe_fd);
		mFreeGlobal();

		exit(EXIT_SUCCESS);
	} 

	/* Parent just moves on */
	/* Write coverage info into pipe output stream for compressor child */

	close(pipe_fd[0]);
	compressor = fdopen(pipe_fd[1], "w");
	global->pipe_fd = pipe_fd;
	return compressor;
}

FILE* mInitOutputStream(const char *filename, int gzip) {
	FILE *output;
	if (gzip) {
		output = mInitCompressedOutputStream(filename);
	} else {
		output = fopen(filename, "w");
	}
	return output;
}

void mFreeCompressedOutputStream(FILE *compressor) {
	int *pipe_fd = global->pipe_fd;
	mFree(pipe_fd);
	fclose(compressor);
}

void mFreeOutputStream(FILE *output, int gzip) {
	if (gzip) {
		mFreeCompressedOutputStream(output);
	} else {
		fclose(output);
	}
}
