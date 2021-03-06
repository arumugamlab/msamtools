#include "msam.h"

/* This main entry file was rewritten after being inspired by the samtools codebase. */

static int usage(FILE *out)
{
	fprintf(out, "\n%s: Metagenomics-related extension to samtools.\n", PROGRAM);

	if (strlen(BUILD) > 0)
		fprintf(out, "%*s: %s %s build %s\n\n", (int) strlen(PROGRAM), "Version", PACKAGE_NAME, PACKAGE_VERSION, BUILD);
	else
		fprintf(out, "%*s: %s\n\n", (int) strlen(PROGRAM), "Version", PACKAGE_VERSION);

	fprintf(out, "Usage:\n------\n\n%s <command> [options]\n\n", PROGRAM);
	fprintf(out, "Command:  filter         filter alignments based on alignment statistics\n");
	fprintf(out, "          profile        estimate relative abundance profile of reference sequences in bam file\n");
	fprintf(out, "          coverage       estimate per-base read coverage of each reference sequence\n");
	fprintf(out, "          summary        summarize alignment statistics per read in a table format\n");
	fprintf(out, "\n");
	return 1;
}

int main(int argc, char* argv[]) {

	if (argc < 2) return usage(stderr);

	     if (strcmp(argv[1], "filter"    ) == 0) return msam_filter_main  (argc-1, argv+1);
	else if (strcmp(argv[1], "profile"   ) == 0) return msam_profile_main (argc-1, argv+1);
	else if (strcmp(argv[1], "coverage"  ) == 0) return msam_coverage_main(argc-1, argv+1);
	else if (strcmp(argv[1], "summary"   ) == 0) return msam_summary_main (argc-1, argv+1);
	else if (strcmp(argv[1], "help"      ) == 0) return usage(stdout);
	else {
		fprintf(stderr, "[msamtools] unrecognized command '%s'\n", argv[1]);
		usage(stderr);
		return 1;
	}

	return 0;
}
