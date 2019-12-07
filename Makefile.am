## Process this file with automake to produce Makefile.in

include versions.txt

SUBDIRS = deps

# installed programs

mymaui_PROGRAMS = seqUtils splitSeq filterSeq msamtools

EXTRA_DIST = README.md

# figure out where to install the programs

if MULTI_ARCH_SUPPORT
mymauidir = $(datadir)/$(host)/maui
else
mymauidir = $(datadir)/maui
endif

# not-installed maui library, but its sources should be distributed

noinst_LIBRARIES = libmaui.a

# Copy the dependencies from the deps/ directory.
# The Makefile.am inside deps/ figures out what to do to
# make sure that these files are available at the top
# level of deps/ directory. The "SUBDIRS = deps" call
# will make sure that make is run inside deps/.
# Here, you just copy it.

argtable_headers = argtable2.h

libargtable2.a $(argtable_headers):
	@cp deps/$@ .

libbam.a libhts.a:
	@cp deps/$@ .

mBamVector.c: libbam.a libhts.a

mCommon.c: libargtable2.a $(argtable_headers)

msamtools_DEPENDENCIES  = libargtable2.a $(argtable_headers) libmaui.a libbam.a libhts.a
dist_libmaui_a_SOURCES  = mCommon.c mCommon.h mDefinitions.h mMatrix.c mMatrix.h mSequence.c mSequence.h mCompress.c mCompress.h zoeTools.c zoeTools.h mBamVector.h mBamVector.c msam.h msam_coverage.c msam_filter.c msam_helper.c msam_profile.c msam_summary.c

# Specify CFLAGS and LDFLAGS

# maui library and argtable are added statically to every program

LDADD = libmaui.a libargtable2.a

# bam and hts library are added statically to msamtools

msamtools_LDADD = libmaui.a libargtable2.a libbam.a libhts.a
AM_CFLAGS = $(MYCFLAGS) -DBUILD=\"@BUILD@\" -DPROGRAM=\"msamtools\" -Ideps/samtools/samtools-$(samtools_version)/htslib-$(samtools_version) -Ideps/samtools/samtools-$(samtools_version) 

uninstall-hook:
	if [ -d $(mymauidir) ]; then rmdir $(mymauidir); fi

clean-local:
	for file in libbam.a libhts.a libargtable2.a $(argtable_headers); do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done

maintainer-clean-local:
	for file in Makefile.in $(argtable_tar_file) $(samtools_tar_file); do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
