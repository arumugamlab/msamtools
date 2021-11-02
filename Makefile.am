# Version info
include versions.txt

# We need samtools dependencies to be built first
SUBDIRS = deps

# These are the installed programs
mymaui_PROGRAMS = msamtools

# We will also distribute a README.md markdown file
EXTRA_DIST = README.md LICENSE

# Where to install the programs
mymauidir = $(bindir)

# maui library will be built and not distributed nor installed, 
# but its sources should be distributed

noinst_LIBRARIES = libmaui.a
dist_libmaui_a_SOURCES = $(srcdir)/mCommon.c \
                         $(srcdir)/mCommon.h \
                         $(srcdir)/mMatrix.c \
                         $(srcdir)/mMatrix.h \
                         $(srcdir)/mCompress.c \
                         $(srcdir)/mCompress.h \
                         $(srcdir)/zoeTools.c \
                         $(srcdir)/zoeTools.h \
                         $(srcdir)/mBamVector.h \
                         $(srcdir)/mBamVector.c \
                         $(srcdir)/msam.h \
                         $(srcdir)/msam_coverage.c \
                         $(srcdir)/msam_filter.c \
                         $(srcdir)/msam_helper.c \
                         $(srcdir)/msam_profile.c \
                         $(srcdir)/msam_summary.c

# Specify CFLAGS and LDFLAGS

AM_CFLAGS = $(MYCFLAGS) \
            -DPROGRAM=\"msamtools\" \
            -DHTSLIB_VERSION=\"$(samtools_version)\" \
            -I$(builddir)/deps/samtools/samtools-$(samtools_version)/htslib-$(samtools_version) \
            -I$(builddir)/deps/samtools/samtools-$(samtools_version)

LDADD = $(builddir)/libmaui.a \
        $(builddir)/deps/samtools/samtools-$(samtools_version)/libbam.a \
        $(builddir)/deps/samtools/samtools-$(samtools_version)/htslib-$(samtools_version)/libhts.a

clean-local:
	for file in libmaui.a; do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done

maintainer-clean-local:
	@echo 'This command is intended for maintainers to use; it'
	@echo 'deletes files that may need special tools to rebuild.'
	for file in Makefile.in; do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
