# Version info
include versions.txt

# Git provenance
BUILT_SOURCES = msam_version.h
DISTCLEANFILES = msam_version.h

.PHONY: FORCE
msam_version.h: FORCE
	@set -e; \
	if git -C "$(srcdir)" rev-parse --git-dir >/dev/null 2>&1; then \
		commit=`git -C "$(srcdir)" rev-parse --short=12 HEAD`; \
		printf '%s\n' \
			'#ifndef MSAM_VERSION_H' \
			'#define MSAM_VERSION_H' \
			"#define MSAM_GIT_COMMIT \"$$commit\"" \
			'#endif' > $@.tmp; \
	elif test -f "$(srcdir)/msam_version.h"; then \
		cat "$(srcdir)/msam_version.h" > $@.tmp; \
	else \
		printf '%s\n' \
			'#ifndef MSAM_VERSION_H' \
			'#define MSAM_VERSION_H' \
			'#define MSAM_GIT_COMMIT "unknown"' \
			'#endif' > $@.tmp; \
	fi; \
	if test -f $@ && cmp -s $@.tmp $@; then \
		rm -f $@.tmp; \
	else \
		mv -f $@.tmp $@; \
	fi

FORCE:

# We need htslib dependencies to be built first
SUBDIRS = deps

# Fast deterministic correctness tests run by `make check`.
TESTS = tests/test_smoke.sh \
        tests/test_profile.sh \
        tests/test_qname_order.sh
TESTS_ENVIRONMENT = MSAMTOOLS='$(abs_builddir)/msamtools'

# These are the installed programs
mymaui_PROGRAMS = msamtools

# We will also distribute a README.md markdown file
EXTRA_DIST = README.md \
             LICENSE \
             zoe-LICENSE \
             THIRD_PARTY_NOTICES.md \
             tests/test_functions.sh \
             tests/test_smoke.sh \
             tests/test_profile.sh \
             tests/test_qname_order.sh \
             tests/fixtures/profile.sam \
             tests/fixtures/qname_coordinate.sam \
             tests/fixtures/qname_reopened.sam

# Install msamtools license information
licensedir = $(datarootdir)/licenses/msamtools
license_DATA = LICENSE \
               zoe-LICENSE \
               THIRD_PARTY_NOTICES.md

# Install licenses for third-party code incorporated during the build.
# HTSlib is downloaded into the build tree by the
# dependency build machinery and are therefore not part of EXTRA_DIST.
install-data-hook:
	$(MKDIR_P) "$(DESTDIR)$(licensedir)"
	$(INSTALL_DATA) \
		"$(builddir)/deps/htslib/htslib-$(htslib_version)/LICENSE" \
		"$(DESTDIR)$(licensedir)/htslib-LICENSE"

uninstall-hook:
	rm -f \
		"$(DESTDIR)$(licensedir)/htslib-LICENSE"
	-rmdir "$(DESTDIR)$(licensedir)"

dist-hook:
	@git -C "$(srcdir)" diff --quiet && \
	 git -C "$(srcdir)" diff --cached --quiet || { \
		echo "ERROR: refusing to create distribution from a dirty Git tree"; \
		exit 1; \
	}
	@commit=`git -C "$(srcdir)" rev-parse --short=12 HEAD` || exit 1; \
	printf '%s\n' \
		'#ifndef MSAM_VERSION_H' \
		'#define MSAM_VERSION_H' \
		"#define MSAM_GIT_COMMIT \"$$commit\"" \
		'#endif' > "$(distdir)/msam_version.h"

# Where to install the programs
mymauidir = $(bindir)

# maui library will be built and not distributed nor installed, 
# but its sources should be distributed

noinst_LIBRARIES = libmaui.a
dist_libmaui_a_SOURCES = $(srcdir)/mCommon.c \
                         $(srcdir)/mCommon.h \
                         $(srcdir)/mMatrix.c \
                         $(srcdir)/mMatrix.h \
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

ARFLAGS = cr

AM_CFLAGS = $(MYCFLAGS) \
            -DPROGRAM=\"msamtools\" \
            -DHTSLIB_VERSION=\"$(htslib_version)\" \
            -I$(builddir)/deps/htslib/htslib-$(htslib_version)

LDADD = $(builddir)/libmaui.a \
        $(builddir)/deps/htslib/htslib-$(htslib_version)/libhts.a

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
