#######################################################
# Initialize autoconf/automake environment
#######################################################


AC_PREREQ([2.69])

AC_INIT(
    [msamtools],
    [1.1.3],
    [https://github.com/arumugamlab/msamtools/issues],
    [msamtools],
    [https://github.com/arumugamlab/msamtools]
)

# we need the aux_dir for installing config.guess to bin/ later.
# So store it in a variable now, and use AC_SUBST later

AC_CONFIG_AUX_DIR([build-aux])
aux_files_dir=build-aux

# automake

AM_INIT_AUTOMAKE([-Wall gnu foreign tar-pax subdir-objects])

# Compilers etc

AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([C])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

#######################################################
# Programs required for dependency download
#######################################################

AC_CHECK_PROG([GZIP], [gzip], [found], [notfound])
AS_IF(
    [test "x$GZIP" = "xnotfound"],
    [AC_MSG_ERROR([gzip not found in PATH])]
)

AC_CHECK_PROG([BZIP2], [bzip2], [found], [notfound])
AS_IF(
    [test "x$BZIP2" = "xnotfound"],
    [AC_MSG_ERROR([bzip2 not found in PATH])]
)

AC_CHECK_PROG([TAR], [tar], [found], [notfound])
AS_IF(
    [test "x$TAR" = "xnotfound"],
    [AC_MSG_ERROR([tar not found in PATH])]
)

AC_CHECK_PROGS([WGET], [wget curl], [notfound])
AS_IF(
    [test "x$WGET" = "xnotfound"],
    [AC_MSG_ERROR([wget or curl not found in PATH])],
    [test "x$WGET" = "xcurl"],
    [WGET_ARGS="--location --progress-bar -o -"],
    [WGET_ARGS="--progress=dot:mega -O -"]
)
AC_SUBST([WGET_ARGS])

AC_ARG_VAR(BUILD, [Build ID for the programs])
AC_SUBST(BUILD, [$BUILD])

#######################################################
# programs for installation
#######################################################

AC_PROG_SED
AC_PROG_GREP
AC_PROG_INSTALL

#######################################################
# C dependencies
#######################################################

# argtable

AC_CHECK_HEADERS(
    [argtable2.h],
    [],
    [AC_MSG_ERROR([argtable2 development headers not found])]
)

AC_CHECK_LIB(
    [argtable2],
    [arg_parse],
    [],
    [AC_MSG_ERROR([argtable2 library not found])]
)

# libz

AC_CHECK_HEADERS(
    [zlib.h],
    [],
    [AC_MSG_ERROR([zlib development headers not found])]
)

AC_CHECK_LIB(
    [z],
    [inflate],
    [],
    [AC_MSG_ERROR([zlib library not found])]
)

# log() from libm

AC_SEARCH_LIBS(
    [log],
    [m],
    [],
    [AC_MSG_ERROR([msamtools requires a working floating-point math library])]
)

# pthread for samtools

AC_SEARCH_LIBS(
    [pthread_create],
    [pthread],
    [],
    [AC_MSG_ERROR([msamtools requires a working pthread library])]
)

#######################################################
# Figure out host type
#######################################################

# We use AC_CANONICAL_HOST to figure it out.
# It requires that we distribute config.guess, config.sub and install.sh with the distro.

AC_CANONICAL_HOST

AC_SUBST([MYCFLAGS], ["-W -Wall -Wwrite-strings -std=gnu99 -pedantic"])

AC_SUBST([MACHINE_TYPE], [$host_cpu])
kernel_name=`uname -s`
AC_SUBST([KERNEL_NAME], [$kernel_name])

#SOFTWAREDIR="$datarootdir"
#AC_SUBST([SOFTWAREDIR], [$SOFTWAREDIR])

#######################################################
# user specified input through environmental variables
#######################################################

# substitutions

# local directories

AC_SUBST([HOST],    [generic])
AC_SUBST([AUX_DIR], [$aux_files_dir])

#######################################################
# Make outputs
#######################################################

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([deps/Makefile])
AC_CONFIG_FILES([deps/samtools/Makefile])
AS_IF([test -f "data/samtools/Makefile"], [rm data/samtools/Makefile])

AC_OUTPUT

#######################################################
# Tell the user you are done!
#######################################################

echo "
-------------------------------------------------------------------------------
msamtools was successfully configured with the following options:

   Host type          : $host
   Install directory  : $prefix

You can build msamtools executables by typing:

   make

You can install msamtools now by typing:

   make install

This will install the executables in $prefix/bin
-------------------------------------------------------------------------------
"
