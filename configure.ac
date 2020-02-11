#######################################################
# Initialize autoconf/automake environment
#######################################################


AC_PREREQ(2.62)
AC_REVISION([$Revision: 1.0 $])
AC_COPYRIGHT()
AC_INIT([msamtools], [0.9], [http://www.arumugamlab.org])

# we need the aux_dir for installing config.guess to bin/ later.
# So store it in a variable now, and use AC_SUBST later

AC_CONFIG_AUX_DIR([build-aux])
aux_files_dir=build-aux

#######################################################
# Define missing things: macro definitions and such
#######################################################

# AS_VAR_APPEND

dnl autoconf versions before 2.63b don't have AS_VAR_APPEND
m4_ifdef([AS_VAR_APPEND],,[as_fn_append () { eval $[1]=\$$[1]\$[2]; }
AC_DEFUN([AS_VAR_APPEND],[as_fn_append $1 $2])])dnl

# automake

AM_INIT_AUTOMAKE([-Wall gnu foreign tar-pax])

#######################################################
#######################################################
# programs for external file download and installation
#######################################################

# Compilers etc

AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([C])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

# network download, uncompress etc

AC_CHECK_PROG([GZIP], [gzip], [found], [notfound])
if test "x$GZIP" = "xnotfound"; then
    AC_MSG_ERROR([gzip not found in path])
fi

AC_CHECK_PROG([TAR], [tar], [found], [notfound])
if test "x$TAR" = "xnotfound"; then
    AC_MSG_ERROR([tar not found in path])
fi

AC_CHECK_PROGS([WGET], [wget curl], [notfound])
if test "x$WGET" = "xnotfound"; then
    AC_MSG_ERROR([wget or curl not found in path])
elif test "x$WGET" = "xcurl"; then
    AC_SUBST([WGET_ARGS], ["--location --progress-bar -o -"])
else
    AC_SUBST([WGET_ARGS], ["--progress=dot:mega -O -"])
fi

AC_ARG_VAR(BUILD, [Build ID for the programs])
AC_SUBST(BUILD, [$BUILD])

# programs for installation

AC_PROG_SED
AC_PROG_GREP
AC_PROG_INSTALL

#######################################################
# C dependencies
#######################################################

AC_CHECK_LIB(z, inflate,  [], [zlib_devel=notfound])
if test "x$zlib_devel" = "xnotfound"; then
    AC_MSG_ERROR([zlib-devel not found])
fi

AC_SEARCH_LIBS([log], [m], [],
  [AC_MSG_ERROR([log() not found

msamtools requires a working floating-point math library.
FAILED.  This error must be resolved in order to build msamtools successfully.])])

AC_SEARCH_LIBS([pthread_create], [pthread], [],
  [AC_MSG_ERROR([pthread_create() not found

msamtools requires a working pthread library.
FAILED.  This error must be resolved in order to build msamtools successfully.])])

#######################################################
# multiple architecture support?
# --enable-multi-arch=yes   installs external software and maui under software_dir/architecture/software
#######################################################

# figure out host type

# We use AC_CANONICAL_HOST to figure it out.
# It requires that we distribute config.guess, config.sub and install.sh with the distro.

AC_CANONICAL_HOST

# We only support static executables in linux
# 10.11.2019. Removed static builds in linux, so removed -static.

case $host_os in
   *linux*) AC_SUBST([MYCFLAGS], ["-W -Wall -Wwrite-strings -std=gnu99 -pedantic"]) ;;
         *) AC_SUBST([MYCFLAGS], ["-W -Wall -Wwrite-strings -std=gnu99 -pedantic"]) ;;
esac

AC_SUBST([MACHINE_TYPE], [$host_cpu])
kernel_name=`uname -s`
AC_SUBST([KERNEL_NAME], [$kernel_name])

AC_ARG_ENABLE([multi-arch], 
  [AS_HELP_STRING([--enable-multi-arch=ARG], [add support to multiple architectures in a common installation point @<:@default=no@:>@])],
  [],
  [enable_multi_arch=no])

case $enable_multi_arch in
   yes|no)  AC_SUBST([MULTIARCH], [$enable_multi_arch]);;
        *)  AC_MSG_ERROR([Unknown option --enable-multi-arch=$enable_multi_arch]);;
esac

AM_CONDITIONAL([MULTI_ARCH_SUPPORT], [test "x$enable_multi_arch" = "xyes"])
if test "x$enable_multi_arch" = "xyes"; then
	SOFTWAREDIR="$datarootdir/$host"
else
	SOFTWAREDIR="$datarootdir"
fi
AC_SUBST([SOFTWAREDIR], [$SOFTWAREDIR])

#######################################################
# user specified input through environmental variables
#######################################################

# substitutions

# local directories

AC_SUBST([HOST],    [generic])
AC_SUBST([AUX_DIR], [$aux_files_dir])

# make outputs

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([deps/Makefile])
AC_CONFIG_FILES([deps/argtable/Makefile])
AC_CONFIG_FILES([deps/samtools/Makefile])
AS_IF([test -f "data/samtools/Makefile"], [rm data/samtools/Makefile])


AC_OUTPUT

echo "
-------------------------------------------------------------------------------
msamtools was successfully configured with the following options:

   Host type          : $host
   Multi-architecture : $enable_multi_arch
   Install directory  : $prefix

You can build msamtools executables by typing:

   make

You can install msamtools now by typing:

   make install
"
