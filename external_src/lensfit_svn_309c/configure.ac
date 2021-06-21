AC_INIT([lensfit], [7.7], [L.Miller@physics.ox.ac.uk])
AM_INIT_AUTOMAKE(foreign)

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
date=`date +%Y-%m-%d`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_SUBST(CFLAGS,$CFLAGS)
AC_SUBST(LDFLAGS, $LDFLAGS)

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stddef.h stdlib.h string.h unistd.h sys/mman.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_MMAP
AC_FUNC_STRFTIME
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([strdup])
AC_CHECK_PROG(HAVE_GSL, gsl-config, [yes], [no])
CFLAGS="$CFLAGS `gsl-config --cflags`"
LDFLAGS="$LDFLAGS `gsl-config --libs`"
AC_CHECK_LIB(m,cos)

# GSL
AC_CHECK_LIB(gslcblas,cblas_dgemm, [], \
                 [echo ; echo Linking against gsl failed.;
                 echo Make sure that GSL with CBLAS is installed and that GSL is in your PATH; echo ;
                 exit 1])
AC_CHECK_LIB(gsl,gsl_blas_dgemm, [], \
                 [echo ; echo Linking against gsl failed.;
                 echo Make sure that GSL with CBLAS is installed and that GSL is in your PATH; echo ;
                 exit 1])
AC_CHECK_HEADER(fitsio.h, [], \
                          [echo ; echo fitsio.h not found.;\
                          echo Make sure that cfitsio is installed.; echo ;
                          exit 1], \
                          [])
# cfitsio
AC_CHECK_LIB(cfitsio, ffopen, [], \
                      [echo ; echo Linking against cfitsio failed.;
                      echo Make sure that cfitsio is installed.; echo ;
                      exit 1], \
                      [-lm])
AC_CHECK_HEADER(fitsio.h, [], \
                          [echo ; echo fitsio.h not found.;\
                          echo Make sure that cfitsio is installed.; echo ;
                          exit 1], \
                          [])
AC_CHECK_LIB(cfitsio, ffopen, [], \
                      [echo ; echo Linking against cfitsio failed.;
                      echo Make sure that cfitsio is installed.; echo ;
                      exit 1], \
                      [-lm])

# fftw
AC_CHECK_HEADER(fftw3.h, [], \
			 [echo ; echo fftw3.h not found.;\
                          echo Make sure that fftw3 is installed.; echo ;
                          exit 1], \
                          [])
AC_CHECK_LIB(fftw3, fftw_malloc, [], \
                      [echo ; echo Linking against libfftw3 failed.;
                      echo Make sure that fftw3 is installed.; echo ;
                      exit 1], \
                      [-lm])


AC_CONFIG_FILES(Makefile \
        swarp/Makefile \
	swarp/fits/Makefile \
	swarp/wcs/Makefile \
	nr/Makefile \
	src/Makefile)
AC_OUTPUT
