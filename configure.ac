AC_PREREQ([2.59])
AC_INIT([rpn], [0.3.3], [claudio.bonati82@gmail.com])
AM_INIT_AUTOMAKE(subdir-objects)
AC_CONFIG_SRCDIR([./include/macro.h])
AC_CONFIG_HEADERS([config.h])


# overwrite default values
if test x${ac_cv_env_CFLAGS_set} != xset; then
  CFLAGS=''
fi 
if test x${ac_cv_env_CPPFLAGS_set} != xset; then
  CPPFLAGS=''
fi 



### GENERAL CHECKS


# Checks for programs.
AC_PROG_AWK
AC_PROG_CC


# Checks for header files.
AC_CHECK_HEADERS([limits.h], [], [AC_MSG_ERROR([Required header limits.h not found.])])

AC_CHECK_HEADERS([math.h], [], [AC_MSG_ERROR([Required header math.h not found.])])

AC_CHECK_HEADERS([openssl/md5.h], [], [AC_MSG_ERROR([Required openssl/md5.h not found.])])

AC_CHECK_HEADERS([stdio.h], [], [AC_MSG_ERROR([Required header stdio.h not found.])])

AC_CHECK_HEADERS([stdlib.h], [], [AC_MSG_ERROR([Required header stdlib.h not found.])])

AC_CHECK_HEADERS([time.h], [], [AC_MSG_ERROR([Required header time.h not found.])])


# Checks for libraries.
AC_SEARCH_LIBS([MD5_Init], [crypto], [], [AC_MSG_ERROR([Required crypto library not found.])])

AC_SEARCH_LIBS([sqrt], [m], [], [AC_MSG_ERROR([Required c mathematical library not found.])]) 

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T



### ENABLING & CO
AC_ARG_ENABLE([use-openmp], 
              [  --enable-use-openmp     to enable the use of openmp], 
              [use_openmp=yes], 
              [use_openmp=no])


### FURTHER CHECKS
if test x${use_openmp} = xyes; then
  #check for openmp options
  AC_OPENMP

  #check for headers
  AC_CHECK_HEADERS([omp.h], [], [AC_MSG_ERROR([Required header omp.h not found.])])
fi



### VARIABLES

AC_ARG_VAR(N_f, the number of flavours (default 2))

AC_ARG_VAR(Num_threads, number of threads to be used in OpenMP (default 1))

AC_ARG_VAR(ST_dim, spacetime dimensionality (default 2))


### DEFINE & CO
if test x${use_openmp} = xyes; then
  AC_DEFINE(OPENMP_MODE,, [Enable the use of openmp])
fi


### TEST

if test x${N_f} = x; then
  AC_MSG_WARN([ 
   *****
   No chosen value for N_f: using the default N_f=2
   ***** 
   ])
  AC_DEFINE_UNQUOTED(NFLAVOUR, 2, [Number of flavours])
else
  AC_DEFINE_UNQUOTED(NFLAVOUR, ${N_f}, [Number of flavours])
fi

if test x${use_openmp} = xyes; then
  if test x${Num_threads} = x; then
    AC_MSG_WARN([
   *****
   No chosen value for Num_threads: using the default value 1
   *****
   ])
    AC_DEFINE_UNQUOTED(NTHREADS, 1, [Number of threads to be used in OpenMP])
  else
    AC_DEFINE_UNQUOTED(NTHREADS, ${Num_threads}, [Number of threads to be used in OpenMP])
  fi
else
  AC_DEFINE_UNQUOTED(NTHREADS, 1, [Number of threads to be used in OpenMP])
fi

if test x${ST_dim} = x; then
  AC_MSG_WARN([ 
   *****
   No chosen value for ST_dim: using the default ST_dim=2
   ***** 
   ])
  AC_DEFINE_UNQUOTED(STDIM, 2, [Spacetime dimensionality])
else
  AC_DEFINE_UNQUOTED(STDIM, ${ST_dim}, [Spacetime dimensionality])
fi


### CONDITIONALS


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
