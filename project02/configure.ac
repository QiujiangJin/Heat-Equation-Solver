#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Project01], [1.0], [qiujiang@ices.utexas.edu])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AX_PATH_MASA([0.30],[yes])
AX_PATH_GRVY([0.32],[yes])
AX_PATH_HDF5([1.8.0],[yes])
AX_PATH_PETSC([no])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Checks for code coverage
AX_CODE_COVERAGE

AC_OUTPUT(Makefile src/Makefile test/Makefile)
