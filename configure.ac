AC_INIT([measim],[1.0],[cjj@uw.edu])
AC_CONFIG_SRCDIR([measim.cc])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_PROG_CXX
AM_PROG_AR
LT_INIT
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
AX_LIB_HDF5
if test $with_hdf5 = no; then
	AC_MSG_ERROR([HDF5 library is required])
fi
AS_VERSION_COMPARE(HDF5_VERSION,[1.8],[
	AC_MSG_ERROR([HDF5 library version >= 1.8 is necessary])
])
AM_CONDITIONAL([HAVE_HDF5], [test "x$no_hdf5" != "xyes"])
PROG_CFLAGS="$HDF5_CFLAGS $HDF5_CPPFLAGS"
PROG_LIBS="$HDF5_LDFLAGS"

AX_CHECK_GLUT
AM_CONDITIONAL([HAVE_GLUT], [test "x$no_glut" != "xyes"])
AM_COND_IF([HAVE_GLUT],
	PROG_CFLAGS="$PROG_CFLAGS $GLUT_CFLAGS"
	PROG_LIBS="$PROG_LIBS $GLUT_LIBS"
)

dnl variables needed
AC_SUBST(PROG_CFLAGS)
AC_SUBST(PROG_LIBS)
AC_CONFIG_FILES([
	Makefile
	cst/Makefile
])
AC_OUTPUT
