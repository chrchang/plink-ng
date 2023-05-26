# Modeled on https://github.com/traversc/qs/blob/master/configure.ac .

AC_INIT([pgenlibr],[0.3.4],[chrchang@alumni.caltech.edu])
AC_PATH_PROG([PKGCONF],[pkg-config],[],[$PATH:/usr/local/bin:ext/bin:ext:/sw/bin:/opt/bin:/opt/local/bin])

########################################################
AC_LANG(C++)
: ${R_HOME=`R RHOME`}
PATH_TO_CPP_COMPILER=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX([$PATH_TO_CPP_COMPILER])
echo "C++ compiler: $PATH_TO_CPP_COMPILER"

########################################################
### Configure args

AC_ARG_WITH([zstd-include],
            AS_HELP_STRING([--with-zstd-include=INCLUDE_PATH],[the location of zstd header files]),
            [zstd_include_path=$withval])

AC_ARG_WITH([zstd-lib],
            AS_HELP_STRING([--with-zstd-lib=LIB_PATH],[the location of zstd library files]),
            [zstd_lib_path=$withval])

########################################################
#### Version value function

getVersion()
{
VERSION_STRING=$1
MAJOR=`echo $VERSION_STRING | cut -d. -f1`
MINOR=`echo $VERSION_STRING | cut -d. -f2`
RELEASE=`echo $VERSION_STRING | cut -d. -f3`
echo $(($MAJOR*100000+$MINOR*100+$RELEASE))
}

########################################################
#### Compile ZSTD checks

ZSTD_INCLUDE_PATH=""
ZSTD_LIBS=""
ZSTD_SHLIB=""
ZSTD_CLEAN=""

if test "xx$zstd_include_path" != "xx"; then
  echo "Using user-defined zstd install paths"
    ZSTD_LIBS="-L${zstd_lib_path}"
    ZSTD_INCLUDE_PATH="-I${zstd_include_path}"
    COMPILE_ZSTD="false"
elif test "xx$PKGCONF" != "xx"; then
  if "${PKGCONF}" --exists libzstd; then
    VERSION_STRING=`${PKGCONF} --modversion libzstd`
    VER=`getVersion ${VERSION_STRING}`
    if test "${VER}" -ge 100505; then
      echo "zstd ${VERSION_STRING} library detected -- skipping zstd compilation"
      ZSTD_LIBS=`"${PKGCONF}" --libs libzstd`
      ZSTD_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I libzstd`
      COMPILE_ZSTD="false"
    else
      echo "zstd ${VERSION_STRING} library detected but is lower than bundled version (1.5.5) -- compiling from source"
      COMPILE_ZSTD="true"
    fi
  else
    echo "zstd library not detected -- compiling from source"
    COMPILE_ZSTD="true"
  fi
else
  echo "pkg-config not detected -- compiling zstd from source"
  COMPILE_ZSTD="true"
fi

if test xx$COMPILE_ZSTD = "xxtrue"; then
  ZSTD_LIBS="${LIBS} -lPGZSTD"
  ZSTD_INCLUDE_PATH="-Izstd/lib -Izstd/lib/common"
  ZSTD_SHLIB="libPGZSTD.a"
  ZSTD_CLEAN="\$(LIBZSTD) libPGZSTD.a"
fi

AC_SUBST([ZSTD_INCLUDE_PATH], $ZSTD_INCLUDE_PATH)
AC_SUBST([ZSTD_LIBS], $ZSTD_LIBS)
AC_SUBST([ZSTD_SHLIB], $ZSTD_SHLIB)
AC_SUBST([ZSTD_CLEAN], $ZSTD_CLEAN)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
