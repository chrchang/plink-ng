# Modeled on https://github.com/traversc/qs/blob/master/configure.ac .

AC_INIT([pgenlibr],[0.3.5],[chrchang@alumni.caltech.edu])
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

AC_ARG_WITH([libdeflate-include],
            AS_HELP_STRING([--with-libdeflate-include=INCLUDE_PATH],[the location of libdeflate header files]),
            [libdeflate_include_path=$withval])

AC_ARG_WITH([libdeflate-lib],
            AS_HELP_STRING([--with-libdeflate-lib=LIB_PATH],[the location of libdeflate library files]),
            [libdeflate_lib_path=$withval])

AC_ARG_WITH([simde-include],
            AS_HELP_STRING([--with-simde-include=INCLUDE_PATH],[the location of simde header files]),
            [simde_include_path=$withval])

########################################################
#### Version value function

getVersion2()
{
VERSION_STRING=$1
MAJOR=`echo $VERSION_STRING | cut -d. -f1`
MINOR=`echo $VERSION_STRING | cut -d. -f2`
echo $(($MAJOR*1000+$MINOR))
}

getVersion3()
{
VERSION_STRING=$1
MAJOR=`echo $VERSION_STRING | cut -d. -f1`
MINOR=`echo $VERSION_STRING | cut -d. -f2`
RELEASE=`echo $VERSION_STRING | cut -d. -f3`
echo $(($MAJOR*100000+$MINOR*100+$RELEASE))
}

########################################################
#### Compile zstd/libdeflate/simde checks

ZSTD_INCLUDE_PATH=""
ZSTD_LIBS=""
ZSTD_SHLIB=""
ZSTD_CLEAN=""
LIBDEFLATE_INCLUDE_PATH=""
LIBDEFLATE_LIBS=""
LIBDEFLATE_SHLIB=""
LIBDEFLATE_CLEAN=""
SIMDE_INCLUDE_PATH=""

if test "xx$zstd_include_path" != "xx"; then
  AC_MSG_NOTICE([Using user-defined zstd install paths])
    ZSTD_LIBS="-L${zstd_lib_path}"
    ZSTD_INCLUDE_PATH="-I${zstd_include_path}"
    COMPILE_ZSTD="false"
elif test "xx$PKGCONF" != "xx"; then
  if "${PKGCONF}" --exists libzstd; then
    VERSION_STRING=`${PKGCONF} --modversion libzstd`
    VER=`getVersion3 ${VERSION_STRING}`
    # The bundled 1.5.5 includes a compression bugfix, but the Python/R
    # interfaces only use zstd for decompression, so we don't need it.
    # The release notes for 1.5.4 mention 32-bit bugfixes, so it is more
    # reasonable to refuse to go below it.
    if test "${VER}" -ge 100504; then
      AC_MSG_NOTICE([zstd ${VERSION_STRING} library detected -- skipping zstd compilation])
      ZSTD_LIBS=`"${PKGCONF}" --libs libzstd`
      ZSTD_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I libzstd`
      COMPILE_ZSTD="false"
    else
      AC_MSG_NOTICE([zstd ${VERSION_STRING} library detected but is too old (<1.5.4) -- compiling from source])
      COMPILE_ZSTD="true"
    fi
  else
    AC_MSG_NOTICE([zstd library not detected -- compiling from source])
    COMPILE_ZSTD="true"
  fi
else
  AC_MSG_NOTICE([pkg-config not detected -- compiling zstd from source])
  COMPILE_ZSTD="true"
fi

if test xx$COMPILE_ZSTD = "xxtrue"; then
  ZSTD_LIBS="${LIBS} -lPGZSTD"
  ZSTD_INCLUDE_PATH="-Izstd/lib -Izstd/lib/common"
  ZSTD_SHLIB="libPGZSTD.a"
  ZSTD_CLEAN="\$(LIBZSTD) libPGZSTD.a"
fi

if test "xx$libdeflate_include_path" != "xx"; then
  AC_MSG_NOTICE([Using user-defined libdeflate install paths])
    LIBDEFLATE_LIBS="-L${libdeflate_lib_path}"
    LIBDEFLATE_INCLUDE_PATH="-I${libdeflate_include_path}"
    COMPILE_LIBDEFLATE="false"
elif test "xx$PKGCONF" != "xx"; then
  if "${PKGCONF}" --exists libdeflate; then
    VERSION_STRING=`${PKGCONF} --modversion libdeflate`
    VER=`getVersion2 ${VERSION_STRING}`
    if test "${VER}" -ge 1014; then
      AC_MSG_NOTICE([libdeflate ${VERSION_STRING} library detected -- skipping libdeflate compilation])
      LIBDEFLATE_LIBS=`"${PKGCONF}" --libs libdeflate`
      LIBDEFLATE_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I libdeflate`
      COMPILE_LIBDEFLATE="false"
    else
      AC_MSG_NOTICE([libdeflate ${VERSION_STRING} library detected but is lower than bundled version (1.14) -- compiling from source])
      COMPILE_LIBDEFLATE="true"
    fi
  else
    AC_MSG_NOTICE([libdeflate library not detected -- compiling from source])
    COMPILE_LIBDEFLATE="true"
  fi
else
  AC_MSG_NOTICE([pkg-config not detected -- compiling libdeflate from source])
  COMPILE_LIBDEFLATE="true"
fi

if test xx$COMPILE_LIBDEFLATE = "xxtrue"; then
  LIBDEFLATE_LIBS="${LIBS} -lPGDEFLATE"
  LIBDEFLATE_INCLUDE_PATH="-Ilibdeflate"
  LIBDEFLATE_SHLIB="libPGDEFLATE.a"
  LIBDEFLATE_CLEAN="\$(LIBDEFLATE) libPGDEFLATE.a"
fi

if test "xx$simde_include_path" != "xx"; then
  AC_MSG_NOTICE([Using user-defined simde install paths])
    SIMDE_INCLUDE_PATH="-I${simde_include_path}"
    USE_BUNDLED_SIMDE="false"
elif test "xx$PKGCONF" != "xx"; then
  if "${PKGCONF}" --exists simde; then
    VERSION_STRING=`${PKGCONF} --modversion simde`
    VER=`getVersion3 ${VERSION_STRING}`
    if test "${VER}" -ge 702; then
      AC_MSG_NOTICE([simde ${VERSION_STRING} library detected -- ignoring bundled simde header files])
      SIMDE_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I simde`
      USE_BUNDLED_SIMDE="false"
    else
      AC_MSG_NOTICE([simde ${VERSION_STRING} library detected but is lower than bundled version (0.7.2) -- using bundled header files])
      USE_BUNDLED_SIMDE="true"
    fi
  else
    AC_MSG_NOTICE([simde library not detected -- using bundled header files])
    USE_BUNDLED_SIMDE="true"
  fi
else
  AC_MSG_NOTICE([pkg-config not detected -- using bundled simde header files])
  USE_BUNDLED_SIMDE="true"
fi

if test xx$USE_BUNDLED_SIMDE = "xxtrue"; then
  SIMDE_INCLUDE_PATH="-Isimde"
fi

AC_SUBST([ZSTD_INCLUDE_PATH], $ZSTD_INCLUDE_PATH)
AC_SUBST([ZSTD_LIBS], $ZSTD_LIBS)
AC_SUBST([ZSTD_SHLIB], $ZSTD_SHLIB)
AC_SUBST([ZSTD_CLEAN], $ZSTD_CLEAN)
AC_SUBST([LIBDEFLATE_INCLUDE_PATH], $LIBDEFLATE_INCLUDE_PATH)
AC_SUBST([LIBDEFLATE_LIBS], $LIBDEFLATE_LIBS)
AC_SUBST([LIBDEFLATE_SHLIB], $LIBDEFLATE_SHLIB)
AC_SUBST([LIBDEFLATE_CLEAN], $LIBDEFLATE_CLEAN)
AC_SUBST([SIMDE_INCLUDE_PATH], $SIMDE_INCLUDE_PATH)


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
