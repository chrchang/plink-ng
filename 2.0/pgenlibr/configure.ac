# Adapted from https://github.com/traversc/qs/blob/master/configure.ac and
# https://github.com/ropensci/git2r/blob/master/configure.ac .

AC_PREREQ([2.69])
AC_INIT([pgenlibr],[0.5.3],[chrchang@alumni.caltech.edu])
AC_PATH_PROG([PKGCONF],[pkg-config],[],[$PATH:/usr/local/bin:ext/bin:ext:/sw/bin:/opt/bin:/opt/local/bin])

# Find the compiler and compiler flags to use for fallback version checks
AC_LANG(C++)
: ${R_HOME=`R RHOME`}
PATH_TO_CPP_COMPILER=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX([$PATH_TO_CPP_COMPILER])
AC_MSG_NOTICE([C++ compiler: $PATH_TO_CPP_COMPILER])
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

# Configure args
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

# pkg-config version value functions
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

# Compile zstd/libdeflate/simde checks
ZSTD_COMPILE_FLAG=""
ZSTD_INCLUDE_PATH=""
ZSTD_LIBS=""
ZSTD_SHLIB=""
ZSTD_CLEAN=""
COMPILE_ZSTD=""
LIBDEFLATE_COMPILE_FLAG=""
LIBDEFLATE_INCLUDE_PATH=""
LIBDEFLATE_LIBS=""
LIBDEFLATE_SHLIB=""
LIBDEFLATE_CLEAN=""
COMPILE_LIBDEFLATE=""
SIMDE_COMPILE_FLAG=""
SIMDE_INCLUDE_PATH=""
USE_BUNDLED_SIMDE=""

if test "xx$zstd_include_path" != "xx"; then
  AC_MSG_NOTICE([Using user-defined zstd install paths])
  ZSTD_INCLUDE_PATH="-I${zstd_include_path}"
  ZSTD_LIBS="-L${zstd_lib_path}"
  COMPILE_ZSTD="false"
else
  if test "xx$PKGCONF" != "xx"; then
    if "${PKGCONF}" --exists libzstd; then
      VERSION_STRING=`${PKGCONF} --modversion libzstd`
      VER=`getVersion3 ${VERSION_STRING}`
      # The bundled 1.5.5-1.5.7 includes compression bugfixes relative to
      # 1.5.2, but the Python/R interfaces currently only use zstd for
      # decompression, so we don't need it.
      # (Tried to bundle 1.5.7 for CRAN, but that triggered a
      # diagnostic-suppression complaint that did not exist for 1.5.5.  So CRAN
      # submissions are probably pinned to 1.5.5 indefinitely.)
      # 1.5.2 is over three years old, and most R build machines have at least
      # that version.
      if test "${VER}" -ge 100502; then
        AC_MSG_NOTICE([zstd ${VERSION_STRING} library detected -- skipping zstd compilation])
        ZSTD_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I libzstd`
        ZSTD_LIBS=`"${PKGCONF}" --libs libzstd`
        COMPILE_ZSTD="false"
      else
        AC_MSG_NOTICE([zstd ${VERSION_STRING} library detected but is too old (<1.5.2) -- compiling from source])
        COMPILE_ZSTD="true"
      fi
    fi
  fi
  if test "xx$COMPILE_ZSTD" = "xx"; then
    AC_MSG_NOTICE([zstd not detected by pkg-config; falling back on system header file check])
    zstd_ver_ok=no
    ${PATH_TO_CPP_COMPILER} ${CPPFLAGS} ${CXXFLAGS} -E tools/zstd_version.cpp >/dev/null 2>&1 && zstd_ver_ok=yes
    AC_MSG_RESULT([${zstd_ver_ok}])
    if test "xx${zstd_ver_ok}" = "xxyes"; then
      ZSTD_LIBS="-lzstd"
      COMPILE_ZSTD="false"
    else
      COMPILE_ZSTD="true"
    fi
  fi
fi

if test xx$COMPILE_ZSTD = "xxtrue"; then
  ZSTD_COMPILE_FLAG="-DZSTD_DISABLE_ASM"
  ZSTD_INCLUDE_PATH="-Izstd/lib -Izstd/lib/common"
  ZSTD_LIBS="${LIBS} -lPGZSTD"
  ZSTD_SHLIB="libPGZSTD.a"
  ZSTD_CLEAN="\$(LIBZSTD) libPGZSTD.a"
else
  ZSTD_COMPILE_FLAG="-DIGNORE_BUNDLED_ZSTD"
fi

if test "xx$libdeflate_include_path" != "xx"; then
  AC_MSG_NOTICE([Using user-defined libdeflate install paths])
    LIBDEFLATE_INCLUDE_PATH="-I${libdeflate_include_path}"
    LIBDEFLATE_LIBS="-L${libdeflate_lib_path}"
    COMPILE_LIBDEFLATE="false"
else
  if test "xx$PKGCONF" != "xx"; then
    if "${PKGCONF}" --exists libdeflate; then
      VERSION_STRING=`${PKGCONF} --modversion libdeflate`
      VER=`getVersion2 ${VERSION_STRING}`
      # 1.14 is over two years old, and includes significant performance
      # improvements over its predecessor.  Subsequent updates have been less
      # significant.
      if test "${VER}" -ge 1014; then
        AC_MSG_NOTICE([libdeflate ${VERSION_STRING} library detected -- skipping libdeflate compilation])
        LIBDEFLATE_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I libdeflate`
        LIBDEFLATE_LIBS=`"${PKGCONF}" --libs libdeflate`
        COMPILE_LIBDEFLATE="false"
      else
        AC_MSG_NOTICE([libdeflate ${VERSION_STRING} library detected but is too old (<1.9) -- compiling from source])
        COMPILE_LIBDEFLATE="true"
      fi
    fi
  fi
  if test "xx$COMPILE_LIBDEFLATE" = "xx"; then
    AC_MSG_NOTICE([libdeflate not detected by pkg-config; falling back on system header file check])
    libdeflate_ver_ok=no
    ${PATH_TO_CPP_COMPILER} ${CPPFLAGS} ${CXXFLAGS} -E tools/libdeflate_version.cpp >/dev/null 2>&1 && libdeflate_ver_ok=yes
    AC_MSG_RESULT([${libdeflate_ver_ok}])
    if test "xx${libdeflate_ver_ok}" = "xxyes"; then
      LIBDEFLATE_LIBS="-ldeflate"
      COMPILE_LIBDEFLATE="false"
    else
      COMPILE_LIBDEFLATE="true"
    fi
  fi
fi

if test xx$COMPILE_LIBDEFLATE = "xxtrue"; then
  LIBDEFLATE_INCLUDE_PATH="-Ilibdeflate"
  LIBDEFLATE_LIBS="${LIBS} -lPGDEFLATE"
  LIBDEFLATE_SHLIB="libPGDEFLATE.a"
  LIBDEFLATE_CLEAN="\$(LIBDEFLATE) libPGDEFLATE.a"
else
  LIBDEFLATE_COMPILE_FLAG="-DIGNORE_BUNDLED_LIBDEFLATE"
fi

if test "xx$simde_include_path" != "xx"; then
  AC_MSG_NOTICE([Using user-defined simde install paths])
  SIMDE_INCLUDE_PATH="-I${simde_include_path}"
  USE_BUNDLED_SIMDE="false"
else
  if test "xx$PKGCONF" != "xx"; then
    if "${PKGCONF}" --exists simde; then
      VERSION_STRING=`${PKGCONF} --modversion simde`
      VER=`getVersion3 ${VERSION_STRING}`
      if test "${VER}" -ge 702; then
        AC_MSG_NOTICE([simde ${VERSION_STRING} library detected -- ignoring bundled simde header files])
        SIMDE_INCLUDE_PATH=`"${PKGCONF}" --cflags-only-I simde`
        USE_BUNDLED_SIMDE="false"
      else
        AC_MSG_NOTICE([simde ${VERSION_STRING} library detected but is too old (<0.7.2) -- using bundled header files])
        USE_BUNDLED_SIMDE="true"
      fi
    fi
  fi
  if test "xx$USE_BUNDLED_SIMDE" = "xx"; then
    AC_MSG_NOTICE([simde not detected by pkg-config; falling back on system header file check])
    simde_ver_ok=no
    ${PATH_TO_CPP_COMPILER} ${CPPFLAGS} ${CXXFLAGS} -E tools/simde_version.cpp >/dev/null 2>&1 && simde_ver_ok=yes
    AC_MSG_RESULT([${simde_ver_ok}])
    if test "xx${simde_ver_ok}" = "xxyes"; then
      USE_BUNDLED_SIMDE="false"
    else
      USE_BUNDLED_SIMDE="true"
    fi
  fi
fi

if test xx$USE_BUNDLED_SIMDE != "xxtrue"; then
  SIMDE_COMPILE_FLAG="-DIGNORE_BUNDLED_SIMDE"
fi


AC_SUBST([ZSTD_COMPILE_FLAG], $ZSTD_COMPILE_FLAG)
AC_SUBST([ZSTD_INCLUDE_PATH], $ZSTD_INCLUDE_PATH)
AC_SUBST([ZSTD_LIBS], $ZSTD_LIBS)
AC_SUBST([ZSTD_SHLIB], $ZSTD_SHLIB)
AC_SUBST([ZSTD_CLEAN], $ZSTD_CLEAN)
AC_SUBST([LIBDEFLATE_COMPILE_FLAG], $LIBDEFLATE_COMPILE_FLAG)
AC_SUBST([LIBDEFLATE_INCLUDE_PATH], $LIBDEFLATE_INCLUDE_PATH)
AC_SUBST([LIBDEFLATE_LIBS], $LIBDEFLATE_LIBS)
AC_SUBST([LIBDEFLATE_SHLIB], $LIBDEFLATE_SHLIB)
AC_SUBST([LIBDEFLATE_CLEAN], $LIBDEFLATE_CLEAN)
AC_SUBST([SIMDE_COMPILE_FLAG], $SIMDE_COMPILE_FLAG)
AC_SUBST([SIMDE_INCLUDE_PATH], $SIMDE_INCLUDE_PATH)


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
