#!/bin/sh -e

##########################################################################
#   Script description:
#       Wrapper to provide build options to generic Makefile
#
#   Arguments:
#       1.  Build template: dev
#       2.  Additional arguments are passed to make
#            (meant to be used for make targets install, clean, ...)
#
#   History:
#   Date        Name        Modification
#   2018-05-06  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

case $(uname) in
    Darwin)
        BLASFLAGS64='-framework Accelerate'
        LDFLAGS='-ldl -lpthread -lz'
        MAKE=make
        ;;

    *)
        MAKE=make
        ;;
esac

# Defaults
BLASFLAGS64="$BLASFLAGS64 -llapack -lcblas -lblas"

if [ $# -ge 1 ]; then
    template=$1
    case $template in
        avx2+atlas)
            CFLAGS="$CFLAGS -mavx2 -mbmi -mbmi2 -mlzcnt"
            BLASFLAGS64="$BLASFLAGS64 -llapack -lf77blas -latlas"
            shift
            ;;

        atlas)
            BLASFLAGS64="$BLASFLAGS64 -llapack -lf77blas -latlas"
            shift
            ;;

        netlib)
            BLASFLAGS64="$BLASFLAGS64 -llapack -lcblas -lblas"
            shift
            ;;
        *)
            # Pass arg 1 to make
            ;;
    esac
fi

CFLAGS="$CFLAGS -DNDEBUG -DZSTD_MULTITHREAD"
CXXFLAGS="$CXXFLAGS -DNDEBUG -DZSTD_MULTITHREAD"
set -x
export CFLAGS CXXFLAGS BLASFLAGS64 LDFLAGS ZLIB BLASFLAGS64
$MAKE -f Makefile "$@"
