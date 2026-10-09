#!/bin/bash
# Activation script for pixi environment
# Ensures the Makefile finds conda-provided headers and libraries

# CPATH is searched by gcc/clang for headers regardless of -I flags set by
# the Makefile (which hard-assigns CFLAGS/CXXFLAGS, ignoring the environment).
export CPATH="${CONDA_PREFIX}/include${CPATH:+:$CPATH}"

# LIBRARY_PATH is searched by the linker for -l flags, even when the
# Makefile hardcodes its own -L paths (e.g. the Darwin override).
export LIBRARY_PATH="${CONDA_PREFIX}/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"

# LDFLAGS is appended to the link command so the binary can find shared
# libraries at runtime (bakes in the rpath).
export LDFLAGS="${LDFLAGS:-} -Wl,-rpath,${CONDA_PREFIX}/lib"

# build_dynamic/Makefile builds its own LINKFLAGS and never reads LDFLAGS, so
# on GNU ld the rpath has to come from LD_RUN_PATH instead.
export LD_RUN_PATH="${CONDA_PREFIX}/lib${LD_RUN_PATH:+:$LD_RUN_PATH}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
