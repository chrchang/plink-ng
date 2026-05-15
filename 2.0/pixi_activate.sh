#!/bin/bash
# Activation script for pixi environment
# Ensures the Makefile finds conda-provided headers and libraries

export CFLAGS="${CFLAGS:-} -I${CONDA_PREFIX}/include"
export CXXFLAGS="${CXXFLAGS:-} -I${CONDA_PREFIX}/include"
export LDFLAGS="${LDFLAGS:-} -L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib"

# LIBRARY_PATH is searched by the linker for -l flags, even when the
# Makefile hardcodes its own -L paths (e.g. the Darwin override).
export LIBRARY_PATH="${CONDA_PREFIX}/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
