plink-ng
========

This is mostly a large collection of report-generating functions which do not
depend on each other.  plink_common contains convenience functions for file
parsing and bit manipulation, plink_matrix encapsulates a few common matrix
operations (which are resolved via LAPACK calls under the hood), and plink.c
handles command-line parsing and initial dispatch; it's often not necessary to
deal with any other files when adding a new command.

When building directly from the code in this repository, you'll probably want
to use Makefile.std instead of Makefile; the latter is really just aimed at my
own OS X and Win64 machines.

Additional discussion is at https://www.cog-genomics.org/plink2/dev .
