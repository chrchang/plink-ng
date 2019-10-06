plink-ng/1.9
============

This is mostly a large collection of report-generating functions which do not
depend on each other.  plink_common contains convenience functions for file
parsing and bit manipulation, plink_matrix encapsulates a few common matrix
operations (which are resolved via LAPACK calls under the hood), and plink.c
handles command-line parsing and initial dispatch; it's often not necessary to
deal with any other files when adding a new command.

Additional discussion is at https://www.cog-genomics.org/plink/1.9/dev .
