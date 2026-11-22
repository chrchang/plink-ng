PROGRAM: PLINK

DESCRIPTION: Whole-genome association analysis toolset

AUTHOR: Shaun Purcell

CONTACT: plink@chgr.mgh.harvard.edu

YEAR: 2006, 2007

LICENSE: Released under GNU General Public License, v2 (see
COPYING.txt)

DOCUMENTATION: http://pngu.mgh.harvard.edu/purcell/plink/

INSTALLATION: If you have download a zip or gzipped archive with an
executable binary, no installation is necessary (except perhaps you
might want to place the executable in your path, see documentation for
details). Otherwise, see notes on compilation below.

COMPILATION: You will need a standard C/C++ compiler such as GNU gcc
(version 3). This is likely available on all Linux/Unix platforms. For
MS-DOS, DJGPP or MinGW are appropriate choices. To help compiling, see
documentation (basically, just be sure to select the correct Makefile
and type make -f Makefile.*)

USAGE: Type "plink" or "./plink" from the command line followed by the
options of choice (see documentation)

EXAMPLE DATA: Two example files test.ped and test.map are included in
the distribution; for example, once PLINK is installed try running:

     plink --file test

     plink --file test --freq

     plink --file test --assoc

     plink --file test --make-bed

     plink --bfile test --assoc

     etc...

SMP, Aug 2006

