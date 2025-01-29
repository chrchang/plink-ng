What this updates: https://zzz.bwh.harvard.edu/plink/

Main methods paper: https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533

PLINK 1.9 and 2.0 user documentation: https://www.cog-genomics.org/plink/1.9/ , https://www.cog-genomics.org/plink/2.0/

Technical support forum: https://groups.google.com/g/plink2-users

The 1.9/ implementation can typically be used as a drop-in replacement for PLINK 1.07 that scales to much larger datasets.  It's technically still a beta version because there are a few rarely-used but possibly-worthwhile PLINK 1.07 commands that are still absent, but active feature development for it ended in 2016.

The 2.0/ implementation is designed to handle VCF files and dosage data, and is under active development.  Most basic features other than non-concatenating merge are now in place.  See its README.md for more details.
