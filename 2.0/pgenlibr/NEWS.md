# pgenlibr 0.4.0 (2025-01-15)
CHROM and POS columns in .pvar are now loaded by default, and they can be
checked with GetVariantChrom() and GetVariantPos().  NewPvar() should now work
properly on .bim files.

# pgenlibr 0.3.7 (2024-06-03)
Undo premature switch to system libdeflate in Makevars.win, which broke R 4.3
Windows build.

# pgenlibr 0.3.6 (2024-05-27)
Patches Windows-build link order, makes .pgen reader forward-compatible with
header/footer extensions in May 2024 specification update.

# pgenlibr 0.3.5 (2023-06-12)
Patches irregularities (e.g. interaction with preinstalled libraries, unaligned
loads) in initial submission.

# pgenlibr 0.3.3 (2023-05-24)
Initial CRAN version.
