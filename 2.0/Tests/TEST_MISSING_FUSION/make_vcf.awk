# Emits a small VCF spanning autosomes, chrX, chrY and chrMT, with missing
# genotypes throughout.  Positions start well clear of PAR1 so the chrX
# variants stay on chrX under --split-par.
BEGIN {
  OFS = "\t";
  ns = 120;
  base = 5000000;
  srand(11);
  print "##fileformat=VCFv4.2";
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">";
  split("1 2 X Y MT", chrs, " ");
  for (i = 1; i <= 5; ++i) {
    print "##contig=<ID=" chrs[i] ">";
  }
  line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (s = 0; s != ns; ++s) {
    line = line "\ts" sprintf("%03d", s);
  }
  print line;
  split("200 150 120 60 40", cts, " ");
  vid = 0;
  for (i = 1; i <= 5; ++i) {
    c = chrs[i];
    for (j = 0; j != cts[i]; ++j) {
      ++vid;
      triallelic = (multi && (j % 4 == 0));
      row = c OFS (base + j * 137) OFS "v" vid OFS "A" OFS (triallelic? "G,T" : "G") OFS "." OFS "." OFS "." OFS "GT";
      for (s = 0; s != ns; ++s) {
        # sample sex alternates, with every 17th sample left unknown
        male = (s % 2 == 0);
        hap = ((c == "X" && male) || c == "Y" || c == "MT");
        if (c == "Y" && !male) {
          row = row OFS "./.";
          continue;
        }
        if (rand() < 0.06) {
          row = row OFS (hap? "." : "./.");
          continue;
        }
        if (hap) {
          if (triallelic) {
            r = rand();
            row = row OFS ((r < 0.34)? "0" : ((r < 0.67)? "1" : "2"));
          } else {
            row = row OFS ((rand() < 0.5)? "0" : "1");
          }
        } else if (triallelic) {
          r = rand();
          row = row OFS ((r < 0.2)? "0/0" : ((r < 0.4)? "0/1" : ((r < 0.6)? "1/1" : ((r < 0.8)? "0/2" : "1/2"))));
        } else {
          r = rand();
          row = row OFS ((r < 0.34)? "0/0" : ((r < 0.67)? "0/1" : "1/1"));
        }
      }
      print row;
    }
  }
}
