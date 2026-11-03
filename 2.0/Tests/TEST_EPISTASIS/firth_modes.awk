# Checks that the default 'firth-fallback' report is made up of the 'firth'
# and 'no-firth' statistics: each pair has one or the other (bit for bit), the
# Firth value is used where it differs from the logistic one, and both
# kinds turn up.
FILENAME == ARGV[1] { if (FNR > 1) { firth[$2 "|" $4] = $5 }; next }
FILENAME == ARGV[2] { if (FNR > 1) { nofirth[$2 "|" $4] = $5 }; next }
FNR > 1 {
    k = $2 "|" $4
    if (!(k in firth)) { print "fallback pair missing from the firth report: " k; exit 1 }
    if ($5 == firth[k]) {
        ++is_firth
        if ((!(k in nofirth)) || (nofirth[k] != $5)) { ++rescued }
    } else if ((k in nofirth) && ($5 == nofirth[k])) {
        ++is_logistic
    } else {
        print "fallback STAT " $5 " matches neither mode on " k; exit 1
    }
    ++n
}
END {
    firthlen = 0
    for (item in firth) {
        firthlen++
    }
    if (n != firthlen) { print "fallback has " n " pairs, firth " firthlen; exit 1 }
    if (!is_firth || !is_logistic || !rescued) {
        print "expected a mix: " is_firth+0 " firth, " is_logistic+0 " logistic, " rescued+0 " rescued"; exit 1
    }
    print n " pairs: " is_logistic " logistic, " is_firth " Firth (" rescued " of them differ from logistic)"
}
