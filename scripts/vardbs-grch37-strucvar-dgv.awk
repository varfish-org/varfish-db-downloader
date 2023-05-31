BEGIN {
    OFS = FS;
    print "#chromosome", "begin", "end", "sv_type", "observed_gains", "observed_losses";
}
($1 !~ /^variantaccession/) {
    if ($16 + $17 > 0) {
        print $2, $3 - 1, $4, $6, $16, $17;
    }
}
