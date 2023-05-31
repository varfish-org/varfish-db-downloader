BEGIN {
    OFS = FS;
    printf("#chromosome\tbegin\tend\tnum_carriers\tsv_type\n");
}
($1 !~ /^#/) {
    print $1, $2 - 1, $3, $4, $5;
}
