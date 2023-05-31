BEGIN {
    OFS = FS;
    print "#chromosome", "begin_outer", "end_outer", "sv_sub_type", "num_carriers";
}
($1 !~ /^variantaccession/) {
    chrom = substr($1, 4)

    match($9, /outer_start=[0-9]+;/);
    outer_start = substr($9, RSTART + 12, RLENGTH - 13);

    match($9, /outer_end=[0-9]+;/);
    outer_end = substr($9, RSTART + 10, RLENGTH - 11);

    match($9, /variant_sub_type=[a-zA-Z_]+;/);
    variant_sub_type = substr($9, RSTART + 17, RLENGTH - 18);

    match($9, /num_samples=[0-9]+;/);
    num_samples = substr($9, RSTART + 12, RLENGTH - 13);

    key = chrom "-" outer_start "-" outer_end "-" variant_sub_type "-" num_samples;
    if (prev != key) {
        print chrom, outer_start - 1, outer_end, variant_sub_type, num_samples;
    }
    prev = key;
}
