BEGIN {
    OFS = FS;
    printf("#chromosome\tbegin\tend\tsv_type\n");
}
($0 ~ /^track name=dupControls/) {
    sv_type = "DUP";
}
($0 ~ /^track name=delControls/) {
    sv_type = "DEL";
}
($1 !~ /^track/) {
    chrom = substr($1, 4);

    start = $2;
    end = $3;

    printf( \
        "%s\t%s\t%s\t%s\n", \
        chrom, \
        start, \
        end, \
        sv_type \
    );
}
