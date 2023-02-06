BEGIN {
    OFS = FS;  # output separator = input separator (TAB)
    print "#chrom", "start", "end", "entrez_id";
}
(NR == FNR && $1 !~ /^#/) {  # first file only
    acc2chrom[$2] = $1
}
(NR != FNR && $3 == "gene") {  # second file
    chrom = acc2chrom[$1];
    start = $4 - 1;
    end = $5;

    match($9, /GeneID:[0-9]+/);
    entrez_id = substr($9, RSTART + 7, RLENGTH - 7);

    print chrom, start, end, entrez_id;
}
