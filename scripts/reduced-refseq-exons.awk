# Create BED file for reduced data exons.
BEGIN {
    OFS = FS;  # output separator = input separator (TAB)
    print "#chrom", "start", "end", "entrez_id", "gene_symbol";
    print gene_symbol_re;
}
(NR == FNR && $1 !~ /^#/) {  # first file only
    acc2chrom[$2] = $1
}
(NR != FNR && $3 == "exon") {  # second file
    chrom = acc2chrom[$1];
    start = $4 - 1;
    end = $5;

    match($9, /GeneID:[0-9]+/);
    entrez_id = substr($9, RSTART + 7, RLENGTH - 7);

    match($9, /gene "[^"]+"/);
    gene_symbol = substr($9, RSTART + 6, RLENGTH - 7);

    if (length(gene_symbol_re) == 0 || match(gene_symbol, gene_symbol_re)) {
        print chrom, start, end, entrez_id, gene_symbol;
    }
}
