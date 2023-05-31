BEGIN {
    OFS = FS;  # output separator = input separator (TAB)
    print "#chrom", "start", "end", "ensembl_gene_id";
}
($3 == "gene") {
    chrom = $1;
    start = $4 - 1;
    end = $5;

    match($9, /ENSG[0-9]+/);
    ensembl_gene_id = substr($9, RSTART, RLENGTH);

    print chrom, start, end, ensembl_gene_id;
}
