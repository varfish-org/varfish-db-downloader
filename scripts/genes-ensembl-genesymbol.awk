# Create the ENSEMBL gene ID to gene symbol table from the ENSEMBL GTF file.
#
# Genes without a symbol are skipped.
BEGIN {
    OFS = FS;  # output separator = input separator (TAB)
}
($3 == "gene") {
    split($9, a, " ");

    ensg = "";
    symbol = "";

    for (i = 1; i <= length(a); i++) {
        if (a[i] == "gene_id") {
            ensg = a[i + 1];
        }
        if (a[i] == "gene_name") {
            symbol = a[i + 1];
        }
    }

    gsub(/[";]/, "", ensg);
    gsub(/[";]/, "", symbol);

    if (ensg != "" && symbol != "") {
        print ensg, symbol;
    }
}
