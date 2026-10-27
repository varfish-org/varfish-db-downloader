# Create the ENSEMBL gene/transcript/Entrez/symbol xlink table.
#
# Input files: the ENSEMBL "entrez.tsv" dump, then the ENSEMBL GTF file.
# Output: one line per transcript and Entrez ID (empty Entrez ID if none).
BEGIN {
    OFS = FS;  # output separator = input separator (TAB)
}
# ENSEMBL entrez.tsv: transcript -> Entrez gene IDs.
(NR == FNR) {
    if ($5 == "EntrezGene" && $2 != "-") {
        if ($2 in entrez) {
            entrez[$2] = entrez[$2] SUBSEP $4;
        } else {
            entrez[$2] = $4;
        }
    }
    next;
}
($3 == "transcript") {
    split($9, a, " ");

    ensg = "";
    enst = "";
    symbol = "";

    for (i = 1; i <= length(a); i++) {
        if (a[i] == "gene_id") {
            ensg = a[i + 1];
        }
        if (a[i] == "transcript_id") {
            enst = a[i + 1];
        }
        if (a[i] == "gene_name") {
            symbol = a[i + 1];
        }
    }

    gsub(/[";]/, "", ensg);
    gsub(/[";]/, "", enst);
    gsub(/[";]/, "", symbol);

    if (ensg != "" && enst != "") {
        if (enst in entrez) {
            n = split(entrez[enst], ids, SUBSEP);
            for (j = 1; j <= n; j++) {
                print ensg, enst, ids[j], symbol;
            }
        } else {
            print ensg, enst, "", symbol;
        }
    }
}
