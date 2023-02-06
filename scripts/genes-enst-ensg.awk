BEGIN {
    OFS = FS;  # output separator = input separator (TAB)
}
($3 == "transcript") {
    split($9, a, " ");

    ensg = "";
    enst = "";

    for (i = 1; i <= length(a); i++) {
        if (a[i] == "gene_id") {
            ensg = a[i + 1];
        }
        if (a[i] == "transcript_id") {
            enst = a[i + 1];
        }
    }

    gsub(/[";]/, "", ensg);
    gsub(/[";]/, "", enst);

    if (ensg != "" && enst != "") {
        print enst, ensg;
    }
}
