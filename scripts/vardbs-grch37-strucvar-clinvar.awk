BEGIN {
    OFS = FS;
    print \
        "#chromosome", "begin", "end", "bin", "reference", \
        "alternative", "clinvar_version", "set_type", "variation_type", \
        "symbols", "hgnc_ids", "vcv", "summary_clinvar_review_status_label", \
        "summary_clinvar_pathogenicity_label", "summary_clinvar_pathogenicity", \
        "summary_clinvar_gold_stars", "summary_paranoid_review_status_label", \
        "summary_paranoid_pathogenicity_label", "summary_paranoid_pathogenicity", \
        "summary_paranoid_gold_stars", "details";
}
(NR > 1) {
    $3 = $3 - 1;
    for (i = 2; i < NF; i++) {
        printf("%s%s", $i, OFS);
    }
    print $NF
}
