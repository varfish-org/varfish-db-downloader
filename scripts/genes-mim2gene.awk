BEGIN {
    OFS = FS;
    print "omim_id\tentrez_id"
}
($3 == "phenotype" && $2 != "-") {
    print $1, $2
}
