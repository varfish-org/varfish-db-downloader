BEGIN {
    OFS = FS;
    print "#chromosome", "begin", "end", "sv_sub_type", "n_homalt", "n_het";
}
($1 ~ /^#CHROM/) {
    n_samples = NF - 9;
}
($1 !~ /^#/) {
    n_homalt = 0;
    n_het = 0;
    n_ref = 0;
    for (i = 10; i <= NF; i++) {
        if ($i == "0|0") {
            n_ref += 1;
        } else if ($i == "1|1" || $i == "2|2" || $i == "3|3" || $i == "4|4" || $i == "5|5" || $i == "5|5" || $i == "6|6" || $i == "7|7" || $i == "8|8" || $i == "9|9") {
            n_homalt += 1;
        } else {
            n_het += 1;
        }
    }

    begin = $2 - 1;

    match($8, /;END=[0-9]+;/);
    end = substr($8, RSTART + 5, RLENGTH - 6);
    if (end == "") {
        end = begin + 1;
    }

    match($8, /SVTYPE=[a-zA-Z0-9:_]+;/);
    sv_type = substr($8, RSTART + 7, RLENGTH - 8);

    sv_sub_type = substr($5, 2, length($5) - 2);
    pat = sv_type "$";
    if (sv_sub_type ~ pat) {
        sv_type = sv_sub_type;
    }

    print $1, begin, end, sv_type, n_homalt, n_het;
}
