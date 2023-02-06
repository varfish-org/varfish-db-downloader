rule grch38_hgmd_public_download:
    output:
        "GRCh38/hgmd_public/ensembl_r104/download/variation_feature.txt.gz",
    log:
        "GRCh38/hgmd_public/ensembl_r104/download/variation_feature.txt.gz_log",
    shell:
        r"""
        wget --no-check-certificate \
            -o {log} \
            -O {output} \
            http://ftp.ensembl.org/pub/release-104/mysql/homo_sapiens_variation_104_38/variation_feature.txt.gz
        """


rule result_grch38_hgmd_public_to_tsv:
    input:
        header="header/hgmd_public.txt",
        txt="GRCh38/hgmd_public/ensembl_r104/download/variation_feature.txt.gz",
    output:
        tsv="GRCh38/hgmd_public/ensembl_r104/HgmdPublicLocus.tsv",
        release_info="GRCh38/hgmd_public/ensembl_r104/HgmdPublicLocus.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            pigz -p 2 -d -c {input.txt} \
            | zgrep -F -w HGMD_MUTATION \
            | awk \
                -F $'\t' \
                'BEGIN {{
                    OFS=FS;
                    map[131550] = "chr1";
                    map[131545] = "chr2";
                    map[131551] = "chr3";
                    map[131552] = "chr4";
                    map[131542] = "chr5";
                    map[131555] = "chr6";
                    map[131559] = "chr7";
                    map[131560] = "chr8";
                    map[131540] = "chr9";
                    map[131544] = "chr10";
                    map[131556] = "chr11";
                    map[131546] = "chr12";
                    map[131541] = "chr13";
                    map[131547] = "chr14";
                    map[131558] = "chr15";
                    map[131549] = "chr16";
                    map[131554] = "chr17";
                    map[131548] = "chr18";
                    map[131537] = "chr19";
                    map[131538] = "chr20";
                    map[131543] = "chr21";
                    map[131557] = "chr22";
                    map[132907] = "chrM";
                    map[131539] = "chrX";
                    map[131553] = "chrY";
                }}
                ($7 == "HGMD_MUTATION") {{
                    chrom = map[$2];
                    print "GRCh38", chrom, $3 - 1, $4"\t", $8;
                }}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHgmdPublicLocus\tensembl_r104\tGRCh38\t" > {output.release_info}
        """
