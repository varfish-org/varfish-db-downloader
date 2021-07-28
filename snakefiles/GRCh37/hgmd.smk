# Download the variation track from GRCh37 ENSEMBL and extract HGMD Public entries, prepare
# TSV file for import.


rule grch37_hgmd_public_download:
    output:
        "GRCh37/hgmd_public/ensembl_r104/download/variation_feature.txt.gz",
    log:
        "GRCh37/hgmd_public/ensembl_r104/download/variation_feature.txt.gz_log",
    shell:
        r"""
        wget \
            -o {log} \
            -O {output} \
            http://ftp.ensembl.org/pub/grch37/release-104/mysql/homo_sapiens_variation_104_37/variation_feature.txt.gz
        """


rule result_grch37_hgmd_public_to_tsv:
    input:
        header="header/hgmd_public.txt",
        txt="GRCh37/hgmd_public/ensembl_r104/download/variation_feature.txt.gz",
    output:
        tsv="GRCh37/hgmd_public/ensembl_r104/HgmdPublicLocus.tsv",
        release_info="GRCh37/hgmd_public/ensembl_r104/HgmdPublicLocus.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            gunzip -c {input.txt} \
            | zgrep -F -w HGMD_MUTATION \
            | awk \
                -F $'\t' \
                'BEGIN {{
                    OFS=FS;
                    map[100965601] = "MT";
                    map[27504] = "11";
                    map[27505] = "21";
                    map[27506] = "7";
                    map[27507] = "Y";
                    map[27508] = "2";
                    map[27509] = "17";
                    map[27510] = "22";
                    map[27511] = "1";
                    map[27512] = "18";
                    map[27513] = "13";
                    map[27514] = "16";
                    map[27515] = "6";
                    map[27516] = "X";
                    map[27517] = "3";
                    map[27518] = "9";
                    map[27519] = "12";
                    map[27520] = "14";
                    map[27521] = "15";
                    map[27522] = "20";
                    map[27523] = "8";
                    map[27524] = "4";
                    map[27525] = "10";
                    map[27526] = "19";
                    map[27527] = "5";
                }}
                ($7 == "HGMD_MUTATION") {{
                    chrom = map[$2];
                    print "GRCh37", chrom, $3 - 1, $4"\t", $8;
                }}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHgmdPublicLocus\tensembl_r104\tGRCh37\t" > {output.release_info}
        """
