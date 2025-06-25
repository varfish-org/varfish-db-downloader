rule grch37_gnomad_constraints_v2_1_1_download:
    output:
        "GRCh37/gnomAD_constraints/v2.1.1/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

        cd $(dirname {output})
        md5sum $(basename {output}) > $(basename {output}).md5
        """


rule result_grch37_gnomad_constraints_v2_1_1_tsv:
    input:
        txt="GRCh37/gnomAD_constraints/v2.1.1/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        header="header/gnomadconstraints.txt",
    output:
        tsv="GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.tsv",
        release_info="GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.txt} \
            | tail -n +2 \
            | sort -u -S {config[sort_memory]} \
            | awk -F $'\t' '
                BEGIN {{ OFS = FS }}
                {{
                    for (i=1; i<=NF; ++i) {{
                        if ($i == "NA") {{
                            $i = ""
                        }}
                    }}
                    print
                }}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadConstraints\tv2.1.1\tGRCh37\t" > {output.release_info}
        """
