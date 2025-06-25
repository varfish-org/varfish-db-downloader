rule result_grch38_gnomad_constraints_v2_1_1_tsv:
    input:
        txt="work/download/genes/gnomad/2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        md5="work/download/genes/gnomad/2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz.md5",
        header="rules/pre-mehari/header/gnomadconstraints.txt",
    output:
        tsv="output/pre-mehari/GRCh38/gnomAD_constraints/v2.1.1/GnomadConstraints.tsv",
        release_info="output/pre-mehari/GRCh38/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.txt} \
            | tail -n +2 \
            | sort -u \
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

        echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadConstraints\tv2.1.1\tGRCh38\t" > {output.release_info}
        """
