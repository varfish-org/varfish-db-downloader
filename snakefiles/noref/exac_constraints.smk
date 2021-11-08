rule noref_exac_constraints_r0_3_1_download:
    output:
        "noref/ExAC_constraints/r0.3.1/download/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output} \
            https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz
        cd $(dirname {output})
        md5sum $(basename {output}) > $(basename {output}).md5
        """


rule result_noref_exac_constraints_r0_3_1_tsv:
    input:
        txt="noref/ExAC_constraints/r0.3.1/download/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz",
        header="header/exacconstraints.txt",
    output:
        tsv="noref/ExAC_constraints/r0.3.1/ExacConstraints.tsv",
        release_info="noref/ExAC_constraints/r0.3.1/ExacConstraints.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.txt} \
            | tail -n +2 \
            | sort -u -S {config[sort_memory]} \
            | cut -f 1-20 \
            | awk -F $'\t' '
                BEGIN {{ OFS = FS }}
                {{
                    sub(/\.[0-9]+$/, "", $1)
                    print $0, ".", "."
                }}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nExacConstraints\tr0.3.1\tnoref\t" > {output.release_info}
        """
