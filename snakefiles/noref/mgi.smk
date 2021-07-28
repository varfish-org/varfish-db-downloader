rule GRChxx_mgi_download:
    output:
        "noref/mgi/{download_date}/download/HOM_MouseHumanSequence.rpt",
    shell:
        r"""
        wget http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt \
            -O {output}
        """


rule result_GRChxx_mgi_tsv:
    input:
        txt="noref/mgi/{download_date}/download/HOM_MouseHumanSequence.rpt",
        header="header/mgihommousehumansequence.txt",
    output:
        tsv="noref/mgi/{download_date}/MgiHomMouseHumanSequence.tsv",
        release_info="noref/mgi/{download_date}/MgiHomMouseHumanSequence.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.txt} \
            | awk -F $'\t' '
                BEGIN {{ OFS=FS }}
                {{
                    for (i=11; i<14; ++i) {{
                        gsub(/,/, "\"\"\",\"\"\"", $i);
                        if ($i=="") {{
                            $i="{{}}"
                        }}
                        else {{
                            $i="{{\"\"\""$i"\"\"\"}}"
                        }}
                    }}
                    print
                }}'
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nMgiHomMouseHumanSequence\t{wildcards.download_date}\t\t" > {output.release_info}
        """
