rule GRChxx_mgi_download:
    output:
        "noref/mgi/latest/download/HOM_MouseHumanSequence.rpt",
    shell:
        r"""
        wget http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt \
            -O {output}
        """


rule GRChxx_mgi_tsv:
    input:
        txt="noref/mgi/latest/download/HOM_MouseHumanSequence.rpt",
        header="header/mgihommousehumansequence.txt",
    output:
        tsv="noref/mgi/latest/MgiHomMouseHumanSequence.tsv",
        release_info="noref/mgi/latest/MgiHomMouseHumanSequence.release_info",
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

        echo -e "table\tversion\tgenomebuild\tnull_value\nMgiHomMouseHumanSequence\tlatest\t{wildcards.genome_build}\t" > {output.release_info}
        """
