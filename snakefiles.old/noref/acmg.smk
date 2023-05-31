rule result_noref_acmg:
    input:
        "tools/data/Acmg.tsv",
    output:
        tsv="noref/acmg/v3.1/Acmg.tsv",
        release_info="noref/acmg/v3.1/Acmg.release_info",
    shell:
        r"""
        cp {input} {output.tsv}
        echo -e "table\tversion\tgenomebuild\tnull_value\nAcmg\tv3.1\t\t" > {output.release_info}
        """
