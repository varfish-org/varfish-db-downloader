rule result_noref_acmg:
    input:
        "tools/data/Acmg.tsv",
    output:
        tsv="noref/acmg/v3.0/Acmg.tsv",
        release_info="noref/acmg/v3.0/Acmg.release_info",
    shell:
        r"""
        cp {input} {output.tsv}
        echo -e "table\tversion\tgenomebuild\tnull_value\nAcmg\tv3.0\t\t" > {output.release_info}
        """
