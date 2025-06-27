rule result_noref_acmg:
    input:
        "data/acmg_sf/{v_acmg_sf}/acmg_sf.tsv",
    output:
        tsv="output/pre-mehari/noref/acmg/{v_acmg_sf}/Acmg.tsv",
        release_info="output/pre-mehari/noref/acmg/{v_acmg_sf}/Acmg.release_info",
    shell:
        r"""
        (
            echo -e "ensembl_gene_id\tsymbol\tentrez_id"
            tail -n +2 {input} | awk -F $"\t" 'BEGIN{{OFS=FS}}{{print $2,$4,$3}}'
        ) > {output.tsv}
        echo -e "table\tversion\tgenomebuild\tnull_value\nAcmg\t{wildcards.v_acmg_sf}\t\t" > {output.release_info}
        """
