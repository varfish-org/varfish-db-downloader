## Rules related to Weghorn (2019) gene annotation.


rule genes_shet:  # -- postprocess file for HGNC gene IDs
    input:
        tsv="bundled-data/weghorn_2019/Weghorn_2019_Supplementary_Table_1.txt.gz",
        xlink=f"output/full/mehari/genes-xlink-{DV.today}/genes-xlink.tsv",
    output:
        tsv="work/genes/shet/2019/shet_weghorn_2019.tsv",
    shell:
        """
        QSV_SKIP_FORMAT_CHECK=1 \
        qsv join -d '\t' \
            'Gene' <(zcat {input.tsv}) \
            gene_symbol {input.xlink} \
        | qsv select 'hgnc_id,low_det' \
        | qsv rename 'hgnc_id,s_het' \
        | qsv sort \
        | tr ',' '\t' \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv}.md5
        """
