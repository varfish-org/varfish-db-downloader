## Rules related to Collins (2022) gene annotation.


rule genes_rcnv_download:  # -- download pHaplo/pTriplo scores
    output:
        tsv="work/download/genes/rcnv/2022/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz",
        tsv_md5="work/download/genes/rcnv/2022/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.tsv} \
            https://zenodo.org/record/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz

        md5sum {output.tsv} > {output.tsv_md5}
        """


rule genes_rcnv_postproces:  # -- postprocess file for HGNC gene IDs
    input:
        tsv="work/download/genes/rcnv/2022/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz",
        xlink=f"output/full/mehari/genes-xlink-{DV.today}/genes-xlink.tsv",
    output:
        tsv="work/genes/rcnv/2022/rcnv_collins_2022.tsv",
        tsv_md5="work/genes/rcnv/2022/rcnv_collins_2022.tsv.md5",
    shell:
        """
        qsv join -d '\t' \
            '#gene' <(zcat {input.tsv}) \
            gene_symbol {input.xlink} \
        | qsv select 'hgnc_id,pHaplo,pTriplo' \
        | qsv rename 'hgnc_id,p_haplo,p_triplo' \
        | qsv sort \
        | tr ',' '\t' \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv}.md5
        """
