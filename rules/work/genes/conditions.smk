## Rules for aggregating per-gene conditions data for annonars.


rule work_conditions_integrate:  # integrate conditions file
    input:
        mim2gene_medgen="work/download/genes/ncbi/{date}/mim2gene_medgen",
        genes_xlink="output/full/mehari/genes-xlink-{date}/genes-xlink.tsv",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        orpha_jsonl="work/genes/orphadata/{date}/orphadata.jsonl",
        panelapp_jsonl="work/download/genes/panelapp/{date}/panelapp.jsonl",
        mondo_obo="work/genes/mondo/{date}/mondo.obo",
        mondo_unmapped_tsv="work/genes/mondo/{date}/omim_unmapped_terms.tsv",
        ctd_tsv="work/download/genes/ctd/{date}/CTD_diseases.tsv.gz",
        do_omim_unmapped="work/download/do/{date}/omim-unmapped.csv",
        do_omim_indo="work/download/do/{date}/OMIMinDO.tsv",
        do_omim_import="work/download/do/{date}/omim_import.obo",
    output:
        jsonl="work/genes/conditions/{v_hpo}+{date}/conditions.jsonl",
    shell:
        r"""
        MIM2GENE_MEDGEN_PATH="{input.mim2gene_medgen}" \
        GENES_XLINK_PATH="{input.genes_xlink}" \
        HPOA_PATH="{input.hpoa}" \
        ORPHA_JSONL_PATH="{input.orpha_jsonl}" \
        PANELAPP_JSONL_PATH="{input.panelapp_jsonl}" \
        MONDO_OBO_PATH="{input.mondo_obo}" \
        MONDO_UNMAPPED_OMIM_PATH="{input.mondo_unmapped_tsv}" \
        CTD_PATH="{input.ctd_tsv}" \
        DO_OMIM_UNMAPPED_PATH="{input.do_omim_unmapped}" \
        DO_OMIM_INDO_PATH="{input.do_omim_indo}" \
        DO_OMIM_IMPORT_PATH="{input.do_omim_import}" \
        python scripts/genes-integrate-diseases.py \
        > {output.jsonl}
        """
