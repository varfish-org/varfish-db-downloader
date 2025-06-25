## Rules related to the HGNC data.


rule genes_hgnc_xlink:  # -- Build HGNC xlink table.
    input:
        json="work/download/annos/hgnc/{quarterly_release_date}/hgnc_complete_set.json",
    output:
        tsv="output/full/mehari/genes-xlink-{quarterly_release_date}/genes-xlink.tsv",
        tsv_md5="output/full/mehari/genes-xlink-{quarterly_release_date}/genes-xlink.tsv.md5",
    shell:
        r"""
        jq \
            --raw-output \
            --from-file scripts/genes-xlink-hgnc.jq \
            {input.json} \
        > {output.tsv}

        md5sum {output.tsv} > {output.tsv_md5}
        """


rule genes_hgnc_gene_info:  # -- Build HGNC gene_info JSONL file.
    input:
        json="work/download/hgnc/{quarterly_release_date}/hgnc_complete_set.json",
    output:
        jsonl="work/genes/hgnc/{quarterly_release_date}/hgnc_info.jsonl",
        jsonl_md5="work/genes/hgnc/{quarterly_release_date}/hgnc_info.jsonl.md5",
    shell:
        r"""
        jq \
            --compact-output \
            --raw-output \
            --from-file scripts/genes-hgnc-info.jq \
            {input.json} \
        > {output.jsonl}

        md5sum {output.jsonl} > {output.jsonl_md5}
        """
