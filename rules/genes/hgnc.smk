## Rules related to the HGNC data.


rule genes_hgnc_download:  # -- Download the HGNC data
    output:
        json="work/download/hgnc/hgnc_complete_set.json",
        json_md5="work/download/hgnc/hgnc_complete_set.json.md5",
        date_txt="work/download/hgnc/date.txt",
    shell:
        r"""
        date > {output.date_txt}

        wget --no-check-certificate \
            -O {output.json} \
            https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json

        md5sum {output.json} > {output.json_md5}
        """


rule genes_hgnc_xlink:  # -- Build HGNC xlink table.
    input:
        json="work/download/hgnc/hgnc_complete_set.json",
    output:
        tsv="work/genes/hgnc/hgnc_xlink.tsv",
        tsv_md5="work/genes/hgnc/hgnc_xlink.tsv.md5",
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
        json="work/download/hgnc/hgnc_complete_set.json",
    output:
        jsonl="work/genes/hgnc/hgnc_info.jsonl",
        jsonl_md5="work/genes/hgnc/hgnc_info.jsonl.md5",
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
