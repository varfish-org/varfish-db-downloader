## Rules related to Orphadata download


rule genes_orphadata_download:  # -- download orphadatas
    output:
        jsonl="work/genes/orphadata/{version}/orphadata.jsonl",
        jsonl_md5="work/genes/orphadata/{version}/orphadata.jsonl.md5",
        tsv="work/genes/orphadata/{version}/orpha_diseases.tsv",  # xxx
    shell:
        """
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        python ./scripts/genes-orpha-diseases.py \
        > {output.jsonl}

        md5sum {output.jsonl} > {output.jsonl}.md5

        echo -e "hgnc_id\torpha_id\tassoc_status\tomim_ids\tdisease_name\tdefinition" \
        > {output.tsv}
        """
