## Rules related to Orphadata download


rule genes_orphadata_download:  # -- download orphadatas
    output:
        jsonl="work/genes/orphadata/{date}/orphadata.jsonl",
        jsonl_md5="work/genes/orphadata/{date}/orphadata.jsonl.md5",
        tsv="work/genes/orphadata/{date}/orpha_diseases.tsv",  # xxx
    shell:
        """
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        python ./scripts/genes-orpha-diseases.py \
        > {output.jsonl}

        md5sum {output.jsonl} > {output.jsonl}.md5

        echo -e "hgnc_id\torpha_id\tassoc_status\tomim_ids\tdisease_name\tdefinition" \
        > {output.tsv}
        """
