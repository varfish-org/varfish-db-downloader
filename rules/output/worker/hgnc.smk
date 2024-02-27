## Convert HGNC info to binary for worker.


rule output_hgnc_xlink_binary:
    input:
        tsv="output/full/mehari/genes-xlink-{date}/genes-xlink.tsv",
    output:
        bin=f"output/full/worker/genes-xlink-{{date}}+{PV.worker}/genes-xlink.bin",
        spec_yaml=f"output/full/worker/genes-xlink-{{date}}+{PV.worker}/genes-xlink.spec.yaml",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        varfish-server-worker db to-bin \
            --input-type xlink \
            --path-input {input.tsv} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/hgnc_xlink.spec.yaml \
            --value today={TODAY} \
            \
            --value date={wildcards.date} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
