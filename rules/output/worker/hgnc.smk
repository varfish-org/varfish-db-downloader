## Convert HGNC info to binary for worker.


rule output_hgnc_xlink_binary:
    input:
        tsv="output/full/mehari/genes-xlink-{hgnc_quarterly}/genes-xlink.tsv",
    output:
        bin="output/full/worker/genes-xlink-{hgnc_quarterly}+{worker}/genes-xlink.bin",
        spec_yaml="output/full/worker/genes-xlink-{hgnc_quarterly}+{worker}/genes-xlink.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type xlink \
            --path-input {input.tsv} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/hgnc_xlink.spec.yaml \
            --value today={TODAY} \
            \
            --value date={wildcards.hgnc_quarterly} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
