## Convert HGNC info to binary for worker.


rule output_hgnc_xlink_binary:
    input:
        tsv="output/full/mehari/genes-xlink-{date}/genes-xlink.tsv",
    output:
        bin=f"output/full/worker/genes-xlink-{{date}}+{PV.worker}/genes-xlink.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type xlink \
            --path-input {input.tsv} \
            --path-output-bin {output.bin}
        """
