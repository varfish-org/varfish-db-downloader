## Convert HGNC info to binary for worker.


rule output_hgnc_xlink_binary:
    input:
        tsv="output/full/worker/genes-xlink-{date}/genes-xlink.tsv",
    output:
        bin="output/full/worker/genes-xlink-{date}/genes-xlink.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type xlink \
            --path-input {input.tsv} \
            --path-output-bin {output.bin}
        """
