## Rules related to Mehari transcripts.


rule genes_mehari_data_tx_download:  # -- Download the HGNC data
    output:
        zstd="work/download/mehari-data-tx/{genome_release}-{version}/mehari-data-txs-{genome_release}-{refseq_or_ensembl}-{version}.bin.zst",
    shell:
        r"""
        wget \
            -O {output.zstd} \
            https://github.com/bihealth/mehari-data-tx/releases/download/v{wildcards.version}/mehari-data-txs-{wildcards.genome_release}-{wildcards.refseq_or_ensembl}-{wildcards.version}.bin.zst
        """


rule genes_mehari_data_tx_copy:  # -- Copy data to output
    input:
        zstd="work/download/mehari-data-tx/{genome_release}-{version}/mehari-data-txs-{genome_release}-{refseq_or_ensembl}-{version}.bin.zst",
    output:
        zstd="output/full/mehari/genes-txs-{genome_release}-{version}/mehari-data-txs-{genome_release}-{refseq_or_ensembl}-{version}.bin.zst",
        zstd_md5="output/full/mehari/genes-txs-{genome_release}-{version}/mehari-data-txs-{genome_release}-{refseq_or_ensembl}-{version}.bin.zst.md5",
    shell:
        r"""
        cp {input.zstd} {output.zstd}

        md5sum {output.zstd} > {output.zstd_md5}
        """
