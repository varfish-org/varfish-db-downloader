## Rules related to AlphaMissense per-gene scores


rule genes_alphamissense_download:  # -- download AlphaMissense per-gene scores
    output:
        tsv="work/download/genes/alphamissense/1/AlphaMissense_gene_hg38.tsv.gz",
        tsv_md5="work/download/genes/alphamissense/1/AlphaMissense_gene_hg38.tsv.gz.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.tsv} \
            https://storage.googleapis.com/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz

        md5sum {output.tsv} > {output.tsv_md5}
        """
