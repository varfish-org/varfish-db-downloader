rule result_grch38_jannovar_download_ensembl:
    output:
        "GRCh38/jannovar/ensembl_91_hg38.ser",
    shell:
        r"wget --no-check-certificate -O {output} https://zenodo.org/record/4916051/files/$(basename {output})?download=1"


rule result_grch38_jannovar_download_refseq:
    output:
        "GRCh38/jannovar/refseq_curated_109_hg38.ser",
    shell:
        r"wget --no-check-certificate -O {output} https://zenodo.org/record/4916051/files/$(basename {output})?download=1"
