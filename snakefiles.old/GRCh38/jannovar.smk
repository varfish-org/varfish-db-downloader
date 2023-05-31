rule result_grch37_jannovar_download_ensembl:
    output:
        "GRCh37/jannovar/ensembl_87_hg19.ser",
    shell:
        r"wget --no-check-certificate -O {output} https://zenodo.org/record/4916051/files/$(basename {output})?download=1"


rule result_grch37_jannovar_download_refseq:
    output:
        "GRCh37/jannovar/refseq_curated_105_hg19.ser",
    shell:
        r"wget --no-check-certificate -O {output} https://zenodo.org/record/4916051/files/$(basename {output})?download=1"
