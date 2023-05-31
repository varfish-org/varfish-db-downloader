rule result_grch37_release_jannovar_db:
    input:
        "GRCh37/jannovar/ensembl_87_hg19.ser",
        "GRCh37/jannovar/refseq_curated_105_hg19.ser",
    output:
        tar="releases/{release_name}/jannovar-db-{release_name}-grch37.tar.gz".format(**config),
        tar_sha256="releases/{release_name}/jannovar-db-{release_name}-grch37.tar.gz.sha256".format(
            **config
        ),
    shell:
        r"""
        tar \
            -czvf {output.tar} \
            --transform "s|^.*/|jannovar-db-{config[release_name]}-grch37/|" \
            {input}
        pushd $(dirname {output.tar})
        sha256sum $(basename {output.tar}) >$(basename {output.tar_sha256})
        """


rule result_grch38_release_jannovar_db:
    input:
        "GRCh38/jannovar/ensembl_91_hg38.ser",
        "GRCh38/jannovar/refseq_curated_109_hg38.ser",
    output:
        tar="releases/{release_name}/jannovar-db-{release_name}-grch38.tar.gz".format(**config),
        tar_sha256="releases/{release_name}/jannovar-db-{release_name}-grch38.tar.gz.sha256".format(
            **config
        ),
    shell:
        r"""
        tar \
            --owner=0 \
            --group=0 \
            -czvf {output.tar} \
            --transform "s|^.*/|jannovar-db-{config[release_name]}-grch38/|" \
            {input}
        pushd $(dirname {output.tar})
        sha256sum $(basename {output.tar}) >$(basename {output.tar_sha256})
        """
