rule result_grch37_release_annotator_db:
    input:
        "GRCh37/varfish-annotator-db/varfish-annotator-db-{release_name}-grch37.h2.db".format(
            **config
        ),
    output:
        gz="releases/{release_name}/varfish-annotator-db-{release_name}-grch37.h2.db.gz".format(
            **config
        ),
        sha256="releases/{release_name}/varfish-annotator-db-{release_name}-grch37.h2.db.gz.sha256".format(
            **config
        ),
    shell:
        r"""
        gzip -c {input} >{output.gz}
        pushd $(dirname {output.gz})
        sha256sum $(basename {output.gz}) >$(basename {output.gz}).sha256
        """


rule result_grch38_release_annotator_db:
    input:
        "GRCh38/varfish-annotator-db/varfish-annotator-db-{release_name}-grch38.h2.db".format(
            **config
        ),
    output:
        gz="releases/{release_name}/varfish-annotator-db-{release_name}-grch38.h2.db.gz".format(
            **config
        ),
        sha256="releases/{release_name}/varfish-annotator-db-{release_name}-grch38.h2.db.gz.sha256".format(
            **config
        ),
    shell:
        r"""
        gzip -c {input} >{output.gz}
        pushd $(dirname {output.gz})
        sha256sum $(basename {output.gz}) >$(basename {output.gz}).sha256
        """
