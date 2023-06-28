## Generate binary version of masked regions for worker.


rule output_masked_repeat:
    input:
        f"output/full/tracks/track-features-ucsc-rmsk-{{genome_release}}-{{version}}+{DV.tracks}/rmsk.bed.gz",
    output:
        bin=f"output/full/worker/masked-repeat-{{genome_release}}-{{version}}+{PV.worker}/masked-repeat.bin",
        spec_yaml=f"output/full/worker/masked-repeat-{{genome_release}}-{{version}}+{PV.worker}/masked-repeat.spec.yaml",
    wildcard_constraints:
        genome_release=RE_GENOME,
        version=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type masked-region \
            --path-input {input} \
            --path-output-bin {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/masked_repeat.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.version} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """


rule output_masked_segdup:
    input:
        f"output/full/tracks/track-features-ucsc-genomicsuperdups-{{genome_release}}-{{version}}+{DV.tracks}/genomicSuperDups.bed.gz",
    output:
        bin=f"output/full/worker/masked-segdup-{{genome_release}}-{{version}}+{PV.worker}/masked-segdup.bin",
        spec_yaml=f"output/full/worker/masked-segdup-{{genome_release}}-{{version}}+{PV.worker}/masked-segdup.spec.yaml",
    wildcard_constraints:
        genome_release=RE_GENOME,
        version=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type masked-region \
            --path-input {input} \
            --path-output-bin {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/masked_segdup.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.version} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
