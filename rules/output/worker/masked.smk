## Generate binary version of masked regions for worker.


rule output_masked_repeat:
    input:
        f"output/full/tracks/track-features-ucsc-rmsk-{{genome_release}}-{{version}}+{DV.tracks}/rmsk.bed.gz",
    output:
        f"output/full/worker/masked-repeat-{{genome_release}}-{{version}}+{PV.worker}/masked-repeat.bin",
    wildcard_constraints:
        genome_release=RE_GENOME,
        version=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type masked-region \
            --path-input {input} \
            --path-output-bin {output}
        """


rule output_masked_segdup:
    input:
        f"output/full/tracks/track-features-ucsc-genomicsuperdups-{{genome_release}}-{{version}}+{DV.tracks}/genomicSuperDups.bed.gz",
    output:
        f"output/full/worker/masked-segdup-{{genome_release}}-{{version}}+{PV.worker}/masked-segdup.bin",
    wildcard_constraints:
        genome_release=RE_GENOME,
        version=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type masked-region \
            --path-input {input} \
            --path-output-bin {output}
        """
