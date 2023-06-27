## Generate binary version of masked regions for worker.


rule output_masked_repeat:
    input:
        "output/full/worker/annos/features/ucsc-rmsk-{genome_release}-{version}/rmsk.bed.gz",
    output:
        "output/full/worker/masked-repeat-{genome_release}-{version}/masked-repeat.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type masked-region \
            --path-input {input} \
            --path-output-bin {output}
        """


rule output_masked_segdup:
    input:
        "output/full/worker/annos/features/ucsc-genomicsuperdups-{genome_release}-{version}/genomicSuperDups.bed.gz",
    output:
        "output/full/worker/masked-segdup-{genome_release}-{version}/masked-segdup.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type masked-region \
            --path-input {input} \
            --path-output-bin {output}
        """
