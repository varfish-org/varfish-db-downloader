## Rules related to exome probesets.


rule output_tracks_exome_probeset:
    output:
        bed="output/full/tracks/track-enrichment-probesets-targets/{name}-{genome_release}-{v_tracks}.bed.gz",
        bed_tbi="output/full/tracks/track-enrichment-probesets-targets/{name}-{genome_release}-{v_tracks}.bed.gz.tbi",
    wildcard_constraints:
        name=r"[a-zA-Z0-9_-]+",
        genome_release=RE_GENOME,
        v_tracks=RE_VERSION,
    shell:
        r"""
        cp bundled-data/track-enrichment-probesetss/{wildcards.name}-targets-{wildcards.genome_release}.bed.gz {output.bed}
        cp bundled-data/track-enrichment-probesetss/{wildcards.name}-targets-{wildcards.genome_release}.bed.gz.tbi {output.bed_tbi}
        """
