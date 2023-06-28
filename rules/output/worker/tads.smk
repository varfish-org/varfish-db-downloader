## Rules related to TADs from worker.


rule output_worker_tads:  # -- generate worker TAD files
    input:
        bed_hesc=f"output/full/tracks/track-tads-{{genome_release}}-dixon2015+{DV.tracks}/hesc.bed",
    output:
        bed_hesc="output/full/worker/tads-{genome_release}-dixon2015/hesc.bed",
        spec_yaml="output/full/worker/tads-{genome_release}-dixon2015/spec.yaml",
    shell:
        r"""
        cp {input.bed_hesc} {output.bed_hesc}

        varfish-db-downloader tpl \
            --template rules/output/worker/tads.spec.yaml \
            --template rules/output/worker/hgnc_xlink.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
