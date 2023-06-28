## Rules to create annonars RocksDB for (Wetzel & Darbro, 2022) data.


rule output_worker_patho_mms:
    output:
        tracks_bed=f"output/full/tracks/track-strucvars-patho-mms-{{genome_release}}-{{v_patho_mms}}+{DV.tracks}/patho-mms.bed",
        tracks_bed_md5=f"output/full/tracks/track-strucvars-patho-mms-{{genome_release}}-{{v_patho_mms}}+{DV.tracks}/patho-mms.bed.md5",
        tracks_spec=f"output/full/tracks/track-strucvars-patho-mms-{{genome_release}}-{{v_patho_mms}}+{DV.tracks}/patho-mms.spec.yaml",
        worker_bed=f"output/full/worker/patho-mms-{{genome_release}}-{{v_patho_mms}}+{PV.worker}/patho-mms.bed",
        worker_bed_md5=f"output/full/worker/patho-mms-{{genome_release}}-{{v_patho_mms}}+{PV.worker}/patho-mms.bed.md5",
        worker_spec=f"output/full/worker/patho-mms-{{genome_release}}-{{v_patho_mms}}+{PV.worker}/patho-mms.spec.yaml",
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_patho_mms=RE_VERSION,
    shell:
        r"""
        base=data/patho-mms/{wildcards.v_patho_mms}/patho-mms-{wildcards.genome_release}

        cp $base.bed       {output.tracks_bed}
        cp $base.bed.md5   {output.tracks_bed_md5}
        cp $base.spec.yaml {output.tracks_spec}

        cp $base.bed       {output.worker_bed}
        cp $base.bed.md5   {output.worker_bed_md5}
        cp $base.spec.yaml {output.worker_spec}
        """
