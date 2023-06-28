## Rules to prepare bgdb files for worker.


rule output_worker_bgdb_g1k:
    input:
        bed=f"output/full/tracks/track-strucvars-g1k-grch37-{{version}}+{DV.tracks}/g1k.bed.gz",
    output:
        bin="output/full/worker/bgdb-g1k-grch37-{version}/bgdb-g1k.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type strucvar-g1k \
            --path-input {input.bed} \
            --path-output-bin {output.bin}
        """


rule output_worker_bgdb_exac:
    input:
        bed=f"output/full/tracks/track-strucvars-exac-grch37-{{version}}+{DV.tracks}/exac.bed.gz",
    output:
        bin="output/full/worker/bgdb-exac-grch37-{version}/bgdb-exac.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type strucvar-exac-cnv \
            --path-input {input.bed} \
            --path-output-bin {output.bin}
        """


rule output_worker_bgdb_gnomad:
    input:
        bed=f"output/full/tracks/track-strucvars-gnomad-grch37-{{version}}+{DV.tracks}/gnomad.bed.gz",
    output:
        bin="output/full/worker/bgdb-gnomad-grch37-{version}/bgdb-gnomad.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type strucvar-gnomad-sv \
            --path-input {input.bed} \
            --path-output-bin {output.bin}
        """


rule output_worker_bgdb_dbvar:
    input:
        bed=f"output/full/tracks/track-strucvars-dbvar-{{genome_release}}-{{version}}+{DV.tracks}/dbvar.bed.gz",
    output:
        bin="output/full/worker/bgdb-dbvar-{genome_release}-{version}/bgdb-dbvar.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type strucvar-db-var \
            --path-input {input.bed} \
            --path-output-bin {output.bin}
        """


rule output_worker_bgdb_dgv:
    input:
        bed=f"output/full/tracks/track-strucvars-dgv-{{genome_release}}-{{version}}+{DV.tracks}/dgv.bed.gz",
    output:
        bin="output/full/worker/bgdb-dgv-{genome_release}-{version}/bgdb-dgv.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type strucvar-dgv \
            --path-input {input.bed} \
            --path-output-bin {output.bin}
        """


rule output_worker_bgdb_dgv_gs:
    input:
        bed=f"output/full/tracks/track-strucvars-dgv-gs-{{genome_release}}-{{version}}+{DV.tracks}/dgv-gs.bed.gz",
    output:
        bin="output/full/worker/bgdb-dgv-gs-{genome_release}-{version}/bgdb-dgv-gs.bin",
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type strucvar-dgv-gs \
            --path-input {input.bed} \
            --path-output-bin {output.bin}
        """
