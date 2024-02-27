## Rules to prepare bgdb files for worker.


rule output_worker_bgdb_g1k:
    input:
        bed=f"output/full/tracks/track-strucvars-g1k-grch37-{{version}}+{DV.tracks}/g1k.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-g1k-grch37-{{version}}+{PV.worker}/bgdb-g1k.bin",
        spec=f"output/full/worker/bgdb-g1k-grch37-{{version}}+{PV.worker}/bgdb-g1k.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-g1k \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=g1k \
            --value title="Thousand Genomes SVs" \
            --value creator="Thousand Genomes Consortium" \
            --value source="https://www.internationalgenome.org/data/" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release=grch37 \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_exac:
    input:
        bed=f"output/full/tracks/track-strucvars-exac-grch37-{{version}}+{DV.tracks}/exac.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-exac-grch37-{{version}}+{PV.worker}/bgdb-exac.bin",
        spec=f"output/full/worker/bgdb-exac-grch37-{{version}}+{PV.worker}/bgdb-exac.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-exac-cnv \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=exac \
            --value title="ExAC CNVs" \
            --value creator="ExAC Consortium" \
            --value source="https://gnomad.broadinstitute.org/downloads#exac" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release=grch37 \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_gnomad_exomes_cnv_grch38:
    input:
        bed=f"output/full/tracks/track-strucvars-gnomad-cnv-grch38-{{version}}+{DV.tracks}/gnomad-cnv.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-gnomad-exomes-cnv-grch38-{{version}}+{PV.worker}/bgdb-gnomad-exomes-cnv-grch38.bin",
        spec=f"output/full/worker/bgdb-gnomad-exomes-cnv-grch38-{{version}}+{PV.worker}/bgdb-gnomad-exomes-cnv-grch38.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-gnomad-cnv4 \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=gnomad-exomes-cnv \
            --value title="gnomAD Exomes CNV" \
            --value creator="gnomAD Consortium" \
            --value source="https://gnomad.broadinstitute.org/downloads#v4-copy-number-variants" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release=grch38 \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_gnomad_sv_grch37:
    input:
        bed=f"output/full/tracks/track-strucvars-gnomad-grch37-{{version}}+{DV.tracks}/gnomad.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-gnomad-grch37-{{version}}+{PV.worker}/bgdb-gnomad.bin",
        spec=f"output/full/worker/bgdb-gnomad-grch37-{{version}}+{PV.worker}/bgdb-gnomad.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-gnomad-sv2 \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=gnomad \
            --value title="gnomAD-SVs" \
            --value creator="gnomAD Consortium" \
            --value source="https://gnomad.broadinstitute.org/downloads#v2-structural-variants" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release=grch37 \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_gnomad_genomes_sv_grch38:
    input:
        bed=f"output/full/tracks/track-strucvars-gnomad-sv-grch38-{{version}}+{DV.tracks}/gnomad-sv.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-gnomad-genomes-sv-grch38-{{version}}+{PV.worker}/bgdb-gnomad-genomes-sv-grch38.bin",
        spec=f"output/full/worker/bgdb-gnomad-genomes-sv-grch38-{{version}}+{PV.worker}/bgdb-gnomad-genomes-sv-grch38.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-gnomad-sv4 \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=gnomad \
            --value title="gnomAD-SVs" \
            --value creator="gnomAD Consortium" \
            --value source="https://gnomad.broadinstitute.org/downloads#v4-structural-variants" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release=grch38 \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_dbvar:
    input:
        bed=f"output/full/tracks/track-strucvars-dbvar-{{genome_release}}-{{version}}+{DV.tracks}/dbvar.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-dbvar-{{genome_release}}-{{version}}+{PV.worker}/bgdb-dbvar.bin",
        spec=f"output/full/worker/bgdb-dbvar-{{genome_release}}-{{version}}+{PV.worker}/bgdb-dbvar.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-db-var \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=dbvar \
            --value title="dbVar" \
            --value creator="NCBI dbVar Team" \
            --value source="https://www.ncbi.nlm.nih.gov/dbvar/" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_dgv:
    input:
        bed=f"output/full/tracks/track-strucvars-dgv-{{genome_release}}-{{version}}+{DV.tracks}/dgv.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-dgv-{{genome_release}}-{{version}}+{PV.worker}/bgdb-dgv.bin",
        spec=f"output/full/worker/bgdb-dgv-{{genome_release}}-{{version}}+{PV.worker}/bgdb-dgv.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-dgv \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=dgv \
            --value title="DGV" \
            --value creator="The Centre for Applied Genomics" \
            --value source="http://dgv.tcag.ca/dgv/app/home" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """


rule output_worker_bgdb_dgv_gs:
    input:
        bed=f"output/full/tracks/track-strucvars-dgv-gs-{{genome_release}}-{{version}}+{DV.tracks}/dgv-gs.bed.gz",
    output:
        bin=f"output/full/worker/bgdb-dgv-gs-{{genome_release}}-{{version}}+{PV.worker}/bgdb-dgv-gs.bin",
        spec=f"output/full/worker/bgdb-dgv-gs-{{genome_release}}-{{version}}+{PV.worker}/bgdb-dgv-gs.spec.yaml",
    shell:
        r"""
        varfish-server-worker strucvars txt-to-bin \
            --input-type strucvar-dgv-gs \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/bgdb.spec.yaml \
            \
            --value db_name=dgv-gs \
            --value title="DGV GS" \
            --value creator="The Centre for Applied Genomics" \
            --value source="http://dgv.tcag.ca/dgv/app/home" \
            \
            --value version={wildcards.version}+{PV.worker} \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec}
        """
