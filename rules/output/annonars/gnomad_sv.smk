## Rules to create annonars RocksDB for gnomAD-SV.


import os


rule output_annonars_gnomad_sv_grch37_exac:  # -- build gnomAD-SV RocksDB with annonars; ExAC CNV variant
    input:
        bed="work/download/annos/grch37/strucvars/exac/0.3.1/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
    output:
        rocksdb_identity=(
            "output/full/annonars/gnomad-sv-exomes-grch37-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY",
        ),
        spec_yaml=("output/full/annonars/gnomad-sv-exomes-grch37-{v_gnomad}+{v_annonars}/spec.yaml"),
        manifest=(
            "output/full/annonars/gnomad-sv-exomes-grch37-{v_gnomad}+{v_annonars}/MANIFEST.txt"
        ),
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        v_gnomad=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars gnomad-sv import \
            --gnomad-kind exomes \
            --genome-release grch37 \
            --path-in-vcf {input.bed} \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --gnomad-version 1.0

        varfish-db-downloader tpl \
            --template rules/output/annonars/gnomad_sv_exomes_grch37.spec.yml \
            --value today={TODAY} \
            --value genome_release=grch37 \
            \
            --value version={wildcards.v_gnomad}+{wildcards.v_annonars} \
            --value v_gnomad={wildcards.v_gnomad} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """


rule output_annonars_gnomad_sv_grch37_gnomad_sv2:  # -- build gnomAD-SV RocksDB with annonars; gnomAD-SV v2 variant
    input:
        vcf=[
            f"work/download/annos/grch37/strucvars/gnomad/2.1.1/gnomad_v2.1_sv.{token}.vcf.gz"
            for token in ("sites", "controls_only.sites", "nonneuro.sites")
        ],
    output:
        rocksdb_identity=(
            "output/full/annonars/gnomad-sv-genomes-grch37-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY",
        ),
        spec_yaml=(
            "output/full/annonars/gnomad-sv-genomes-grch37-{v_gnomad}+{v_annonars}/spec.yaml"
        ),
        manifest=(
            "output/full/annonars/gnomad-sv-genomes-grch37-{v_gnomad}+{v_annonars}/MANIFEST.txt"
        ),
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        v_gnomad=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars gnomad-sv import \
            --gnomad-kind genomes \
            --genome-release grch37 \
            $(for vcf in {input.vcf}; do echo --path-in-vcf $vcf; done) \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --gnomad-version {wildcards.v_gnomad}

        varfish-db-downloader tpl \
            --template rules/output/annonars/gnomad_sv_genomes_grch37.spec.yml \
            --value today={TODAY} \
            --value genome_release=grch37 \
            \
            --value version={wildcards.v_gnomad}+{wildcards.v_annonars} \
            --value v_gnomad={wildcards.v_gnomad} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """


rule output_annonars_gnomad_sv_grch38_gnomad_cnv4:  # -- build gnomAD-SV RocksDB with annonars; gnomAD-CNV v4 variant
    input:
        vcf=[
            f"work/download/annos/grch38/strucvars/gnomad_cnv/{{v_gnomad}}/gnomad.v{{v_gnomad}}.cnv.{token}.vcf.gz"
            for token in ("all", "non_neuro", "non_neuro_controls")
        ],
    output:
        rocksdb_identity=(
            "output/full/annonars/gnomad-sv-exomes-grch38-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY",
        ),
        spec_yaml=("output/full/annonars/gnomad-sv-exomes-grch38-{v_gnomad}+{v_annonars}/spec.yaml"),
        manifest=(
            "output/full/annonars/gnomad-sv-exomes-grch38-{v_gnomad}+{v_annonars}/MANIFEST.txt"
        ),
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        v_gnomad=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars gnomad-sv import \
            --gnomad-kind exomes \
            --genome-release grch38 \
            $(for vcf in {input.vcf}; do echo --path-in-vcf $vcf; done) \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --gnomad-version {wildcards.v_gnomad}

        varfish-db-downloader tpl \
            --template rules/output/annonars/gnomad_sv_exomes_grch38.spec.yml \
            --value today={TODAY} \
            --value genome_release=grch38 \
            \
            --value version={wildcards.v_gnomad}+{wildcards.v_annonars} \
            --value v_gnomad={wildcards.v_gnomad} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """


rule output_annonars_gnomad_sv_grch38_gnomad_sv4:  # -- build gnomAD-SV RocksDB with annonars; gnomAD-SV v4 variant
    input:
        vcf="work/download/annos/grch38/strucvars/gnomad_sv/{v_gnomad}/gnomad.v{v_gnomad}.sv.sites.vcf.gz",
    output:
        rocksdb_identity=(
            "output/full/annonars/gnomad-sv-genomes-grch38-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY",
        ),
        spec_yaml=(
            "output/full/annonars/gnomad-sv-genomes-grch38-{v_gnomad}+{v_annonars}/spec.yaml"
        ),
        manifest=(
            "output/full/annonars/gnomad-sv-genomes-grch38-{v_gnomad}+{v_annonars}/MANIFEST.txt"
        ),
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_gnomad=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars gnomad-sv import \
            --gnomad-kind genomes \
            --genome-release grch38 \
            --path-in-vcf {input.vcf} \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --gnomad-version {wildcards.v_gnomad}

        varfish-db-downloader tpl \
            --template rules/output/annonars/gnomad_sv_genomes_grch38.spec.yml \
            --value today={TODAY} \
            --value genome_release=grch38 \
            \
            --value version={wildcards.v_gnomad}+{wildcards.v_annonars} \
            --value v_gnomad={wildcards.v_gnomad} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """
