## Rules to create build annonars functional annotation database..


rule work_annonars_functional_download:
    output:
        "work/download/refseq/{genomebuild}/{version}/{assembly}_genomic.gff.gz",
    params:
        url_version=lambda wildcards: (
            f"{DV.refseq_ref_38}-{DV.refseq_38}"
            if DV.refseq_38.startswith("RS_")
            else f"{DV.refseq_38}/{DV.refseq_ref_38_assembly}"
        )
        if wildcards.genomebuild == "grch38"
        else (
            f"{DV.refseq_ref_37}-{DV.refseq_37}"
            if DV.refseq_37.startswith("RS_")
            else f"{DV.refseq_37}/{DV.refseq_ref_37_assembly}"
        ),
    shell:
        r"""
        wget -O {output} \
            {DV.refseq_base_url}/{params.url_version}/{wildcards.assembly}_genomic.gff.gz
        """


def output_annonars_functional_input(wildcards):
    if wildcards.genome_release == "grch37":
        return (
            f"work/download/refseq/grch37/{DV.refseq_37}/{DV.refseq_ref_37_assembly}_genomic.gff.gz"
        )
    else:
        return (
            f"work/download/refseq/grch38/{DV.refseq_38}/{DV.refseq_ref_38_assembly}_genomic.gff.gz"
        )


rule output_annonars_functional:  # -- build annonars functional RocksDB file
    input:
        output_annonars_functional_input,
        validate_script="scripts/validate_rocksdb.sh",
    output:
        rocksdb_dir=directory("output/full/annonars/functional-{genome_release}-{v_refseq}+{v_annonars}/rocksdb"),
        rocksdb_identity=(
            "output/full/annonars/functional-{genome_release}-{v_refseq}+{v_annonars}/"
            "rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/full/annonars/functional-{genome_release}-{v_refseq}+{v_annonars}/spec.yaml"
        ),
        manifest=(
            "output/full/annonars/functional-{genome_release}-{v_refseq}+{v_annonars}/MANIFEST.txt"
        ),
    wildcard_constraints:
        v_refseq=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        zgrep '^#\|RefSeqFE' {input} > $TMPDIR/tmp.gff

        annonars functional import -vvv \
            --genome-release {wildcards.genome_release} \
            --path-in-gff $TMPDIR/tmp.gff \
            --path-out-rocksdb {output.rocksdb_dir}

        bash {input.validate_script} "{output.rocksdb_dir}"

        varfish-db-downloader tpl \
            --template rules/output/annonars/functional.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.v_refseq}+{wildcards.v_annonars} \
            --value v_refseq={wildcards.v_refseq} \
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
