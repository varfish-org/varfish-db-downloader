## Rules to create annonars RocksDB for AlphaMissense.


def input_output_annonars_alphamissense(wildcards):
    """Input function for ``rule output_annonars_alphamissense``."""
    if wildcards.genome_release == "grch37":
        genome = "hg19"
    else:
        genome = "hg38"
    return f"work/download/annos/alphamissense/1/{genome}/AlphaMissense_{genome}.tsv.gz"


rule output_annonars_alphamissense:  # -- build AlphaMissense RocksDB with annonars
    input:
        input_output_annonars_alphamissense,
    output:
        rocksdb_identity=(
            "output/full/annonars/alphamissense-{genome_release}-{v_alphamissense}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/full/annonars/alphamissense-{genome_release}-{v_alphamissense}+{v_annonars}/spec.yaml"
        ),
        manifest=(
            "output/full/annonars/alphamissense-{genome_release}-{v_alphamissense}+{v_annonars}/MANIFEST.txt"
        ),
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_alphamissense=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars tsv import \
            --path-in-tsv {input} \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            \
            --col-chrom Chrom \
            --col-start Pos \
            --col-ref Ref \
            --col-alt Alt \
            \
            --db-name AlphaMissense \
            --db-version {wildcards.v_alphamissense} \
            --genome-release {wildcards.genome_release} \
            \
            --inference-row-count 100000 \
            --skip-row-count 3 \
            --add-default-null-values

        varfish-db-downloader tpl \
            --template rules/output/annonars/alphamissense.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_alphamissense}+{wildcards.v_annonars} \
            --value v_alphamissense={wildcards.v_alphamissense} \
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
