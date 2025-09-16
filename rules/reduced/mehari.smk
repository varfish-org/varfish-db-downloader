## Rules to create build mehari database subsets (dev/exomes).
#
# We will copy the full HPO (text and binary) but reduce the simulation count.


def input_subset_mehari(wildcards):
    """Input function for ``rule subset_mehari``."""
    if wildcards.genome_release == "grch37":
        refseq_version = DV.refseq_37
    else:
        refseq_version = DV.refseq_38
    result = {
        "bed": (
            f"output/reduced-{wildcards.set_name}/targets/{wildcards.genome_release}/"
            f"refseq/{refseq_version}/refseq_target_exons.bed"
        ),
        "rocksdb_identity": (
            f"output/full/mehari/freqs-{wildcards.genome_release}-"
            f"{wildcards.version_multi}/rocksdb/IDENTITY"
        ),
        "spec_yaml": (
            f"output/full/mehari/freqs-{wildcards.genome_release}-"
            f"{wildcards.version_multi}/spec.yaml"
        ),
    }
    return result


rule subset_mehari:  # -- create exomes subset
    input:
        unpack(input_subset_mehari),
    output:
        rocksdb_identity="output/reduced-{set_name}/mehari/freqs-{genome_release}-{version_multi}/rocksdb/IDENTITY",
        spec_yaml="output/reduced-{set_name}/mehari/freqs-{genome_release}-{version_multi}/spec.yaml",
        manifest="output/reduced-{set_name}/mehari/freqs-{genome_release}-{version_multi}/MANIFEST.txt",
    wildcard_constraints:
        genome_release=RE_GENOME,
        versions=RE_VERSION_MULTI,
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    shell:
        r"""
         if [ "${{CI:-false}}" = "true" ]; then
            echo "Skipping subset mehari CI environment."
            mkdir -p $(dirname {output.rocksdb_identity})
            touch {output.rocksdb_identity} {output.spec_yaml} {output.manifest}
            exit 0
        fi

        annonars db-utils copy \
            --path-in $(dirname {input.rocksdb_identity}) \
            --path-out $(dirname {output.rocksdb_identity}) \
            --path-beds {input.bed}

        cp {input.spec_yaml} {output.spec_yaml}

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """
