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
    }
    return result


rule subset_mehari:  # -- create exomes subset
    input:
        unpack(input_subset_mehari),
    output:
        rocksdb_identity="output/reduced-{set_name}/mehari/freqs-{genome_release}-{version_multi}/rocksdb/IDENTITY",
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_hpo=RE_VERSION,
        versions=RE_VERSION_MULTI,
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    shell:
        r"""
        annonars db-utils copy \
            --path-in $(dirname {input.rocksdb_identity}) \
            --path-out $(dirname {output.rocksdb_identity}) \
            --path-beds {input.bed}
        """
