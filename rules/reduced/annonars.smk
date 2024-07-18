## Rules to create build annonars database subsets (dev/exomes).
#
# We will copy the full HPO (text and binary) but reduce the simulation count.


def input_subset_annonars(wildcards):
    """Input function for ``rule subset_annonars``."""
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
            f"output/full/annonars/{wildcards.name}-{wildcards.genome_release}-"
            f"{wildcards.version_multi}/rocksdb/IDENTITY"
        ),
        "spec_yaml": (
            f"output/full/annonars/{wildcards.name}-{wildcards.genome_release}-"
            f"{wildcards.version_multi}/spec.yaml"
        ),
    }
    return result


rule subset_annonars:  # -- create exomes subset
    input:
        unpack(input_subset_annonars),
    output:
        rocksdb_identity="output/reduced-{set_name}/annonars/{name}-{genome_release}-{version_multi}/rocksdb/IDENTITY",
        spec_yaml="output/reduced-{set_name}/annonars/{name}-{genome_release}-{version_multi}/spec.yaml",
    wildcard_constraints:
        name=RE_NAME,
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
            --skip-cfs dbsnp_by_rsid \
            --skip-cfs clinvar_by_accession \
            --path-in $(dirname {input.rocksdb_identity}) \
            --path-out $(dirname {output.rocksdb_identity}) \
            --path-beds {input.bed}

        cp {input.spec_yaml} {output.spec_yaml}
        """
