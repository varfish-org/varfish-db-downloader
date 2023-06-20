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
                f"reduced-{wildcards.set_name}/targets/{wildcards.genome_release}/"
                f"refseq/{refseq_version}/refseq_target_exoms.bed"
            ),
            "rocksdb_identity": (
                f"output/annonars/{wildcards.name}-{wildcards.genome_release}-"
                f"{wildcards.version_multi}/rocksdb/IDENTITY"
            ),
        }
    return result


rule subset_annonars:  # -- create exomes subset
    input:
        unpack(input_subset_annonars),
    output:
        rocksdb_identity="reduced-{set_name}/annonars/{name}-{genome_release}-{version_multi}/rocksdb/IDENTITY",
    wildcard_constraints:
        name=RE_NAME,
        genome_release=RE_GENOME,
        v_hpo=RE_VERSION,
        versions=RE_VERSION_MULTI,
    shell:
        r"""
        annonars db-utils copy \
            --path-in $(dirname {input.rocksdb_identity}) \
            --path-out $(dirname {output.rocksdb_identity}) \
            --path-beds {input.bed}
        """
