## Rules to create annonars RocksDB for dbNSFP.

import os


def input_output_annonars_dbnsfp(wildcards):
    """Input function for ``rule output_annonars_dbnsfp``."""
    return [entry.format(version=wildcards.v_dbnsfp, variant="") for entry in var_tsv_dbnsfp()]


rule output_annonars_dbnsfp:  # -- build dbNSFP RocksDB with annonars
    input:
        input_output_annonars_dbnsfp,
    output:
        rocksdb_identity=(
            "output/full/annonars/dbnsfp-{genome_release}-{v_dbnsfp}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=("output/full/annonars/dbnsfp-{genome_release}-{v_dbnsfp}+{v_annonars}/spec.yaml"),
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_dbnsfp=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars tsv import \
            --db-name dbNSFP \
            --db-version {wildcards.v_dbnsfp} \
            --genome-release {wildcards.genome_release} \
            --null-values=. \
            --inference-row-count 100000 \
            \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --path-schema-json rules/output/annonars/dbnsfp-schema-{wildcards.v_dbnsfp}.json \
            \
            $(if [[ "{wildcards.genome_release}" == "grch37" ]]; then \
                echo --col-chrom 'hg19_chr'; \
                echo --col-start 'hg19_pos(1-based)'; \
            else \
                echo --col-chrom 'chr'; \
                echo --col-start 'pos(1-based)'; \
            fi) \
            --col-ref 'ref' \
            --col-alt 'alt' \
            \
            $(for path in {input}; do echo --path-in-tsv $path; done)

        varfish-db-downloader tpl \
            --template rules/output/annonars/dbnsfp.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_dbnsfp}+{wildcards.v_annonars} \
            --value v_dbnsfp={wildcards.v_dbnsfp} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
