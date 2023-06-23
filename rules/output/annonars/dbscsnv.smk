## Rules to create annonars RocksDB for dbscSNV.

import os


def input_output_annonars_dbscsnv(wildcards):
    """Input function for ``rule output_annonars_dbscsnv``."""
    return [entry for entry in var_tsv_files_dbscsnv(wildcards.v_dbscsnv)]


rule output_annonars_dbscsnv:  # -- build dbscSNV RocksDB with annonars
    input:
        input_output_annonars_dbscsnv,
    output:
        rocksdb_identity=(
            "output/full/annonars/dbscsnv-{genome_release}-{v_dbscsnv}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=("output/full/annonars/dbscsnv-{genome_release}-{v_dbscsnv}+{v_annonars}/spec.yaml"),
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_dbscsnv=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars tsv import \
            --db-name dbscSNV \
            --db-version {wildcards.v_dbscsnv} \
            --genome-release {wildcards.genome_release} \
            --null-values=. \
            --inference-row-count 100000 \
            --path-schema-json rules/output/annonars/dbscsnv-schema.json \
            \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            \
            $(if [[ "{wildcards.genome_release}" == "grch37" ]]; then \
                echo --col-chrom 'chr'; \
                echo --col-start 'pos'; \
            else \
                echo --col-chrom 'hg38_chr'; \
                echo --col-start 'hg38_pos'; \
            fi) \
            --col-ref 'ref' \
            --col-alt 'alt' \
            \
            $(for path in {input}; do echo --path-in-tsv $path; done)

        varfish-db-downloader tpl \
            --template rules/output/annonars/dbscsnv.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_dbscsnv}+{wildcards.v_annonars} \
            --value v_dbscsnv={wildcards.v_dbscsnv} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
