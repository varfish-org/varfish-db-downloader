## Rules to create annonars RocksDB for CADD.

import os


def input_output_annonars_cadd(wildcards):
    """Input function for ``rule output_annonars_cadd``."""
    if wildcards.genome_release == "grch37":
        return {
            "indels": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "InDels_inclAnno.tsv.gz"
            ),
            "indels_tbi": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "InDels_inclAnno.tsv.gz.tbi"
            ),
            "snvs": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "whole_genome_SNVs_inclAnno.tsv.gz"
            ),
            "snvs_tbi": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "whole_genome_SNVs_inclAnno.tsv.gz.tbi"
            ),
        }
    else:
        return {
            "indels": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "gnomad.genomes.r3.0.indel_inclAnno.tsv.gz"
            ),
            "indels_tbi": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "gnomad.genomes.r3.0.indel_inclAnno.tsv.gz.tbi"
            ),
            "snvs": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "whole_genome_SNVs_inclAnno.tsv.gz"
            ),
            "snvs_tbi": (
                f"work/download/annos/{wildcards.genome_release}/seqvars/cadd/{wildcards.v_cadd}/"
                "whole_genome_SNVs_inclAnno.tsv.gz.tbi"
            ),
        }


rule output_annonars_cadd:  # -- build CADD RocksDB with annonars
    input:
        unpack(input_output_annonars_cadd),
    output:
        "output/annonars/cadd-{genome_release}-{v_cadd}+{v_annonars}/rocksdb/IDENTITY",
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_cadd=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars tsv import \
            --path-in-tsv {input.indels} \
            --path-in-tsv {input.snvs} \
            --path-out-rocksdb $(dirname {output}) \
            \
            --col-chrom Chrom \
            --col-start Pos \
            --col-ref Ref \
            --col-alt Alt \
            \
            --db-name CADD \
            --db-version {wildcards.v_cadd} \
            --genome-release {wildcards.genome_release} \
            \
            --inference-row-count 100000 \
            --skip-row-count 1 \
            --add-default-null-values \
            --path-schema-json rules/output/annonars/cadd-schema-{wildcards.genome_release}.json
        """
