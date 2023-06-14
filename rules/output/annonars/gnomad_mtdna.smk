## Rules to create annonars RocksDB for gnomAD-mtDNA.

import os


rule output_annonars_gnomad_mtdna:  # -- build gnomAD-mtDNA RocksDB with annonars
    input:
        vcf="work/annos/{genome_release}/seqvars/gnomad_mtdna/{v_gnomad}/gnomad_mtdna.vcf.gz",
    output:
        "output/annonars/gnomad-mtdna-{genome_release}-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY",
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
        annonars gnomad-mtdna import \
            --path-in-vcf {input.vcf} \
            --path-out-rocksdb $(dirname {output}) \
            --genome-release {wildcards.genome_release} \
            --gnomad-version {wildcards.v_gnomad}
        """
