## Rules to create annonars RocksDB for dbSNP.

import os


rule output_annonars_dbsnp:  # -- build dbSNP RocksDB with annonars
    input:
        vcf="work/download/annos/{genome_release}/seqvars/dbsnp/{v_dbsnp}/dbsnp.vcf.gz",
    output:
        "output/annonars/dbsnp-{genome_release}-{v_dbsnp}+{v_annonars}/rocksdb/IDENTITY",
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_dbsnp=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars dbsnp import \
            --path-in-vcf {input.vcf} \
            --path-out-rocksdb $(dirname {output}) \
            --genome-release {wildcards.genome_release}
        """
