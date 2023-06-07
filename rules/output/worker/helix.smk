## Rules to create annonars RocksDB for HelixMtDb.

import os


rule output_worker_helixmtdb:  # -- build HelixMtDb RocksDB with annonars
    input:
        vcf="work/annos/{genome_release}/seqvars/helixmtdb/{v_helixmtdb}/helixmtdb.vcf.gz",
    output:
        "output/worker/annos/seqvars/helixmtdb-{genome_release}-{v_helixmtdb}+{v_annonars}/rocksdb/IDENTITY",
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_helixmtdb=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars helixmtdb import \
            --path-in-vcf {input.vcf} \
            --path-out-rocksdb $(dirname {output}) \
            --genome-release {wildcards.genome_release}
        """
