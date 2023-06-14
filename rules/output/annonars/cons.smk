## Rules to create annonars RocksDB for UCSC conservation track.

import os


rule output_annonars_cons:  # -- build UCSC conservation track RocksDB with annonars
    input:
        tsv="work/annos/{genome_release}/features/cons/{v_cons}/ucsc_conservation.tsv",
    output:
        "output/annonars/annos/seqvars/cons-{genome_release}-{v_cons}+{v_annonars}/rocksdb/IDENTITY",
    threads: int(os.environ.get("THREADS_ANNONARS_IMPORT", "96"))
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb_per_cpu=2000,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_cons=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars cons import \
            --path-in-tsv {input.tsv} \
            --path-out-rocksdb $(dirname {output}) \
            --genome-release {wildcards.genome_release}
        """
