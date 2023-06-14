## Rules to create annonars RocksDB for UCSC conservation track.

import os


rule output_annonars_cons:  # -- build UCSC conservation track RocksDB with annonars
    input:
        tsv="work/annos/{genome_release}/features/cons/{v_cons}/ucsc_conservation.tsv",
    output:
        rocksdb_identity=(
            "output/annonars/cons-{genome_release}-{v_cons}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=("output/annonars/cons-{genome_release}-{v_cons}+{v_annonars}/spec.yaml"),
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
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --genome-release {wildcards.genome_release}

        varfish-db-downloader tpl \
            --template rules/output/annonars/cons.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_cons}+{wildcards.v_annonars} \
            --value v_cons={wildcards.v_cons} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
