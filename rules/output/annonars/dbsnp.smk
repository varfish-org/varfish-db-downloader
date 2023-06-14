## Rules to create annonars RocksDB for dbSNP.

import os


rule output_annonars_dbsnp:  # -- build dbSNP RocksDB with annonars
    input:
        vcf="work/download/annos/{genome_release}/seqvars/dbsnp/{v_dbsnp}/dbsnp.vcf.gz",
    output:
        rocksdb_identity=(
            "output/annonars/dbsnp-{genome_release}-{v_dbsnp}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=("output/annonars/dbsnp-{genome_release}-{v_dbsnp}+{v_annonars}/spec.yaml"),
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
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --genome-release {wildcards.genome_release}

        varfish-db-downloader tpl \
            --template rules/output/annonars/dbsnp.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_dbsnp}+{wildcards.v_annonars} \
            --value v_dbsnp={wildcards.v_dbsnp} \
            \
            --value v_annonars={wildcards.v_annonars} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
