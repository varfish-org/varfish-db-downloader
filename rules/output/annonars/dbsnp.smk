## Rules to create annonars RocksDB for dbSNP.

import os


rule output_annonars_dbsnp:  # -- build dbSNP RocksDB with annonars
    input:
        vcf="work/download/annos/{genome_release}/seqvars/dbsnp/{v_dbsnp}/dbsnp.vcf.gz",
    output:
        rocksdb_identity=(
            "output/full/annonars/dbsnp-{genome_release}-{v_dbsnp}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=("output/full/annonars/dbsnp-{genome_release}-{v_dbsnp}+{v_annonars}/spec.yaml"),
        manifest=("output/full/annonars/dbsnp-{genome_release}-{v_dbsnp}+{v_annonars}/MANIFEST.txt"),
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
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

        export TMPDIR=$(mktemp -d)
        pushd $(dirname {output.spec_yaml})
        rm -f MANIFEST.txt
        hashdeep -l -r . >$TMPDIR/MANIFEST.txt
        CHECKSUM=$(sha256sum $TMPDIR/MANIFEST.txt | cut -d ' ' -f 1)
        echo "## EOF SHA256=$CHECKSUM" >> $TMPDIR/MANIFEST.txt
        cp $TMPDIR/MANIFEST.txt MANIFEST.txt
        popd
        """
