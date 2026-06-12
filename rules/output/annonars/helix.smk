## Rules to create annonars RocksDB for HelixMtDb.

import os


rule output_annonars_helixmtdb:  # -- build HelixMtDb RocksDB with annonars
    input:
        vcf="work/annos/{genome_release}/seqvars/helixmtdb/{v_helixmtdb}/helixmtdb.vcf.gz",
        validate_script="scripts/validate_rocksdb.sh",
    output:
        rocksdb_dir=directory(
            "output/full/annonars/helixmtdb-{genome_release}-{v_helixmtdb}+{v_annonars}/rocksdb"
        ),
        rocksdb_identity=(
            "output/full/annonars/helixmtdb-{genome_release}-{v_helixmtdb}+{v_annonars}/rocksdb/IDENTITY",
        ),
        spec_yaml=(
            "output/full/annonars/helixmtdb-{genome_release}-{v_helixmtdb}+{v_annonars}/spec.yaml",
        ),
        manifest=(
            "output/full/annonars/helixmtdb-{genome_release}-{v_helixmtdb}+{v_annonars}/MANIFEST.txt",
        ),
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_helixmtdb=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        annonars helixmtdb import \
            --path-in-vcf {input.vcf} \
            --path-out-rocksdb {output.rocksdb_dir} \
            --genome-release {wildcards.genome_release}

        bash {input.validate_script} "{output.rocksdb_dir}"

        varfish-db-downloader tpl \
            --template rules/output/annonars/helix.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_helixmtdb}+{wildcards.v_annonars} \
            --value v_helixmtdb={wildcards.v_helixmtdb} \
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
