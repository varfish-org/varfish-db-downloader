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
        rocksdb_identity=(
            "output/full/annonars/cadd-{genome_release}-{v_cadd}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=("output/full/annonars/cadd-{genome_release}-{v_cadd}+{v_annonars}/spec.yaml"),
        manifest=("output/full/annonars/cadd-{genome_release}-{v_cadd}+{v_annonars}/MANIFEST.txt"),
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_cadd=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        if [ "${{CI:-false}}" = "true" ]; then
            echo "Skipping annonars CADD import in CI environment."
            mkdir -p $(dirname {output.rocksdb_identity})
            touch {output.rocksdb_identity} {output.spec_yaml} {output.manifest}
            exit 0
        fi
        annonars tsv import \
            --path-in-tsv {input.indels} \
            --path-in-tsv {input.snvs} \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
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

        varfish-db-downloader tpl \
            --template rules/output/annonars/cadd.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_cadd}+{wildcards.v_annonars} \
            --value v_cadd={wildcards.v_cadd} \
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
