## Rules to create annonars RocksDB for gnomAD-genomes.

import os


def input_gnomad_genomes(wildcards):
    return expand(
        "work/download/annos/{genome_release}/seqvars/gnomad_genomes/{v}/gnomad.genomes.v{v}.sites.chr{c}.vcf.bgz",
        v=gnomad_versions[wildcards.genome_release],
        c=CHROMS,
    )


def input_gnomad_genomes_tbi(wildcards):
    return expand(
        "work/download/annos/{genome_release}/seqvars/gnomad_genomes/{v}/gnomad.genomes.v{v}.sites.chr{c}.vcf.bgz.tbi",
        v=gnomad_versions[wildcards.genome_release],
        c=CHROMS,
    )


rule output_annonars_gnomad_genomes:  # -- build gnomAD-genomes RocksDB with annonars
    input:
        vcf=input_gnomad_genomes,
        tbi=input_gnomad_genomes_tbi,
        validate_script="scripts/validate_rocksdb.sh",
    output:
        rocksdb_dir=directory(
            "output/full/annonars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/rocksdb"
        ),
        rocksdb_identity=(
            "output/full/annonars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/full/annonars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/spec.yaml"
        ),
        manifest=(
            "output/full/annonars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/MANIFEST.txt"
        ),
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    wildcard_constraints:
        genome_release=RE_GENOME,
        v_gnomad=RE_VERSION,
        v_annonars=RE_VERSION,
    shell:
        r"""
        if [[ "${{CI:-false}}" == "true" ]]; then
            echo "Skipping gnomad in CI environment."
            mkdir -p {output.rocksdb_dir}
            touch {output.rocksdb_identity} {output.spec_yaml} {output.manifest}
            exit 0
        fi

        annonars gnomad-nuclear import \
            $(for file in {input.vcf}; do echo --path-in-vcf $file; done) \
            --import-fields-json '{{
                "vep": true,
                "var_info": true,
                "effect_info": true,
                "global_cohort_pops": true,
                "all_cohorts": true,
                "rf_info": false,
                "quality": true,
                "age_hists": true,
                "depth_details": false,
                "liftover": false
            }}' \
            --path-out-rocksdb {output.rocksdb_dir} \
            --gnomad-kind genomes \
            --genome-release {wildcards.genome_release} \
            --gnomad-version {wildcards.v_gnomad}

        bash {input.validate_script} "{output.rocksdb_dir}"

        varfish-db-downloader tpl \
            --template rules/output/annonars/gnomad_genomes.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.v_gnomad}+{wildcards.v_annonars} \
            --value v_gnomad={wildcards.v_gnomad} \
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
