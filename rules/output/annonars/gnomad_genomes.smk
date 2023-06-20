## Rules to create annonars RocksDB for gnomAD-genomes.

import os


rule output_annonars_gnomad_genomes:  # -- build gnomAD-genomes RocksDB with annonars
    input:
        vcf="work/download/annos/{genome_release}/seqvars/gnomad_genomes/{v_gnomad}/.done",
    output:
        rocksdb_identity=(
            "output/annonars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY"
        ),
        spec_yaml=(
            "output/annonars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/spec.yaml"
        ),
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
        annonars gnomad-nuclear import \
            $(for path in $(dirname {input.vcf})/*.bgz; do \
                echo --path-in-vcf $path; \
            done) \
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
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --gnomad-kind genomes \
            --genome-release {wildcards.genome_release} \
            --gnomad-version {wildcards.v_gnomad}

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
        """