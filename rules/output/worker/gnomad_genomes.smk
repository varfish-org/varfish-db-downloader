## Rules to create annonars RocksDB for gnomAD-gnomes.

import os


rule output_worker_gnomad_genomes:  # -- build gnomAD-genomes RocksDB with annonars
    input:
        vcf="work/download/annos/{genome_release}/seqvars/gnomad_genomes/{v_gnomad}/.done",
    output:
        "output/worker/annos/seqvars/gnomad-genomes-{genome_release}-{v_gnomad}+{v_annonars}/rocksdb/IDENTITY",
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
            --path-out-rocksdb $(dirname {output}) \
            --gnomad-kind genomes \
            --genome-release {wildcards.genome_release} \
            --gnomad-version {wildcards.v_gnomad}
        """
