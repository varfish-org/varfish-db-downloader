## Rules to create build worker phenotypes database.

import os


#: Number of threads to use for simulating in viguno.
VIGUNO_SIMULATE_THREADS = int(os.environ.get("VIGUNO_SIMULATE_THREADS", 96))


rule output_worker_pheno:  # -- copy HPO and simulate
    input:
        obo="work/download/hpo/{v_hpo}/hp.obo",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        genes_to_phenotype="work/download/hpo/{v_hpo}/phenotype_to_genes.txt",
    output:
        obo="output/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        rocksdb_identity="output/viguno/hpo-{v_hpo}+{v_viguno}/resnik/IDENTITY",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    threads: VIGUNO_SIMULATE_THREADS
    resources:
        mem_mb_per_cpu="2GB",
    shell:
        r"""
        cp -a {input.obo} {output.obo}
        cp -a {input.hpoa} {output.hpoa}

        awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $3, $4, $1, $2, $6 }}' \
            {input.genes_to_phenotype} \
        > {output.phenotype_to_genes}

        viguno simulate \
            --ic-base gene \
            --similarity fun-sim-avg \
            --combiner resnik \
            --path-hpo-dir $(dirname {input.obo}) \
            --path-out-rocksdb $(dirname {output.rocksdb_identity})/scores-fun-sim-avg-resnik-gene \
            --min-terms 1 \
            $(if [[ "{RUNS_IN_CI}" == "True" ]]; then \
                echo --max-terms 1; \
                echo --num-simulations 10; \
                echo --only-gene ARID1B; \
            else \
                echo --max-terms 10; \
            fi) \
            --seed 42
        """


rule global_hpo_to_bin:  # -- convert to .bin
    input:
        obo="output/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        rocksdb_identity="output/viguno/hpo-{v_hpo}+{v_viguno}/resnik/IDENTITY",
    output:
        obo="output/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    shell:
        r"""
        viguno convert \
            --path-hpo-dir $(dirname {input.obo}) \
            --path-out-bin {output.bin}
        """
