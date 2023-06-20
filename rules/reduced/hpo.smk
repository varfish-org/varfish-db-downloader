## Rules to create build viguno phenotypes database subsets (dev/exomes).
#
# We will copy the full HPO (text and binary) but reduce the simulation count.


rule subset_viguno_pheno_exomes:  # -- create exomes subset
    input:
        obo="output/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        rocksdb_identity="output/viguno/hpo-{v_hpo}+{v_viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        bin="output/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/viguno/hpo-{v_hpo}-{v_viguno}/spec.yaml",
    output:
        obo="reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        rocksdb_identity="reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        bin="reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="reduced-exomes/viguno/hpo-{v_hpo}-{v_viguno}/spec.yaml",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    shell:
        r"""
        cp -a {input.obo} {output.obo}
        cp -a {input.hpoa} {output.hpoa}
        cp -a {input.phenotype_to_genes} {output.phenotype_to_genes}
        cp -a {input.bin} {output.bin}
        cp -a {input.spec_yaml} {output.spec_yaml}

        cp -ar $(dirname {input.rocksdb_identity})/. $(dirname {output.rocksdb_identity})/.
        """


rule subset_worker_pheno_dev:  # -- create development subset
    input:
        obo="output/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        rocksdb_identity="output/viguno/hpo-{v_hpo}+{v_viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        bin="output/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/viguno/hpo-{v_hpo}-{v_viguno}/spec.yaml",
    output:
        obo="reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        rocksdb_identity="reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/scores-fun-sim-avg-resnik-gene/IDENTITY",
        bin="reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="reduced-dev/viguno/hpo-{v_hpo}-{v_viguno}/spec.yaml",
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
        cp -a {input.phenotype_to_genes} {output.phenotype_to_genes}
        cp -a {input.bin} {output.bin}
        cp -a {input.spec_yaml} {output.spec_yaml}

        viguno simulate \
            --ic-base gene \
            --similarity resnik \
            --combiner fun-sim-avg \
            --path-hpo-dir $(dirname {input.obo}) \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --min-terms 1 \
            $(if [[ "{RUNS_IN_CI}" == "True" ]]; then \
                echo --max-terms 1; \
                echo --num-simulations 10; \
                echo --only-gene ARID1B; \
            else \
                echo --max-terms 10; \
                echo --num-simulations 100; \
            fi) \
            --seed 42
        """
