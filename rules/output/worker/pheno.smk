## Rules to create build worker phenotypes database.


rule output_worker_pheno:  # -- build genes protobuf file
    input:
        obo="work/download/hpo/{v_hpo}/hp.obo",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        phenotype_to_genes="work/download/hpo/{v_hpo}/phenotype_to_genes.txt",
    output:
        obo="output/worker/pheno-{v_hpo}+{v_worker}/hp.obo",
        hpoa="output/worker/pheno-{v_hpo}+{v_worker}/phenotype.hpoa",
        phenotype_to_genes="output/worker/pheno-{v_hpo}+{v_worker}/phenotype_to_genes.txt",
        rocksdb_identity="output/worker/pheno-{v_hpo}+{v_worker}/resnik/IDENTITY",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_worker=RE_VERSION,
    shell:
        r"""
        cp -a {input.obo} {output.obo}
        cp -a {input.hpoa} {output.hpoa}
        cp -a {input.phenotype_to_genes} {output.phenotype_to_genes}

        varfish-server-worker pheno prepare \
            --path-hpo-dir $(dirname {input.obo}) \
            --path-out-rocksdb $(dirname {output.rocksdb_identity}) \
            --min-terms 10 \
            --max-terms 10 \
            --seed 42 \
            --num-simulations 10
        """
