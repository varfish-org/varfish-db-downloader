## Rules to create build viguno phenotypes database subsets (dev/exomes).
#
# We will copy the full HPO (text and binary) but reduce the simulation count.


rule subset_viguno_pheno_exomes:  # -- create exomes subset
    input:
        obo="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/full/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
    output:
        obo="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/reduced-exomes/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
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
        obo="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/full/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
    output:
        obo="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
        bin="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml="output/reduced-dev/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml",
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
        """
