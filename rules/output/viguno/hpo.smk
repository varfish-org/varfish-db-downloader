## Rules to create build viguno phenotypes database.

import os



rule output_viguno_pheno:  # -- copy HPO
    input:
        obo="work/download/hpo/{v_hpo}/hp.obo",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        genes_to_phenotype="work/download/hpo/{v_hpo}/phenotype_to_genes.txt",
    output:
        obo="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hp.obo",
        hpoa="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype.hpoa",
        phenotype_to_genes="output/full/viguno/hpo-{v_hpo}+{v_viguno}/phenotype_to_genes.txt",
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    threads: THREADS
    resources:
        mem_mb=MEMORY,
    shell:
        r"""
        cp -a {input.obo} {output.obo}
        cp -a {input.hpoa} {output.hpoa}

        awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $3, $4, $1, $2, $6 }}' \
            {input.genes_to_phenotype} \
        > {output.phenotype_to_genes}
        """


rule global_hpo_to_bin:  # -- convert to .bin
    input:
        obo="work/download/hpo/{v_hpo}/hp.obo",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        genes_to_phenotype="work/download/hpo/{v_hpo}/phenotype_to_genes.txt",
    output:
        bin="output/full/viguno/hpo-{v_hpo}+{v_viguno}/hpo.bin",
        spec_yaml=("output/full/viguno/hpo-{v_hpo}+{v_viguno}/spec.yaml"),
    wildcard_constraints:
        v_hpo=RE_VERSION,
        v_viguno=RE_VERSION,
    shell:
        r"""
        viguno convert \
            --path-hpo-dir $(dirname {input.obo}) \
            --path-out-bin {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/viguno/hpo.spec.yaml \
            --value today={TODAY} \
            \
            --value version={wildcards.v_hpo}+{wildcards.v_viguno} \
            --value v_hpo={wildcards.v_hpo} \
            \
            --value v_viguno={wildcards.v_viguno} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
