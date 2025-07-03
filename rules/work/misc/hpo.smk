## Rules related to HPO download.


rule global_hpo_download:  # -- download HPO files
    output:
        obo="work/download/hpo/{v_hpo}/hp.obo",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        phenotype_to_genes="work/download/hpo/{v_hpo}/phenotype_to_genes.txt",
        genes_to_phenotype="work/download/hpo/{v_hpo}/genes_to_phenotype.txt",
    shell:
        r"""
        wget -O {output.obo} https://github.com/obophenotype/human-phenotype-ontology/releases/download/{wildcards.v_hpo}/hp.obo
        wget -O {output.hpoa} https://github.com/obophenotype/human-phenotype-ontology/releases/download/{wildcards.v_hpo}/phenotype.hpoa
        wget -O {output.phenotype_to_genes} https://github.com/obophenotype/human-phenotype-ontology/releases/download/{wildcards.v_hpo}/phenotype_to_genes.txt
        wget -O {output.genes_to_phenotype} https://github.com/obophenotype/human-phenotype-ontology/releases/download/{wildcards.v_hpo}/genes_to_phenotype.txt
        """
