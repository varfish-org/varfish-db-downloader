## Rules related to HPO download.


rule global_hpo_download:  # -- download HPO files
    output:
        obo="work/download/hpo/{v_hpo}/hp.obo",
        hpoa="work/download/hpo/{v_hpo}/phenotype.hpoa",
        phenotype_to_genes="work/download/hpo/{v_hpo}/phenotype_to_genes.txt",
    shell:
        r"""
        tag=v$(echo {wildcards.v_hpo} | cut -b 1-4)-$(echo {wildcards.v_hpo} | cut -b 5-6)-$(echo {wildcards.v_hpo} | cut -b 7-8)

        wget -O {output.obo} https://github.com/obophenotype/human-phenotype-ontology/releases/download/$tag/hp.obo
        wget -O {output.hpoa} https://github.com/obophenotype/human-phenotype-ontology/releases/download/$tag/phenotype.hpoa
        wget -O {output.phenotype_to_genes} https://github.com/obophenotype/human-phenotype-ontology/releases/download/$tag/phenotype_to_genes.txt
        """
