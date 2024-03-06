## Convert gene regions to binary for worker.


def input_genes_regions_worker_convert(wildcards):
    """Input function for rule genes_regions_worker_convert."""
    return {
        "bed": f"work/annos/{wildcards.genome_release}/features/{wildcards.source}/{wildcards.version}/{wildcards.source}_genes.bed.gz",
    }


rule genes_regions_worker_convert:
    input:
        unpack(input_genes_regions_worker_convert),
    output:
        bin=f"output/full/worker/genes-regions-{{genome_release}}-{{version}}+{PV.worker}/{{source}}_genes.bin",
        spec_yaml=f"output/full/worker/genes-regions-{{genome_release}}-{{version}}+{PV.worker}/{{source}}_genes.spec.yaml",
    shell:
        r"""
        varfish-server-worker \
            db to-bin \
            --input-type gene-region \
            --path-input {input.bed} \
            --path-output {output.bin}

        varfish-db-downloader tpl \
            --template rules/output/worker/{wildcards.source}_genes.spec.yaml \
            --value today={TODAY} \
            --value genome_release={wildcards.genome_release} \
            \
            --value version={wildcards.version} \
            \
            --value v_worker={PV.worker} \
            --value v_downloader={PV.downloader} \
        > {output.spec_yaml}
        """
