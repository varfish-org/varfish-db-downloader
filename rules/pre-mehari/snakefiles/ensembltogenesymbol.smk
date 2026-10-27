rule grchxx_ensembl_to_genesymbol_download:
    input:
        download_gtf=lambda wildcards: (
            f"work/genes/ensembl/{wildcards.genomebuild.lower()}/{wildcards.ensembl}/download/"
            f"Homo_sapiens.{wildcards.genomebuild}.{wildcards.ensembl}.chr_patch_hapl_scaff.gtf.gz"
        ),
    output:
        tsv="output/pre-mehari/{genomebuild}/ensembltogenesymbol/{ensembl}/EnsemblToGeneSymbol.tsv",
        release_info="output/pre-mehari/{genomebuild}/ensembltogenesymbol/{ensembl}/EnsemblToGeneSymbol.release_info",
    shell:
        r"""
        echo -e "ensembl_gene_id\tgene_symbol" > {output.tsv}
        awk \
            -F $'\t' \
            -f scripts/genes-ensembl-genesymbol.awk \
            <(zcat {input.download_gtf}) \
        | sort -u \
        >> {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblToGeneSymbol\t{wildcards.ensembl}\t{wildcards.genomebuild}\t" > {output.release_info}
        """
