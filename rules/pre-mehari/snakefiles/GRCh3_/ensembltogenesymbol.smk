
rule grchxx_ensembl_to_genesymbol_download:
    output:
        tsv="output/pre-mehari/{genomebuild}/ensembltogenesymbol/{ensembl}/EnsemblToGeneSymbol.tsv",
        release_info="output/pre-mehari/{genomebuild}/ensembltogenesymbol/{ensembl}/EnsemblToGeneSymbol.release_info",
    params:
        ensembl_archive_url=lambda wildcards: DV.ensembl_38_archive_url if wildcards.genomebuild == "GRCh38" else DV.ensembl_37_archive_url,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)

        echo -e "ensembl_gene_id\tgene_symbol" > {output.tsv}
        wget --no-check-certificate \
            -O $TMPDIR/tmp \
            '{params.ensembl_archive_url}/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>'
        awk -F $'\t' 'BEGIN {{ OFS=FS }} ($NF != "")' $TMPDIR/tmp | sort -u >> {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblToGeneSymbol\t{wildcards.ensembl}\t{wildcards.genomebuild}\t" > {output.release_info}
        """
