rule result_grchxx_ensembl_to_genesymbol:
    input:
        "header/ensembltogenesymbol.txt",
    output:
        tsv="{genome_build}/ensembltogenesymbol/{download_date}/EnsemblToGeneSymbol.tsv",
        release_info=(
            "{genome_build}/ensembltogenesymbol/{download_date}/EnsemblToGeneSymbol.release_info"
        ),
    shell:
        r"""
        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            prefix=grch37.
        else
            prefix=
        fi

        (
            cat {input} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            wget --no-check-certificate \
                -O- \
                http://${{prefix}}ensembl.org'/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>' \
        ) \
        | awk -F $'\t' 'BEGIN {{ OFS=FS }} ($NF != "")' \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblToGeneSymbol\t{wildcards.download_date}\t{wildcards.genome_build}\t" > {output.release_info}
        """
