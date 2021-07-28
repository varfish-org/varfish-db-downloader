# Obtain mapping of Ensembl Gene ID - Ensembl Transcript ID - Entrez ID
# If you are just interested in Ensembl Gene Id and Entrez ID this will result in duplicate rows.


rule grchxx_ensembl_to_refseq_download:
    output:
        "{genome_build}/ensembltorefseq/{download_date}/download/ensembl_to_refseq.tsv",
    shell:
        r"""
        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            prefix=grch37.
        else
            prefix=
        fi

        wget \
            "http://${prefix}ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?> <!DOCTYPE Query> <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" > <Dataset name = "hsapiens_gene_ensembl" interface = "default" > <Attribute name = "ensembl_gene_id" /> <Attribute name = "ensembl_transcript_id" /> <Attribute name = "entrezgene" /> </Dataset> </Query>" \
            -O {output}
        """


rule result_grchxx_ensembl_to_refseq_tsv:
    input:
        header="header/ensembltorefseq.txt",
        tsv="{genome_build}/ensembltorefseq/{download_date}/download/ensembl_to_refseq.tsv",
    output:
        tsv="{genome_build}/ensembltorefseq/{download_date}/EnsemblToRefseq.tsv",
        release_info="{genome_build}/ensembltorefseq/{download_date}/EnsemblToRefseq.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            cat {input.tsv}
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nEnsemblToRefseq\t$(date +%Y/%m/%d)\t{wildcards.genome_build}\t" > {output.release_info}
        """
