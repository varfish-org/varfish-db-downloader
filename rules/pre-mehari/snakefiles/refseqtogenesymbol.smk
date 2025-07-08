rule noref_refseq_to_genesymbol_dl:
    output:
        "work/download/pre-mehari/noref/refseqtogenesymbol/{download_date}/gene2accession.gz",
    shell:
        r"""
        wget --no-check-certificate -O {output} http://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
        """


rule result_noref_refseq_to_genesymbol_tsv:
    input:
        gz="work/download/pre-mehari/noref/refseqtogenesymbol/{download_date}/gene2accession.gz",
    output:
        tsv="output/pre-mehari/noref/refseqtogenesymbol/{download_date}/RefseqToGeneSymbol.tsv",
        release_info="output/pre-mehari/noref/refseqtogenesymbol/{download_date}/RefseqToGeneSymbol.release_info",
    shell:
        r"""
        (
            echo -e "entrez_id\tgene_symbol"
            zcat {input.gz} \
            | awk -F $'\t' '
                BEGIN {{
                    OFS=FS
                }}
                $1 == 9606 {{
                    print $2,$16
                }}' \
            | sort -u
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToGeneSymbol\t{wildcards.download_date}\t\t" > {output.release_info}
        """
