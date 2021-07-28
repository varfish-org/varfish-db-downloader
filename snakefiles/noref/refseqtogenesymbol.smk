rule noref_refseq_to_genesymbol_dl:
    output:
        "noref/refseqtogenesymbol/latest/download/gene2accession.gz",
    shell:
        r"""
        wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
        """


rule noref_refseq_to_genesymbol_tsv:
    input:
        header="header/refseqtogenesymbol.txt",
        gz="noref/refseqtogenesymbol/latest/download/gene2accession.gz",
    output:
        tsv="noref/refseqtogenesymbol/latest/RefseqToGeneSymbol.tsv",
        release_info="noref/refseqtogenesymbol/latest/RefseqToGeneSymbol.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
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

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToGeneSymbol\tlatest\t{wildcards.genome_build}\t" > {output.release_info}
        """
