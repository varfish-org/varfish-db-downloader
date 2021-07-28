rule noref_refseq_to_ensembl_download:
    output:
        "{noref}/refseqtoensembl/latest/download/gene2ensembl.gz",
    shell:
        r"""
        wget \
            'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz' \
            -O {output}
        """


rule noref_refseq_to_ensembl_tsv:
    input:
        header="header/refseqtoensembl.txt",
        gz="{noref}/refseqtoensembl/latest/download/gene2ensembl.gz",
    output:
        tsv="{noref}/refseqtoensembl/latest/RefseqToEnsembl.tsv",
        release_info="{noref}/refseqtoensembl/latest/RefseqToEnsembl.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            gunzip -c {input.gz} \
            | awk -F $"\t" 'BEGIN{{OFS=FS}}$1=="9606"{{split($5,a,".");print $2,$3,a[1]}}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToEnsembl\t$(date +%Y/%m/%d)\t{wildcards.genome_build}\t-" > {output.release_info}
        """
