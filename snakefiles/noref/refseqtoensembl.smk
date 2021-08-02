rule noref_refseq_to_ensembl_download:
    output:
        "noref/refseqtoensembl/{download_date}/download/gene2ensembl.gz",
    shell:
        r"""
        wget \
            -O {output} \
            'http://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz'
        """


rule result_noref_refseq_to_ensembl_tsv:
    input:
        header="header/refseqtoensembl.txt",
        gz="noref/refseqtoensembl/{download_date}/download/gene2ensembl.gz",
    output:
        tsv="noref/refseqtoensembl/{download_date}/RefseqToEnsembl.tsv",
        release_info="noref/refseqtoensembl/{download_date}/RefseqToEnsembl.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            gunzip -c {input.gz} \
            | awk -F $"\t" 'BEGIN{{OFS=FS}}$1=="9606"{{split($5,a,".");print $2,$3,a[1]}}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToEnsembl\t$(date +%Y/%m/%d)\t\t-" > {output.release_info}
        """
