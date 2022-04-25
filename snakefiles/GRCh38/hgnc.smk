# Obtain current dump of HGNC gene information info.


rule grch38_hgnc_download:
    output:
        txt="GRCh38/hgnc/{download_date}/download/hgnc_complete_set.txt",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.txt} \
            http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
        """


rule grch38_refseq_to_hgnc_download:
    output:
        "GRCh38/hgnc/{download_date}/download/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output} \
            http://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
        """


rule result_grch38_hgnc_to_tsv:
    input:
        header="header/hgnc.txt",
        reference="GRCh38/reference/hs38/hs38.fa",
        txt="GRCh38/hgnc/{download_date}/download/hgnc_complete_set.txt",
    output:
        tsv="GRCh38/hgnc/{download_date}/Hgnc.tsv",
        release_info="GRCh38/hgnc/{download_date}/Hgnc.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.txt} \
            | awk -F $'\t' '
                BEGIN {{
                    OFS = FS
                }}
                {{
                    $NF = $NF "\t" substr($22, 1, length($22) - 2)
                    print
                }}'
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHgnc\t$(date +%Y/%m/%d)\tGRCh38\t" > {output.release_info}
        """


rule result_grch38_refseq_to_hgnc_to_tsv:
    input:
        header="header/refseqtohgnc.txt",
        gff="GRCh38/hgnc/{download_date}/download/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
    output:
        tsv="GRCh38/hgnc/{download_date}/RefseqToHgnc.tsv",
        release_info="GRCh38/hgnc/{download_date}/RefseqToHgnc.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.gff} \
            | perl -ne 'if (/GeneID:([^,;]+)[,;]HGNC:([^,;]+)[,;]/) {{print $1,"\tHGNC:",$2,"\n"}}' \
            | sed -e 's/HGNC:HGNC:/HGNC:/g' \
            | sort -u
        ) \
        > {output.tsv}
        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToHgnc\t$(date +%Y/%m/%d)\tGRCh38\t" > {output.release_info}
        """
