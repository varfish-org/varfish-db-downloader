# Obtain current dump of HGNC gene information info.


rule grch37_hgnc_download:
    output:
        txt="GRCh37/hgnc/{download_date}/download/hgnc_complete_set.txt",
    shell:
        r"""
        wget \
            -O {output.txt} \
                        http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
        """


rule grch37_refseq_to_hgnc_download:
    output:
        "GRCh37/hgnc/{download_date}/download/ref_GRCh37.p13_top_level.gff3.gz",
    shell:
        r"""
        wget \
            -O {output} \
            http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz \             
        """


rule result_grch37_hgnc_to_tsv:
    input:
        header="header/hgnc.txt",
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        txt="GRCh37/hgnc/{download_date}/download/hgnc_complete_set.txt",
    output:
        tsv="GRCh37/hgnc/{download_date}/Hgnc.tsv",
        release_info="GRCh37/hgnc/{download_date}/Hgnc.release_info",
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

        echo -e "table\tversion\tgenomebuild\tnull_value\nHgnc\t$(date +%Y/%m/%d)\tGRCh37\t" > {output.release_info}
        """


rule result_grch37_refseq_to_hgnc_to_tsv:
    input:
        header="header/refseqtohgnc.txt",
        gff="GRCh37/hgnc/{download_date}/download/ref_GRCh37.p13_top_level.gff3.gz",
    output:
        tsv="GRCh37/hgnc/{download_date}/RefseqToHgnc.tsv",
        release_info="GRCh37/hgnc/{download_date}/RefseqToHgnc.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.gff} \
            | perl -ne 'if (/GeneID:([^,;]+)[,;]HGNC:([^,;]+)[,;]/) {{print $1,"\tHGNC:",$2,"\n"}}' \
            | sort -u
        ) \
        > {output.tsv}
        echo -e "table\tversion\tgenomebuild\tnull_value\nRefseqToHgnc\t$(date +%Y/%m/%d)\tGRCh37\t" > {output.release_info}
        """
