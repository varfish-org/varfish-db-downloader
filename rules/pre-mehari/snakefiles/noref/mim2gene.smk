# Obtain current dump of mim2gene_medgen file from NCBI.


rule grchXX_mim2gene_medgen_download:
    output:
        txt="noref/mim2gene/{download_date}/download/mim2gene_medgen",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.txt} \
            http://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen
        """


rule result_grchXX_mim2gene_medgen_to_tsv:
    input:
        header="header/mim2gene.txt",
        txt="noref/mim2gene/{download_date}/download/mim2gene_medgen",
    output:
        tsv="noref/mim2gene/{download_date}/Mim2geneMedgen.tsv",
        release_info="noref/mim2gene/{download_date}/Mim2geneMedgen.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.txt}
        ) \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nMim2geneMedgen\t$(date +%Y/%m/%d)\t\t-" > {output.release_info}
        """
