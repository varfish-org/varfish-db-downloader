## Rules related to CTD download


rule genes_ctd_download:  # -- download CTD
    output:
        tsv="work/download/genes/ctd/{date}/CTD_diseases.tsv.gz",
    shell:
        """
        wget -O {output.tsv} https://ctdbase.org/reports/CTD_diseases.tsv.gz
        """
