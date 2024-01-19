## Rules related to CTD download


rule genes_ctd_download:  # -- download CTD
    output:
        tsv="work/download/genes/ctd/{date}/CTD_diseases.tsv.gz",
    shell:
        """
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        wget -O {output.tsv} https://ctdbase.org/reports/CTD_diseases.tsv.gz
        """
