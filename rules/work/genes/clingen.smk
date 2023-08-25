## Rules related to ClinGen curation download.


rule genes_clingen_download:  # -- download ClinGen curations
    output:
        csv="work/download/genes/clingen/{date}/clingen.csv",
        csv_md5="work/download/genes/clingen/{date}/clingen.csv.md5",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        wget --no-check-certificate \
            -O {output.csv} \
            https://search.clinicalgenome.org/kb/reports/curation-activity-summary-report

        md5sum {output.csv} > {output.csv_md5}
        """
