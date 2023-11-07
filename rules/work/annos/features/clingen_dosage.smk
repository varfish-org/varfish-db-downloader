## Work rules related to ClinGen dosage sensitivity regions.


rule annos_features_clingen_dosage_download:  # -- download ClinGen dosage sensitivity
    output:
        tsv="work/download/annos/{genome_release}/features/clingen_dosage/{date}/clingen_region_curation_list.tsv",
        tsv_md5="work/download/annos/{genome_release}/features/clingen_dosage/{date}/clingen_region_curation_list.tsv.md5",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        if [[ "{wildcards.genome_release}" == "grch37" ]]; then
            URL_RELEASE=GRCh37
        else
            URL_RELEASE=GRCh38
        fi

        wget --no-check-certificate \
            -O {output.tsv} \
            ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_${{URL_RELEASE}}.tsv

        md5sum {output.tsv} > {output.tsv_md5}
        """
