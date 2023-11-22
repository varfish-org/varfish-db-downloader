## Rules related to ClinGen gene dosage pathogenicity annotation.


rule clingen_gene_download:  # -- download files
    output:
        tsv="work/genes/clingen/{date}/ClinGen_gene_curation_list_{genome_release}.tsv",
        tsv_md5="work/genes/clingen/{date}/ClinGen_gene_curation_list_{genome_release}.tsv.md5",
    shell:
        r"""
        if [[ "$(date +%Y%m%d)" != "{wildcards.date}" ]] && [[ "{FORCE_TODAY}" != "True" ]]; then
            >&2 echo "{wildcards.date} is not today"
            exit 1
        fi

        wget -O {output.tsv} \
            ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_{wildcards.genome_release}.tsv

        md5sum {output.tsv} > {output.tsv}.md5
        """
