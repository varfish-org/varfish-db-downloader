## Download of clinvar files for worker.


rule annos_clinvar_download:  # -- download/extract ClinVar files
    output:
        "work/download/annos/clinvar/clinvar-data-extract-vars-{clinvar_version}+{annonars_version}.tar.gz",
    wildcard_constraints:
        clinvar_version=RE_VERSION_MULTI,
    shell:
        r"""
        wget --no-check-certificate \
            -O {output} \
            https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-{wildcards.clinvar_version}/clinvar-data-extract-vars-{wildcards.clinvar_version}+{wildcards.annonars_version}.tar.gz
        """
