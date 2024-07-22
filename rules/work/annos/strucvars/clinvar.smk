## Download of clinvar files for worker.


rule annos_strucvars_clinvar_download:  # -- download/extract ClinVar files
    output:
        jsonl="work/download/annos/{genome_release}/strucvars/clinvar/{clinvar_version}/clinvar-variants-{genome_release}-strucvars.jsonl.gz",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_version=RE_VERSION_MULTI,
    shell:
        r"""
        clinvar_version=$(echo "{wildcards.clinvar_version}" | sed -e 's/-//g' | cut -d '+' -f 1)

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        wget --no-check-certificate \
            -O /tmp/clinvar-data-extract-vars-{wildcards.clinvar_version}.tar.gz \
            https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-$clinvar_version/clinvar-data-extract-vars-{wildcards.clinvar_version}.tar.gz

        tar -C $TMPDIR -xvf /tmp/clinvar-data-extract-vars-{wildcards.clinvar_version}.tar.gz
        cp $TMPDIR/clinvar-data-extract-vars-{wildcards.clinvar_version}/clinvar-variants-{wildcards.genome_release}-strucvars.jsonl.gz {output.jsonl}
        """
