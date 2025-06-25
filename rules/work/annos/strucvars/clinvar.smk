## Download of clinvar files for worker.


rule annos_strucvars_clinvar_extract:  # -- download/extract ClinVar files
    input:
        gz="work/download/annos/clinvar/clinvar-data-extract-vars-{clinvar_version}+{annonars_version}.tar.gz",
    output:
        jsonl="work/download/annos/{genome_release}/strucvars/clinvar/{clinvar_version}+{annonars_version}/clinvar-variants-{genome_release}-strucvars.jsonl.gz",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_version=RE_VERSION_MULTI,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        tar -C $TMPDIR -xvf {input.gz}
        cp \
            $TMPDIR/clinvar-data-extract-vars-{wildcards.clinvar_version}+{wildcards.annonars_version}/clinvar-variants-{wildcards.genome_release}-strucvars.jsonl.gz \
            {output.jsonl}
        """
