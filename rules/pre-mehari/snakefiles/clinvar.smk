# Build ClinVar TSV files.


rule grch3x_clinvar_hgnc_download:
    output:
        tsv="work/pre-mehari/{genomebuild}/clinvar/{quarterly_release_date}/download/hgnc_complete_set.tsv",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.tsv} \
            https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_{wildcards.quarterly_release_date}.tsv
        """


rule grch3x_clinvar_extract:  # -- download/extract ClinVar files
    input:
        gz="work/download/annos/clinvar/clinvar-data-extract-vars-{clinvar_version}+{annonars_version}.tar.gz",
    output:
        jsonl="work/pre-mehari/{genomebuild}/clinvar/{clinvar_version}+{annonars_version}/download/clinvar-variants.jsonl.gz",
    params:
        genomebuild_sm=lambda wildcards: wildcards.genomebuild.lower(),
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        tar -C $TMPDIR -xvf {input.gz}
        cp \
            $TMPDIR/clinvar-data-extract-vars-{wildcards.clinvar_version}+{wildcards.annonars_version}/clinvar-variants-{params.genomebuild_sm}-seqvars.jsonl.gz \
            {output.jsonl}
        """


rule grch3x_clinvar_convert:
    input:
        tsv="work/pre-mehari/{genomebuild}/clinvar/{quarterly_release_date}/download/hgnc_complete_set.tsv",
        jsonl="work/pre-mehari/{genomebuild}/clinvar/{clinvar_version}+{annonars_version}/download/clinvar-variants.jsonl.gz",
    output:
        tsv="output/pre-mehari/{genomebuild}/clinvar/{quarterly_release_date}+{clinvar_version}+{annonars_version}/Clinvar.tsv",
        release_info="output/pre-mehari/{genomebuild}/clinvar/{quarterly_release_date}+{clinvar_version}+{annonars_version}/Clinvar.release_info",
    params:
        version=lambda wildcards: wildcards.clinvar_version,
    script:
        "scripts/clinvar-to-tsv.py"
