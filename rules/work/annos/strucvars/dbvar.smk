## Rules related to structural variants in dbVar.


rule annos_strucvars_dbvar_download:  # -- download dbVar files
    output:
        "work/download/annos/{genome_release}/strucvars/dbvar/{version}/{genome_release_nolower}.nr_{sv_type}.tsv.gz",
    shell:
        r"""
        # Ensure that the version is the latest in the release notes.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        wget --no-check-certificate \
            -O $TMPDIR/release_notes_listing \
            https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/release_notes

        grep NR_stats.2 $TMPDIR/release_notes_listing \
        | tail -n 1 \
        | grep "$(echo {wildcards.version} | sed 's/-//g')" \
        || {{ >&2 echo "Version {wildcards.version} is not latest in dbVar release notes"; exit 1; }}

        # Actual download.
        wget --no-check-certificate \
            -O {output} \
            https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/{wildcards.sv_type}/{wildcards.genome_release_nolower}.nr_{wildcards.sv_type}.tsv.gz
        """


def input_annos_strucvars_dbvar_process(wildcards):
    """Input functions for rule ``rule annos_strucvars_dbvar_process``."""
    mapping = {
        "grch37": "GRCh37",
        "grch38": "GRCh38",
    }
    tpl = (
        "work/download/annos/{genome_release}/strucvars/dbvar/{version}/"
        "{genome_release_nolower}.nr_{typ}.tsv.gz"
    )
    return {
        "download": [
            tpl.format(
                genome_release=wildcards.genome_release,
                version=wildcards.version,
                genome_release_nolower=mapping[wildcards.genome_release],
                typ=typ,
            )
            for typ in ("deletions", "duplications", "insertions")
        ]
    }


rule annos_strucvars_dbvar_process:  # -- process dbVar files
    input:
        unpack(input_annos_strucvars_dbvar_process),
    output:
        bed="output/full/worker/annos/strucvars/dbvar-{genome_release}-{version}/dbvar.bed.gz",
        bed_md5="output/full/worker/annos/strucvars/dbvar-{genome_release}-{version}/dbvar.bed.gz.md5",
        bed_tbi="output/full/worker/annos/strucvars/dbvar-{genome_release}-{version}/dbvar.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/annos/strucvars/dbvar-{genome_release}-{version}/dbvar.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/vardbs-strucvar-dbvar.awk \
            <(zcat {input.download}) \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
