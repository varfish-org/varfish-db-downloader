## Rules related to structural variants in dbVar.


rule annos_strucvars_dbvar_grch37_download:  # -- download dbVar files
    output:
        "work/download/annos/grch37/strucvars/dbvar/GRCh37.nr_{sv_type}.tsv.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output} \
            https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/{wildcards.sv_type}/GRCh37.nr_{wildcards.sv_type}.tsv.gz
        """


rule annos_strucvars_dbvar_grch37_process:  # -- process dbVar files
    input:
        download=expand(
            "work/download/annos/grch37/strucvars/dbvar/GRCh37.nr_{type}.tsv.gz",
            type=["deletions", "duplications", "insertions"],
        ),
    output:
        bed="work/annos/grch37/strucvars/dbvar/dbvar.bed.gz",
        bed_md5="work/annos/grch37/strucvars/dbvar/dbvar.bed.gz.md5",
        bed_tbi="work/annos/grch37/strucvars/dbvar/dbvar.bed.gz.tbi",
        bed_tbi_md5="vardbs/grch37/strucvars/dbvar/dbvar.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-dbvar.awk \
            <(zcat {input.download}) \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
