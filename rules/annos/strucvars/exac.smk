## Rules related to ExAC-CNV.


rule annos_strucvars_exac_grch37_download:  # -- process ExAC-CNV files
    output:
        bed="work/download/annos/grch37/strucvars/exac/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.bed} \
            ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed
        """


rule annos_strucvars_exac_grch37_process:  # -- process ExAC-CNV files
    input:
        bed="work/download/annos/grch37/strucvars/exac/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
    output:
        bed="work/annos/grch37/strucvars/exac/exac.bed.gz",
        bed_md5="work/annos/grch37/strucvars/exac/exac.bed.gz.md5",
        bed_tbi="work/annos/grch37/strucvars/exac/exac.bed.gz.tbi",
        bed_tbi_md5="work/annos/grch37/strucvars/exac/exac.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -f scripts/vardbs-grch37-strucvar-exac.awk \
            {input.bed} \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
