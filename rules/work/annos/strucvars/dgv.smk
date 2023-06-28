## Rules related to structural variants in DGV and DGV-GS.


rule annos_strucvars_dgv_download:  # -- download DGV files
    output:
        txt="work/download/annos/{genome_release}/strucvars/dgv/{version}/{filename}",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.txt} \
            http://dgv.tcag.ca/dgv/docs/{wildcards.filename}
        """


def input_annos_strucvars_dgv_process(wildcards):
    """Input function for ``rule annos_strucvars_dgv_process``."""
    mapping = {
        "grch37": {
            "genome_release_nolower": "GRCh37",
            "genome_release_ucsc": "hg19",
        },
        "grch38": {
            "genome_release_nolower": "GRCh38",
            "genome_release_ucsc": "hg38",
        },
    }
    tpl = (
        "work/download/annos/{genome_release}/strucvars/dgv/{version}/"
        "{genome_release_nolower}_{genome_release_ucsc}_variants_{version_dashes}.txt"
    )
    return {
        "txt": tpl.format(
            version_dashes="-".join(
                [wildcards.version[0:4], wildcards.version[4:6], wildcards.version[6:8]]
            ),
            **mapping[wildcards.genome_release],
            **wildcards,
        )
    }


rule annos_strucvars_dgv_process:  # -- download DGV files
    input:
        unpack(input_annos_strucvars_dgv_process),
    output:
        bed="output/full/worker/track-strucvars-dgv-{genome_release}-{version}/dgv.bed.gz",
        bed_md5="output/full/worker/track-strucvars-dgv-{genome_release}-{version}/dgv.bed.gz.md5",
        bed_tbi="output/full/worker/track-strucvars-dgv-{genome_release}-{version}/dgv.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/track-strucvars-dgv-{genome_release}-{version}/dgv.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/vardbs-strucvar-dgv.awk \
            {input.txt} \
        | grep -v _gl \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule annos_strucvars_dgv_gs_download:  # -- download DGV GS files
    output:
        gff3="work/download/annos/{genome_release}/strucvars/dgv_gs/{version}/{filename}",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.gff3} \
            http://dgv.tcag.ca/dgv/docs/{wildcards.filename}
        """


def input_annos_strucvars_dgv_gs_process(wildcards):
    """Input function for ``rule annos_strucvars_dgv_gs_process``."""
    mapping = {
        "grch37": "DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
        "grch38": "DGV.GS.hg38.gff3",
    }
    return {
        "gff3": (
            "work/download/annos/{genome_release}/strucvars/dgv_gs/{version}/{filename}"
        ).format(filename=mapping[wildcards.genome_release], **wildcards)
    }


rule annos_strucvars_dgv_gs_process:  # -- download DGV GS files
    input:
        unpack(input_annos_strucvars_dgv_gs_process),
    output:
        bed="output/full/worker/track-strucvars-dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz",
        bed_md5="output/full/worker/track-strucvars-dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz.md5",
        bed_tbi="output/full/worker/track-strucvars-dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz.tbi",
        bed_tbi_md5="output/full/worker/track-strucvars-dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/vardbs-strucvar-dgv_gs.awk \
            {input.gff3} \
        | grep -v _gl \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
