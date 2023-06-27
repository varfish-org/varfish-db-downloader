## Rules related to structural variants from ClinVar.


rule annos_strucvars_clinvar_download:  # -- download/extract ClinVar files
    output:
        tar=f"work/download/annos/{{genome_release}}/strucvars/clinvar/{{clinvar_version}}/clinvar-strucvar-{{genome_release}}-{DV.clinvar_release}.tar.gz",
        tsv="work/download/annos/{genome_release}/strucvars/clinvar/{clinvar_version}/clinvar_strucvar.tsv.gz",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_version=RE_VERSION,
    shell:
        r"""
        clinvar_version=$(echo "{wildcards.clinvar_version}" | sed -e 's/-//g' | cut -d '+' -f 1)

        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" ERR EXIT

        wget --no-check-certificate \
            -O {output.tar} \
            https://github.com/bihealth/annonars-data-clinvar/releases/download/clinvar-weekly-$clinvar_version/$(basename {output.tar})

        if [[ {wildcards.genome_release} == grch37 ]]; then
            release=GRCh37
        else
            release=GRCh38
        fi

        tar -C $TMPDIR -xvf $(readlink -f {output.tar})
        cut -f 2- $TMPDIR/clinvar-strucvar-{wildcards.genome_release}*/output.tsv \
        | awk \
            -F $'\t' \
            -v release=$release \
            -v clinvar_version={wildcards.clinvar_version} \
            \
            '
            BEGIN {{ OFS=FS }}
            (NR == 1) {{ print "#" $0; }}
            (NR > 1) {{ $8 = clinvar_version; print $0; }}
            ' \
        | gzip -c \
        > {output.tsv}
        """


def input_annos_srucvar_clinvar_convert(wildcards):
    """Return input files for ``rule annos_strucvar_clinvar_convert``."""
    clinvar_version = wildcards.clinvar_release.replace("-", "").split("+")[0]
    return {
        "tsv": (
            f"work/download/annos/{wildcards.genome_release}/strucvars/"
            f"clinvar/{clinvar_version}/clinvar_strucvar.tsv.gz"
        )
    }


rule annos_strucvar_clinvar_convert:
    input:
        unpack(input_annos_srucvar_clinvar_convert),
    output:
        bin="output/full/worker/annos/strucvars/clinvar-{genome_release}-{clinvar_release}/clinvar_strucvars.bin",
    wildcard_constraints:
        genome_release=RE_GENOME,
        clinvar_release=RE_VERSION,
    shell:
        r"""
        varfish-server-worker db to-bin \
            --input-type clinvar-sv \
            --path-input {input.tsv} \
            --path-output-bin {output.bin}
        """


# def input_annos_strucvars_dgv_process(wildcards):
#     """Input function for ``rule annos_strucvars_dgv_process``."""
#     mapping = {
#         "grch37": {
#             "genome_release_nolower": "GRCh37",
#             "genome_release_ucsc": "hg19",
#         },
#         "grch38": {
#             "genome_release_nolower": "GRCh38",
#             "genome_release_ucsc": "hg38",
#         },
#     }
#     tpl = (
#         "work/download/annos/{genome_release}/strucvars/dgv/{version}/"
#         "{genome_release_nolower}_{genome_release_ucsc}_variants_{version_dashes}.txt"
#     )
#     return {
#         "txt": tpl.format(
#             version_dashes="-".join([wildcards.version[0:4], wildcards.version[4:6], wildcards.version[6:8]]),
#             **mapping[wildcards.genome_release],
#             **wildcards,
#         )
#     }


# rule annos_strucvars_dgv_process:  # -- download DGV files
#     input:
#         unpack(input_annos_strucvars_dgv_process),
#     output:
#         bed="output/full/worker/annos/strucvars/dgv-{genome_release}-{version}/dgv.bed.gz",
#         bed_md5="output/full/worker/annos/strucvars/dgv-{genome_release}-{version}/dgv.bed.gz.md5",
#         bed_tbi="output/full/worker/annos/strucvars/dgv-{genome_release}-{version}/dgv.bed.gz.tbi",
#         bed_tbi_md5="output/full/worker/annos/strucvars/dgv-{genome_release}-{version}/dgv.bed.gz.tbi.md5",
#     shell:
#         r"""
#         awk \
#             -F $'\t' \
#             -f scripts/vardbs-strucvar-dgv.awk \
#             {input.txt} \
#         | grep -v _gl \
#         | sort-bed - \
#         | bgzip -c \
#         > {output.bed}

#         tabix -p bed -S 1 -f {output.bed}

#         md5sum {output.bed} >{output.bed_md5}
#         md5sum {output.bed_tbi} >{output.bed_tbi_md5}
#         """


# rule annos_strucvars_dgv_gs_download:  # -- download DGV GS files
#     output:
#         gff3="work/download/annos/{genome_release}/strucvars/dgv_gs/{version}/{filename}",
#     shell:
#         r"""
#         wget --no-check-certificate \
#             -O {output.gff3} \
#             http://dgv.tcag.ca/dgv/docs/{wildcards.filename}
#         """


# def input_annos_strucvars_dgv_gs_process(wildcards):
#     """Input function for ``rule annos_strucvars_dgv_gs_process``."""
#     mapping = {
#         "grch37": "DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
#         "grch38": "DGV.GS.hg38.gff3",
#     }
#     return {
#         "gff3": (
#             "work/download/annos/{genome_release}/strucvars/dgv_gs/{version}/{filename}"
#         ).format(filename=mapping[wildcards.genome_release], **wildcards)
#     }
# rule annos_strucvars_dgv_gs_process:  # -- download DGV GS files
#     input:
#         unpack(input_annos_strucvars_dgv_gs_process),
#     output:
#         bed="output/full/worker/annos/strucvars/dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz",
#         bed_md5="output/full/worker/annos/strucvars/dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz.md5",
#         bed_tbi="output/full/worker/annos/strucvars/dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz.tbi",
#         bed_tbi_md5="output/full/worker/annos/strucvars/dgv-gs-{genome_release}-{version}/dgv-gs.bed.gz.tbi.md5",
#     shell:
#         r"""
#         awk \
#             -F $'\t' \
#             -f scripts/vardbs-strucvar-dgv_gs.awk \
#             {input.gff3} \
#         | grep -v _gl \
#         | sort-bed - \
#         | bgzip -c \
#         > {output.bed}
#         tabix -p bed -S 1 -f {output.bed}
#         md5sum {output.bed} >{output.bed_md5}
#         md5sum {output.bed_tbi} >{output.bed_tbi_md5}
#         """
