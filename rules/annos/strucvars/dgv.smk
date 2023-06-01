## Rules related to structural variants in DGV and DGV-GS.


rule annos_strucvars_dgv_grch37_download:  # -- download DGV files
    output:
        txt="work/download/annos/grch37/strucvars/dgv/GRCh37_hg19_variants_2020-02-25.txt",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.txt} \
            http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt
        """


rule annos_strucvars_dgv_grch37_process:  # -- download DGV files
    input:
        txt="work/download/annos/grch37/strucvars/dgv/GRCh37_hg19_variants_2020-02-25.txt",
    output:
        bed="work/annos/grch37/strucvars/dgv/dgv.bed.gz",
        bed_md5="work/annos/grch37/strucvars/dgv/dgv.bed.gz.md5",
        bed_tbi="work/annos/grch37/strucvars/dgv/dgv.bed.gz.tbi",
        bed_tbi_md5="work/annos/grch37/strucvars/dgv/dgv.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-dgv.awk \
            {input.txt} \
        | grep -v _gl \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """


rule annos_strucvars_dgv_gs_grch37_download:  # -- download DGV GS files
    output:
        gff3="work/download/annos/grch37/strucvars/dgv_gs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.gff3} \
            http://dgv.tcag.ca/dgv/docs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3
        """


rule annos_strucvars_dgv_gs_grch37_process:  # -- download DGV GS files
    input:
        gff3="work/download/annos/grch37/strucvars/dgv_gs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3",
    output:
        bed="work/annos/grch37/strucvars/dgv_gs/dgv_gs.bed.gz",
        bed_md5="work/annos/grch37/strucvars/dgv_gs/dgv_gs.bed.gz.md5",
        bed_tbi="work/annos/grch37/strucvars/dgv_gs/dgv_gs.bed.gz.tbi",
        bed_tbi_md5="work/annos/grch37/strucvars/dgv_gs/dgv_gs.bed.gz.tbi.md5",
    shell:
        r"""
        awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-dgv_gs.awk \
            {input.gff3} \
        | grep -v _gl \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
