## Rules related to Thousand Genomes SVs.


rule annos_strucvars_g1k_grch37_download:  # -- download Thousand Genomes SVs
    output:
        vcf="work/download/annos/grch37/strucvars/g1k/phase3-v2/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.vcf} \
            https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz
        """


rule annos_strucvars_g1k_grch37_process:  # -- process Thousand Genomes SVs
    input:
        vcf="work/download/annos/grch37/strucvars/g1k/phase3-v2/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
    output:
        bed=f"output/full/tracks/track-strucvars-g1k-grch37-phase3-v2+{DV.tracks}/g1k.bed.gz",
        bed_md5=f"output/full/tracks/track-strucvars-g1k-grch37-phase3-v2+{DV.tracks}/g1k.bed.gz.md5",
        bed_tbi=f"output/full/tracks/track-strucvars-g1k-grch37-phase3-v2+{DV.tracks}/g1k.bed.gz.tbi",
        bed_tbi_md5=f"output/full/tracks/track-strucvars-g1k-grch37-phase3-v2+{DV.tracks}/g1k.bed.gz.tbi.md5",
    shell:
        r"""
        zcat {input.vcf} \
        | awk \
            -F $'\t' \
            -f scripts/vardbs-grch37-strucvar-g1k.awk \
        | sort-bed - \
        | bgzip -c \
        > {output.bed}

        tabix -p bed -S 1 -f {output.bed}

        md5sum {output.bed} >{output.bed_md5}
        md5sum {output.bed_tbi} >{output.bed_tbi_md5}
        """
