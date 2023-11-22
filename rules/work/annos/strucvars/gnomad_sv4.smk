## Rules related to gnomAD-SV and gnomAD-CNV v4.


rule annos_strucvars_gnomad_sv_4_grch38_download:  # -- download gnomAD-SV 4.0 files
    output:
        vcf="work/download/annos/grch38/strucvars/gnomad_sv/4.0/gnomad.v4.0.sv.chr{chrom}.vcf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/genome_sv/gnomad.v4.0.sv.chr{wildcards.chrom}.vcf.gz
        """


rule annos_strucvars_gnomad_cnv_4_grch38_download:  # -- download gnomAD-CNV 4.0 files
    output:
        vcf="work/download/annos/grch38/strucvars/gnomad_cnv/4.0/gnomad.v4.0.cnv.{token}.vcf.gz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/exome_cnv/gnomad.v4.0.cnv.{wildcards.token}.vcf.gz
        """
