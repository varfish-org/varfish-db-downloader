## Rules related to gnomAD exomes/genomes sequence variants.


# URL prefix of gnomAD downloads.
GNOMAD_PREFIX = "https://gnomad-public-us-east-1.s3.amazonaws.com/release"


rule annos_gnomad_nuclear_download_grch37:  # -- download gnomAD v2 exomes/genomes for GRCh37
    output:
        vcf="work/download/annos/grch37/seqvars/gnomad_{kind}/gnomad.{kind}.r{version}.sites.{chrom}.vcf.bgz",
        vcf_tbi="work/download/annos/grch37/seqvars/gnomad_{kind}/gnomad.{kind}.r{version}.sites.{chrom}.vcf.bgz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.vcf.bgz

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf_tbi} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.vcf.bgz.tbi
        """


rule annos_gnomad_nuclear_download_grch38_liftover_v2:  # -- download gnomAD v2 exomes lift-over for GRCh38
    output:
        vcf="work/download/annos/grch38/seqvars/gnomad_{kind}/gnomad.{kind}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz",
        vcf_tbi="work/download/annos/grch38/seqvars/gnomad_{kind}/gnomad.{kind}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/liftover_grch38/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.liftover_grch38.vcf.bgz

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf_tbi} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/liftover_grch38/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r{wildcards.version}.sites.{wildcards.chrom}.liftover_grch38.vcf.bgz.tbi
        """


rule annos_gnomad_nuclear_download_grch38_v3:  # -- download gnomAD genomes v3
    output:
        vcf="work/download/annos/grch38/seqvars/gnomad_{kind}/gnomad.{kind}.v{version}.sites.chr{chrom}.vcf.bgz",
        vcf_tbi="work/download/annos/grch38/seqvars/gnomad_{kind}/gnomad.{kind}.v{version}.sites.chr{chrom}.vcf.bgz.tbi",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.v{wildcards.version}.sites.chr{wildcards.chrom}.vcf.bgz

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf_tbi} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            {GNOMAD_PREFIX}/{wildcards.version}/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.v{wildcards.version}.sites.chr{wildcards.chrom}.vcf.bgz.tbi
        """


def input_annos_gnomad_nuclear_grch37(wildcards):
    """Input files for gnomAD exomes/genomes GRCh37."""
    chroms = list(range(1, 23)) + ["X"]
    # chrY is only available for GRCh37 genomes
    if wildcards.kind == "exomes":
        chroms.append("Y")
    tpl = "work/download/annos/grch37/seqvars/gnomad_{kind}/gnomad.{kind}.r{version}.sites.{chrom}.vcf.bgz"
    return [tpl.format(kind=wildcards.kind, version=DV.gnomad_v2, chrom=chrom) for chrom in chroms]


def input_annos_gnomad_nuclear_grch38(wildcards):
    """Input files for gnomAD exomes/genomes GRCh38."""
    chroms = list(range(1, 23)) + ["X", "Y"]
    if wildcards.kind == "exomes":
        tpl = "work/download/annos/grch38/seqvars/gnomad_{kind}/gnomad.{kind}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz"
        return [
            tpl.format(kind=wildcards.kind, version=DV.gnomad_v2, chrom=chrom) for chrom in chroms
        ]
    else:
        tpl = "work/download/annos/grch38/seqvars/gnomad_{kind}/gnomad.{kind}.v{version}.sites.chr{chrom}.vcf.bgz"
        return [
            tpl.format(kind=wildcards.kind, version=DV.gnomad_v3, chrom=chrom) for chrom in chroms
        ]


rule annos_seqvars_gnomad_nuclear_grch37:  # -- collect gnomAD exomes/genomes for GRCh37
    input:
        input_annos_gnomad_nuclear_grch37,
    output:
        touch("work/annos/grch37/seqvars/gnomad_{kind}/.done"),


rule annos_seqvars_gnomad_nuclear_grch38:  # -- collect gnomAD exomes/genomes for GRCh38
    input:
        input_annos_gnomad_nuclear_grch38,
    output:
        touch("work/annos/grch38/seqvars/gnomad_{kind}/.done"),
