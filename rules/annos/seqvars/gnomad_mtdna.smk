## Rules related to gnomAD mtDNA.


rule annos_gnomad_mtdna_download:  # -- download gnomAD mtDNA
    output:
        dl="work/download/annos/{genome_build}/seqvars/gnomad_mtdna/{version}/gnomad.genomes.v{version}.sites.chrM.vcf.bgz",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.dl} \
            https://datasetgnomad.blob.core.windows.net/dataset/release/{wildcards.version}/vcf/genomes/gnomad.genomes.v{wildcards.version}.sites.chrM.vcf.bgz
        """


rule annos_gnomad_mtdna_process:  # -- process gnomAD mtDNA
    input:
        dl="work/download/annos/{genome_build}/seqvars/gnomad_mtdna/{version}/gnomad.genomes.v{version}.sites.chrM.vcf.bgz",
    output:
        vcf="work/annos/{genome_build}/seqvars/gnomad_mtdna/{version}/gnomad_mtdna.vcf.gz",
        vcf_md5="work/annos/{genome_build}/seqvars/gnomad_mtdna/{version}/gnomad_mtdna.vcf.gz.md5",
        vcf_tbi="work/annos/{genome_build}/seqvars/gnomad_mtdna/{version}/gnomad_mtdna.vcf.gz.tbi",
        vcf_tbi_md5="work/annos/{genome_build}/seqvars/gnomad_mtdna/{version}/gnomad_mtdna.vcf.gz.tbi.md5",
    shell:
        r"""
        if [[ {wildcards.genome_build} == grch37 ]]; then
            zcat {input.dl} \
            | sed \
                -e 's/chrM/MT/g' \
                -e 's/GRCh38_MT/GRCh37/g'  \
                -e 's/,GERP_DIST/\&GERP_DIST/g' \
                -e 's/,BP_DIST/\&BP_DIST/g' \
                -e 's/,DIST_FROM_LAST_EXON/\&DIST_FROM_LAST_EXON/g' \
                -e 's/,50_BP_RULE/\&50_BP_RULE/g' \
                -e 's/,PHYLOCSF_TOO_SHORT/\&PHYLOCSF_TOO_SHORT/g' \
            | bgzip -c \
            > {output.vcf}
        else
            cp {input.dl} {output.vcf}
        fi

        tabix -f {output.vcf}

        md5sum {output.vcf} > {output.vcf_md5}
        md5sum {output.vcf_tbi} > {output.vcf_tbi_md5}
        """
