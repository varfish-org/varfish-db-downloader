## Rules related to HelixMtDb.


rule annos_seqvars_helixmtdb_download:  # -- download HelixMtDb data
    output:
        tsv="work/download/annos/{genome_build}/seqvars/helixmtdb/20200327/helixmtdb.tsv",
    shell:
        r"""
        wget \
            --no-check-certificate \
            -O {output} \
            https://helix-research-public.s3.amazonaws.com/mito/HelixMTdb_20200327.tsv
        """


rule annos_seqvars_helixmtdb_convert:  # -- process HelixMtDb data
    input:
        tsv="work/download/annos/{genome_build}/seqvars/helixmtdb/20200327/helixmtdb.tsv",
    output:
        vcf="work/annos/{genome_build}/seqvars/helixmtdb/20200327/helixmtdb.vcf.gz",
        vcf_tbi="work/annos/{genome_build}/seqvars/helixmtdb/20200327/helixmtdb.vcf.gz.tbi",
    shell:
        r"""
        cat {input.tsv} \
        | python3  scripts/helix-to-vcf.py \
        > {output.vcf}.tmp

        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            sed -e 's/chrM/MT/g' {output.vcf}.tmp \
            | bgzip -c \
            > {output.vcf}
        else
            bgzip -c {output.vcf}.tmp >{output.vcf}
        fi

        tabix -f {output.vcf}

        rm -f {output.vcf}.tmp
        """
