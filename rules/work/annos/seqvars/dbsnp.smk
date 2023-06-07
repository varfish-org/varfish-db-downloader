## Rules related to dbSNP.


rule annos_dbsnp_download:  # -- download dbSNP data
    output:
        vcf="work/download/annos/{genome_release}/seqvars/dbsnp/{version}/dbsnp.vcf.gz",
        vcf_tbi="work/download/annos/{genome_release}/seqvars/dbsnp/{version}/dbsnp.vcf.gz.tbi",
    shell:
        r"""
        # Check the version.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        if [[ "{wildcards.genome_release}" == grch37 ]]; then
            reference=human_9606_b151_GRCh37p13
        else
            reference=human_9606_b151_GRCh38p7
        fi

        # Perform the actual download.
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            https://ftp.ncbi.nih.gov/snp/organisms/$reference/VCF/00-All.vcf.gz
        tabix -f {output.vcf}
        """
