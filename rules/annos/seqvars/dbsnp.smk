## Rules related to dbSNP.


rule annos_dbsnp_download_grch37:  # -- download dbSNP data for GRCh38
    output:
        vcf="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.vcf.gz",
        vcf_tbi="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.vcf.gz.tbi",
    shell:
        r"""
        # Check the version.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        wget -O $TMPDIR/listing \
            https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF
        version=$(grep -w '00-All.vcf.gz"' $TMPDIR/listing | awk '{{ print $3 }}')
        if [[ "$version" != "{wildcards.version}" ]]; then
            >&2 echo "Version mismatch for 00-All.vcf.gz: expected {wildcards.version}, got $version"
            exit 1
        fi

        # Perform the actual download.
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
        tabix -f {output.vcf}
        """
