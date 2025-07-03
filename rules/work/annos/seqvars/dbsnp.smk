## Rules related to dbSNP.


rule annos_dbsnp_download:  # -- download dbSNP data
    output:
        vcf="work/download/annos/{genome_release}/seqvars/dbsnp/{version}/dbsnp.vcf.gz",
        vcf_tbi="work/download/annos/{genome_release}/seqvars/dbsnp/{version}/dbsnp.vcf.gz.tbi",
    params:
        reference=lambda wildcards: DV.refseq_ref_38 if wildcards.genome_release == "grch37" else DV.refseq_ref_37
    shell:
        r"""
        # Check the version.
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT

        # Perform the actual download.
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            https://ftp.ncbi.nih.gov/snp/archive/{wildcards.version}/VCF/{params.reference}.gz
        tabix -f {output.vcf}
        """


rule annos_dbsnp_assembly_release:
    output:
        txt="work/download/annos/{genome_release}/seqvars/dbsnp/{version}/assembly_report.txt",
    params:
        assembly=lambda wildcards: (
            DV.refseq_ref_37_assembly if wildcards.genome_release == "grch37"
            else DV.refseq_ref_38_assembly
        ),
        version=lambda wildcards: DV.refseq_37 if wildcards.genome_release == "grch37" else DV.refseq_38
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.txt} \
            {DV.refseq_base_url}/{params.version}/{params.assembly}/{params.assembly}_assembly_report.txt
        """

