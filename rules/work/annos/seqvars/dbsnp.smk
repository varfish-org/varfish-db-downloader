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
            reference={DV.refseq_ref_37}
        else
            reference={DV.refseq_ref_38}
        fi

        # Perform the actual download.
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.vcf} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            https://ftp.ncbi.nih.gov/snp/archive/{wildcards.version}/VCF/$reference.gz
        tabix -f {output.vcf}
        """


rule annos_dbsnp_assembly_release:
    output:
        txt="work/download/annos/{genome_release}/seqvars/dbsnp/{version}/assembly_report.txt",
    shell:
        r"""
        if [[ "{wildcards.genome_release}" == grch37 ]]; then
            reference={DV.refseq_ref_37_assembly}
            reference_report={DV.refseq_ref_37_date}
        else
            reference={DV.refseq_ref_38_assembly}
            reference_report={DV.refseq_ref_38_date}
        fi

        wget --no-check-certificate \
            -O {output.txt} \
            {DV.refseq_base_url}/$reference_report/${{reference}}_assembly_report.txt
        """

