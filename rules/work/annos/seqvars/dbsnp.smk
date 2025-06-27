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
            reference={DV.dbsnp_reference_37}
        else
            reference={DV.dbsnp_reference_38}
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
            reference={DV.dbsnp_reference_37_ext}
            reference_report={DV.dbsnp_reference_37_report}
        else
            reference={DV.dbsnp_reference_38_ext}
            reference_report={DV.dbsnp_reference_38_report}
        fi

        wget --no-check-certificate \
            -O {output.txt} \
            https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/$reference_report/${{reference}}_assembly_report.txt
        """


# TODO
# Note: annonars uses an archived version, b151
# VarFish had a rolling release of dbSNP, with the latest version being b155.
#
# Download links used for annonars for b151 are:
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
# However, the files in this folders never received any update.
#
# Download links for the latest version used for VarFish are:
# https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
# https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz
#
# Best would be to use a tagged version:
# https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.25.gz
# https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz
