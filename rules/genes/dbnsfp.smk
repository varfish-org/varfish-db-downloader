## Rules related to dbNSFP gene information.

#: Download URL for dbNSFP 4.4a
DBNSFP_ACADEMIC_URL = "https://usf.box.com/shared/static/bvfzmkpgtphvbmmrvb2iyl2jl21o49kc"
#: Download URL for dbNSFP 4.4c
DBNSFP_COMMMERCIAL_URL = "https://usf.box.com/shared/static/a84zcdlkx2asq2nxh6xr2gdb4csmyvhk"


def files_dbnsfp():
    """Helper that returns the files within the dbNSFP archive."""
    lst = [
        "dbNSFP{version}{variant}.readme.txt",
        "dbNSFP{version}{variant}_variant.chr1.gz",
        "dbNSFP{version}{variant}_variant.chr10.gz",
        "dbNSFP{version}{variant}_variant.chr11.gz",
        "dbNSFP{version}{variant}_variant.chr12.gz",
        "dbNSFP{version}{variant}_variant.chr13.gz",
        "dbNSFP{version}{variant}_variant.chr14.gz",
        "dbNSFP{version}{variant}_variant.chr15.gz",
        "dbNSFP{version}{variant}_variant.chr16.gz",
        "dbNSFP{version}{variant}_variant.chr17.gz",
        "dbNSFP{version}{variant}_variant.chr18.gz",
        "dbNSFP{version}{variant}_variant.chr19.gz",
        "dbNSFP{version}{variant}_variant.chr2.gz",
        "dbNSFP{version}{variant}_variant.chr20.gz",
        "dbNSFP{version}{variant}_variant.chr21.gz",
        "dbNSFP{version}{variant}_variant.chr22.gz",
        "dbNSFP{version}{variant}_variant.chr3.gz",
        "dbNSFP{version}{variant}_variant.chr4.gz",
        "dbNSFP{version}{variant}_variant.chr5.gz",
        "dbNSFP{version}{variant}_variant.chr6.gz",
        "dbNSFP{version}{variant}_variant.chr7.gz",
        "dbNSFP{version}{variant}_variant.chr8.gz",
        "dbNSFP{version}{variant}_variant.chr9.gz",
        "dbNSFP{version}{variant}_variant.chrM.gz",
        "dbNSFP{version}{variant}_variant.chrX.gz",
        "dbNSFP{version}{variant}_variant.chrY.gz",
        "dbNSFP{version}_gene.complete.gz",
        "dbNSFP{version}_gene.gz",
        "LICENSE.txt",
        "try.vcf",
        "tryhg18.in",
        "tryhg19.in",
        "tryhg38.in",
    ]
    return ["work/download/grch37/seqvars/dbnsfp-{version}{variant}/%s" % e for e in lst]


rule genes_dbnsfp_download:  # -- download dbNSFP ZIP file
    output:
        zip="work/download/grch37/seqvars/dbnsfp-{version}{variant}/dbNSFP{version}{variant}.zip",
    shell:
        r"""
        if [[ "{wildcards.variant}" == a ]]; then
            url={DBNSFP_ACADEMIC_URL}
        else
            url={DBNSFP_COMMMERCIAL_URL}
        fi

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.zip} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            $url
        """


rule genes_dbnsfp_extract:  # -- extract dbNSFP ZIP file
    input:
        zip="work/download/grch37/seqvars/dbnsfp-{version}{variant}/dbNSFP{version}{variant}.zip",
    output:
        files_dbnsfp(),
    shell:
        r"""
        unzip -d $(dirname {input.zip}) {input.zip}
        """


rule genes_dbnsfp_genes_copy:  # -- copy over dbNSFP genes file
    input:
        tsv=f"work/download/grch37/seqvars/dbnsfp-{DV.dbnsfp}a/dbNSFP{DV.dbnsfp}_gene.complete.gz",
    output:
        tsv="work/genes/dbnsfp/genes.tsv.gz",
        tsv_md5="work/genes/dbnsfp/genes.tsv.gz.md5",
    shell:
        r"""
        zcat {input.tsv} \
        | pigz -c \
        > {output.tsv}

        md5sum {output.tsv} >{output.tsv_md5}
        """
