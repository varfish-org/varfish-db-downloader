## Rules related to dbNSFP.


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
    return ["work/download/annos/grch37/seqvars/dbnsfp/{version}{variant}/%s" % e for e in lst]


def var_tsv_dbnsfp():
    """Return paths to variant TSV files."""
    return [path for path in files_dbnsfp() if "_variant" in path and path.endswith(".gz")]


rule annos_seqvars_dbnsfp_download:  # -- download dbNSFP ZIP file
    output:
        zip="work/download/annos/grch37/seqvars/dbnsfp/{version}{variant}/dbNSFP{version}{variant}.zip",
    wildcard_constraints:
        version=r"\d\.\d",
    threads: 8
    shell:
        r"""
        url_a=https://usf.box.com/shared/static/2hzcx5s6p1xui7oen16xqzndfrkt8l9l
        url_c=https://usf.box.com/shared/static/03xsrpna0nzgrytfo2pzk326t8jad4oc

        case "{wildcards.variant}" in
          a) url=$url_a ;;
          c) url=$url_c ;;
          *) echo "Unsupported dbNSFP variant '{wildcards.variant}'" >&2; exit 1 ;;
        esac

        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.zip} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            $url
        """


rule annos_seqvars_dbnsfp_extract:  # -- extract dbNSFP ZIP file
    input:
        zip="work/download/annos/grch37/seqvars/dbnsfp/{version}{variant}/dbNSFP{version}{variant}.zip",
    output:
        files_dbnsfp(),
    wildcard_constraints:
        version=r"\d\.\d",
    shell:
        r"""
        unzip -d $(dirname {input.zip}) {input.zip}
        """
