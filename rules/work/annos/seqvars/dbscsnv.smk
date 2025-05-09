## Rules related to dbscSNV.


def files_dbscsnv(version: str = "1.1"):
    """Files contained in the dbscSNV ZIP file."""
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    return [
        f"work/download/annos/grch37/seqvars/dbscsnv/{version}/dbscSNV{version}.chr{chrom}"
        for chrom in chroms
    ]


def var_tsv_files_dbscsnv(version: str = "1.1"):
    """Variant TSV files contained in the dbscSNV ZIP file."""
    return [entry for entry in files_dbscsnv(version) if ".chr" in entry]


rule annos_seqvars_dbscsnv_download:  # -- download dbscSNV ZIP file
    output:
        zip="work/download/annos/grch37/seqvars/dbscsnv/1.1/dbscSNV1.1.zip",
    shell:
        r"""
        aria2c \
            --check-certificate=false \
            --file-allocation=trunc \
            --out={output.zip} \
            --split=8 \
            --max-concurrent-downloads=8 \
            --max-connection-per-server=8 \
            https://usf.box.com/shared/static/ffwlywsat3q5ijypvunno3rg6steqfs8
        """


rule annos_seqvars_dbscsnv_extract:  # -- extract dbscSNV ZIP file
    input:
        zip="work/download/annos/grch37/seqvars/dbscsnv/1.1/dbscSNV1.1.zip",
    output:
        files_dbscsnv(),
    shell:
        r"""
        unzip -d $(dirname {input.zip}) {input.zip}
        """
