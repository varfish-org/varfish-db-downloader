## Download of hgnc_complete_set


rule annos_hgnc_complete_set_download:
    output:
        tsv="work/download/annos/hgnc/{quarterly_release_date}/hgnc_complete_set.tsv",
        json="work/download/annos/hgnc/{quarterly_release_date}/hgnc_complete_set.json",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.tsv} \
            https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_{wildcards.quarterly_release_date}.tsv

        wget --no-check-certificate \
            -O {output.json} \
            https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/json/hgnc_complete_set_{wildcards.quarterly_release_date}.json
        """
