## Rules related to the AlphaMissense scores


rule annos_alphamissense_download:  # -- download AlphaMissense data
    output:
        tsv="work/download/annos/alphamissense/1/{genome}/AlphaMissense_{genome}.tsv.gz",
        tsv_md5="work/download/annos/alphamissense/1/{genome}/AlphaMissense_{genome}.tsv.gz.md5",
    shell:
        r"""
        wget --no-check-certificate \
            -O {output.tsv} \
            https://storage.googleapis.com/dm_alphamissense/AlphaMissense_{wildcards.genome}.tsv.gz

        md5sum {output.tsv} >{output.tsv_md5}
        """
