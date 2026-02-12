## Rules related to the CADD score and its feature annotations.

#: Prefix to the CADD download URL.
CADD_PREFIX = f"https://kircherlab.bihealth.org/download/CADD"


rule annos_seqvars_cadd_download:  # -- download CADD data
    output:
        tsv="work/download/annos/{genome_release}/seqvars/cadd/{version}/{filename}.tsv.gz",
        tsv_tbi="work/download/annos/{genome_release}/seqvars/cadd/{version}/{filename}.tsv.gz.tbi",
    shell:
        r"""
        # Remove .bihealth.org from no_proxy variables to allow direct access
        if [ -n "${{no_proxy:-}}" ]; then
            export no_proxy=$(echo "$no_proxy" | sed -E 's/(^|,|;)\.?bihealth\.org(,|;|$)/\1/g; s/(^|,|;)(,|;)/\1/g; s/^,//; s/,$//; s/^;//; s/;$//')
        fi
        if [ -n "${{NO_PROXY:-}}" ]; then
            export NO_PROXY=$(echo "$NO_PROXY" | sed -E 's/(^|,|;)\.?bihealth\.org(,|;|$)/\1/g; s/(^|,|;)(,|;)/\1/g; s/^,//; s/,$//; s/^;//; s/;$//')
        fi

        for path in {output};
        do
            aria2c \
                --check-certificate=false \
                --file-allocation=trunc \
                --out=$path \
                --split=16 \
                --max-concurrent-downloads=16 \
                --max-connection-per-server=16 \
                {CADD_PREFIX}/v{wildcards.version}/$(echo {wildcards.genome_release} | sed -e 's/grch/GRCh/')/$(basename $path)
        done
        """


# rule annos_cadd_process_37:  # -- process CADD data for GRCh37
#     input:
#         "work/download/annos/grch37/seqvars/cadd/v{version}/whole_genome_SNVs_inclAnno.tsv.gz",
#         "work/download/annos/grch37/seqvars/cadd/v{version}/whole_genome_SNVs_inclAnno.tsv.gz.tbi",
#         "work/download/annos/grch37/seqvars/cadd/v{version}/InDels_inclAnno.tsv.gz",
#         "work/download/annos/grch37/seqvars/cadd/v{version}/InDels_inclAnno.tsv.gz.tbi",
#     output:
#         touch("annos/grch37/cadd/.done"),
# rule annos_cadd_process_38:  # -- process CADD data for GRCh38
#     input:
#         "work/download/annos/grch38/seqvars/cadd/v{version}/whole_genome_SNVs_inclAnno.tsv.gz",
#         "work/download/annos/grch38/seqvars/cadd/v{version}/whole_genome_SNVs_inclAnno.tsv.gz.tbi",
#         "work/download/annos/grch38/seqvars/cadd/v{version}/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz",
#         "work/download/annos/grch38/seqvars/cadd/v{version}/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz.tbi",
#     output:
#         touch("annos/grch38/cadd/.done"),
