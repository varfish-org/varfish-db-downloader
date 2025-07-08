# Obtain "knownGeneAA" track from ENSEMBL and convert into TSV via VCF.

# rule grchxx_knowngeneaa_download:
#     output:
#         fa="{genome_build}/knowngeneaa/{download_date}/download/knownGene.exonAA.fa.gz",
#     shell:
#         r"""
#         if [[ {wildcards.genome_build} == GRCh37 ]]; then
#             ucsc_name=hg19
#         else
#             ucsc_name=hg38
#         fi

#         wget --no-check-certificate \
#             -O {output.fa} \
#             "http://hgdownload.cse.ucsc.edu/goldenpath/${{ucsc_name}}/multiz100way/alignments/knownGene.exonAA.fa.gz"
#         """


def input_grchxx_knowngeneaa_to_vcf(wildcards):
    genome_build = wildcards.genome_build.lower()
    return {
        "reference": f"work/reference/{genome_build}/reference.fa",
        "fa": f"work/download/annos/{genome_build}/features/cons/{wildcards.version}/knownGene.exonAA.fa.gz"
    }


rule grchxx_knowngeneaa_to_vcf:
    input:
        unpack(input_grchxx_knowngeneaa_to_vcf),
    output:
        vcf="work/download/pre-mehari/{genome_build}/knowngeneaa/{version}/knownGeneAA.vcf.gz",
        tbi="work/download/pre-mehari/{genome_build}/knowngeneaa/{version}/knownGeneAA.vcf.gz.tbi",
    shell:
        r"""
        python rules/pre-mehari/tools/knowngeneaa.py \
            {input.reference} \
            {input.fa} \
            --output /dev/stdout \
        | bcftools sort \
            -O z \
            -o {output.vcf}
        tabix -f {output.vcf}
        """


# def input_result_grchxx_knowngeneaa_to_tsv(wildcards):
#     """Input function for ``rule result_grchxx_knowngeneaa_to_tsv``."""
#     return {
#         "header": "rules/pre-mehari/header/knowngeneaa.txt",
#         "vcf": (
#             f"work/download/annos/{wildcards.genome_build.lower()}/features/cons/"
#             f"{wildcards.version}/ucsc_conservation.vcf.gz"
#         ),
#         "tbi": (
#             f"work/download/annos/{wildcards.genome_build.lower()}/features/cons/"
#             f"{wildcards.version}/ucsc_conservation.vcf.gz.tbi"
#         ),
#     }


rule result_grchxx_knowngeneaa_to_tsv:
    input:
        header="rules/pre-mehari/header/knowngeneaa.txt",
        vcf="work/download/pre-mehari/{genome_build}/knowngeneaa/{version}/knownGeneAA.vcf.gz",
        tbi="work/download/pre-mehari/{genome_build}/knowngeneaa/{version}/knownGeneAA.vcf.gz.tbi",
    output:
        tsv="output/pre-mehari/{genome_build}/knowngeneaa/{version}/KnowngeneAA.tsv",
        release_info="output/pre-mehari/{genome_build}/knowngeneaa/{version}/KnowngeneAA.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query \
                -f "{wildcards.genome_build}\t%CHROM\t%POS\t%END\t\t%UCSC_GENE\t\t%ALIGNMENT\n" \
                {input.vcf} \
            | awk -F $'\t' '
                BEGIN {{ OFS = FS }}
                {{
                    $7 = substr($6, 1, length($6) - 2)
                    print
                }}'
        ) \
        | python rules/pre-mehari/tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nKnowngeneAA\t{wildcards.version}\t{wildcards.genome_build}\t" > {output.release_info}
        """
