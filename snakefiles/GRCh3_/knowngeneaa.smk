# Obtain "knownGeneAA" track from ENSEMBL and convert into TSV via VCF.


rule grchxx_knowngeneaa_download:
    output:
        fa="{genome_build}/knowngeneaa/{download_date}/download/knownGene.exonAA.fa.gz",
    shell:
        r"""
        if [[ {wildcards.genome_build} == GRCh37 ]]; then
            ucsc_name=hg19
        else
            ucsc_name=hg38
        fi

        wget \
            -O {output.fa} \
            "http://hgdownload.cse.ucsc.edu/goldenpath/${ucsc_name}/multiz100way/alignments/knownGene.exonAA.fa.gz"
        """


def input_grchxx_knowngeneaa_to_vcf(wildcards):
    if wildcards.genome_build == "GRCh37":
        return {
            "reference": "GRCh37/reference/hs37d5/hs37d5.fa",
            "fa": "{genome_build}/knowngeneaa/{download_date}/download/knownGene.exonAA.fa.gz".format(
                **wildcards
            ),
        }
    else:
        return {
            "reference": "GRCh38/reference/hs38/hs38.fa",
            "fa": "{genome_build}/knowngeneaa/{download_date}/download/knownGene.exonAA.fa.gz".format(
                **wildcards
            ),
        }


rule grchxx_knowngeneaa_to_vcf:
    input:
        unpack(input_grchxx_knowngeneaa_to_vcf),
    output:
        vcf="{genome_build}/knowngeneaa/{download_date}/knownGeneAA.vcf.gz",
        tbi="{genome_build}/knowngeneaa/{download_date}/knownGeneAA.vcf.gz.tbi",
    shell:
        r"""
        python tools/knowngeneaa.py \
            {input.reference} \
            {input.fa} \
            --output /dev/stdout \
        | bcftools sort \
            -O z \
            -o {output.vcf}
        tabix -f {output.vcf}
        """


rule result_grchxx_knowngeneaa_to_tsv:
    input:
        header="header/knowngeneaa.txt",
        vcf="{genome_build}/knowngeneaa/{download_date}/knownGeneAA.vcf.gz",
    output:
        tsv="{genome_build}/knowngeneaa/{download_date}/KnowngeneAA.tsv",
        release_info="{genome_build}/knowngeneaa/{download_date}/KnowngeneAA.release_info",
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
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nKnowngeneAA\t$(date +%Y/%m/%d)\t{wildcards.genome_build}\t" > {output.release_info}
        """
