rule mitomap_download:
    output:
        "work/download/pre-mehari/grch37/MITOMAP/{download_date}/polymorphisms.vcf",
    shell:
        r"""
        wget --no-check-certificate \
            --no-check-certificate \
            -O {output} \
            https://mitomap.org/cgi-bin/polymorphisms.cgi?format=vcf
        """


rule GRChXX_mitomap_normalize:
    input:
        vcf="work/download/pre-mehari/grch37/MITOMAP/{download_date}/polymorphisms.vcf",
        ref="work/reference/{genomebuild}/reference.fa",
    output:
        vcf="work/download/pre-mehari/{genomebuild}/MITOMAP/{download_date}/polymorphisms.normalized.annotated.vcf",
        tmp_vcf=temp("work/download/pre-mehari/{genomebuild}/MITOMAP/{download_date}/mt_tmp.vcf"),
        norm=temp("work/download/pre-mehari/{genomebuild}/MITOMAP/{download_date}/polymorphisms.normalized.vcf"),
        txt_tmp=temp("work/download/pre-mehari/{genomebuild}/MITOMAP/{download_date}/query-result.txt"),
        ann=temp("work/download/pre-mehari/{genomebuild}/MITOMAP/{download_date}/annotate.bed.gz"),
        anntbi=temp("work/download/pre-mehari/{genomebuild}/MITOMAP/{download_date}/annotate.bed.gz.tbi"),
    shell:
        r"""
        if [[ "{wildcards.genomebuild}" == "grch38" ]]; then
            sed 's/MT\b/chrM/' {input.vcf} > {output.tmp_vcf}
        else
            cp {input.vcf} {output.tmp_vcf}
        fi
        perl -p -e 's/;HGFL=[^;\s]*;?//g' {output.tmp_vcf} | perl -p -e 's/;FreqCR=[^;\s]*;?//g' \
        | bcftools norm \
            -m -any \
            -c w \
            -f {input.ref} \
            -o {output.norm} \
            -
        AN=$(bcftools view -h {output.norm} | sed -n 's/^[^"]\+"Allele count in GenBank out of \([[:digit:]]\+\) .*sequences">/\1/p')
        bcftools query \
            -f "%CHROM\t%POS\t%END\t%REF\t%ALT\t%AC\t$AN\n" \
            -o {output.txt_tmp} \
            {output.norm}
        uniq {output.txt_tmp} | awk -F$'\t' 'BEGIN{{OFS=FS}}{{$2-=1; print $0,$6/$7}}' | bgzip -c > {output.ann}
        tabix -p bed {output.ann}
        bcftools annotate \
            -a {output.ann} \
            -c CHROM,FROM,-,REF,ALT,-,INFO/AN,INFO/AF \
            -h <(echo -e "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number in Genbank\">\n##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">") \
            -o {output.vcf} \
            {output.norm}
        """


def input_mitomap_tsv(wildcards):
    return {
        "vcf": f"work/download/pre-mehari/{wildcards.genomebuild.lower()}/MITOMAP/{wildcards.download_date}/polymorphisms.normalized.annotated.vcf",
    }


rule result_GRChXX_mitomap_tsv:
    input:
        unpack(input_mitomap_tsv)
    output:
        tsv="output/pre-mehari/{genomebuild}/MITOMAP/{download_date}/Mitomap.tsv",
        release_info="output/pre-mehari/{genomebuild}/MITOMAP/{download_date}/Mitomap.release_info",
    shell:
        r"""
        (
            echo -e "release\tchromosome\tstart\tend\tbin\treference\talternative\tac\tan\taf"
            bcftools query \
                -f "{wildcards.genomebuild}\t%CHROM\t%POS\t%END\t\t%REF\t%ALT\t%AC\t%AN\t%AF\n" \
                {input.vcf}
        ) \
        | python rules/pre-mehari/tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nMitomap\t{wildcards.download_date}\t{wildcards.genomebuild}\t." > {output.release_info}
        """
