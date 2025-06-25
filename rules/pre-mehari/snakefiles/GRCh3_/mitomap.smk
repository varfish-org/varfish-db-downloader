rule GRCh37_mitomap_download:
    output:
        "GRCh37/MITOMAP/{download_date}/download/polymorphisms.vcf",
    shell:
        r"""
        wget --no-check-certificate \
            --no-check-certificate \
            -O {output} \
            https://mitomap.org/cgi-bin/polymorphisms.cgi?format=vcf
        """


rule GRChXX_mitomap_normalize:
    input:
        vcf="GRCh37/MITOMAP/{download_date}/download/polymorphisms.vcf",
        ref="GRCh37/reference/hs37d5/hs37d5.fa",
    output:
        vcf="{reference}/MITOMAP/{download_date}/download/polymorphisms.normalized.annotated.vcf",
        norm=temp("{reference}/MITOMAP/{download_date}/download/polymorphisms.normalized.vcf"),
        txt_tmp=temp("{reference}/MITOMAP/{download_date}/download/query-result.txt"),
        ann=temp("{reference}/MITOMAP/{download_date}/download/annotate.bed.gz"),
        anntbi=temp("{reference}/MITOMAP/{download_date}/download/annotate.bed.gz.tbi"),
    shell:
        r"""
        perl -p -e 's/;HGFL=[^;\s]*;?//g' {input.vcf} | perl -p -e 's/;FreqCR=[^;\s]*;?//g' \
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
        ls -lh {output.ann}
        tabix -p bed {output.ann}
        bcftools annotate \
            -a {output.ann} \
            -c CHROM,FROM,-,REF,ALT,-,INFO/AN,INFO/AF \
            -h <(echo -e "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number in Genbank\">\n##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">") \
            -o {output.vcf} \
            {output.norm}
        """


rule result_GRChXX_mitomap_tsv:
    input:
        vcf=(
            "{genome_build}/MITOMAP/{download_date}/download/polymorphisms.normalized.annotated.vcf"
        ),
        header="header/mitomap.txt",
    output:
        tsv="{genome_build}/MITOMAP/{download_date}/Mitomap.tsv",
        release_info="{genome_build}/MITOMAP/{download_date}/Mitomap.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query \
                -f "{wildcards.genome_build}\t%CHROM\t%POS\t%END\t\t%REF\t%ALT\t%AC\t%AN\t%AF\n" \
                {input.vcf}
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nMitomap\t$(date +%Y/%m/%d)\t{wildcards.genome_build}\t." > {output.release_info}
        """
