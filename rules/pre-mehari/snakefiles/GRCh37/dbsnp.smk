# Download dbSNP, map chr names, normalize, and convert into TSV for import.


rule grch37_dbsnp_map_chr:
    input:
        vcf="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.vcf.gz",
        tbi="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.vcf.gz.tbi",
        report="work/download/annos/grch37/seqvars/dbsnp/{version}/assembly_report.txt",
    output:
        map=temp("work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.map_chr"),
        vcf="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.map_chr.gz",
        tbi="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.map_chr.gz.tbi",
    shell:
        r"""
        awk -v RS="(\r)?\n" 'BEGIN {{ FS="\t" }} !/^#/ {{ if ($10 != "na") print $7,$10; else print $7,$5 }}' \
            {input.report} \
        | sed -e 's/chrM/chrMT/g' \
        | sed -e 's/chr//g' \
        > {output.map}

        bcftools annotate --threads=2 --rename-chrs {output.map} {input.vcf} -O z -o {output.vcf}
        tabix -f {output.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch37_dbsnp_normalize:
    input:
        reference="work/reference/grch37/reference.fa",
        vcf="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.map_chr.gz",
    output:
        vcf="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.normalized.{chrom}.vcf.gz",
        tbi="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.normalized.{chrom}.vcf.gz.tbi",
    shell:
        r"""
        bcftools norm \
            --check-ref s \
            --regions "{wildcards.chrom}" \
            --threads 16 \
            --multiallelics - \
            --fasta-ref {input.reference} \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        tabix -f {output.vcf}
        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        """


rule result_grch37_dbsnp_tsv:
    input:
        header="rules/pre-mehari/header/dbsnp.txt",
        vcf="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.normalized.{chrom}.vcf.gz",
        tbi="work/download/annos/grch37/seqvars/dbsnp/{version}/dbsnp.normalized.{chrom}.vcf.gz.tbi",
    output:
        release_info="GRCh37/dbSNP/b155/Dbsnp.{chrom}.release_info",
        tsv="GRCh37/dbSNP/b155/Dbsnp.{chrom}.tsv",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query {input.vcf} \
                -f 'GRCh37\t%CHROM\t%POS\t%END\t\t%REF\t%ALT\t%ID\n' \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} ((length($6) <= 512) && (length($7) <= 512)) {{ print }}' \
            | sort -u -t $'\t' -k 2,2 -k 3,3 -k 6,6 -k 7,7
        ) \
        | python rules/pre-mehari/tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nDbsnp\t{wildcards.version}\tGRCh37\t" > {output.release_info}
        """
