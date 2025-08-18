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

        bcftools annotate --threads 16 --rename-chrs {output.map} {input.vcf} -O z -o {output.vcf}
        tabix -f {output.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch38_dbsnp_map_chr:
    input:
        vcf="work/download/annos/grch38/seqvars/dbsnp/{version}/dbsnp.vcf.gz",
        tbi="work/download/annos/grch38/seqvars/dbsnp/{version}/dbsnp.vcf.gz.tbi",
        report="work/download/annos/grch38/seqvars/dbsnp/{version}/assembly_report.txt",
    output:
        map=temp("work/download/annos/grch38/seqvars/dbsnp/{version}/dbsnp.map_chr"),
        vcf="work/download/annos/grch38/seqvars/dbsnp/{version}/dbsnp.map_chr.gz",
        tbi="work/download/annos/grch38/seqvars/dbsnp/{version}/dbsnp.map_chr.gz.tbi",
    threads: THREADS
    shell:
        r"""
        awk -v RS="(\r)?\n" 'BEGIN {{ FS="\t" }} !/^#/ {{ if ($10 != "na") print $7,$10; else print $7,$5 }}' \
            {input.report} \
        > {output.map}

        bcftools annotate --threads {threads} --rename-chrs {output.map} {input.vcf} -O z -o {output.vcf}
        tabix -f {output.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch3x_dbsnp_normalize:
    input:
        reference="work/reference/{genomebuild}/reference.fa",
        vcf="work/download/annos/{genomebuild}/seqvars/dbsnp/{version}/dbsnp.map_chr.gz",
    output:
        vcf="work/download/annos/{genomebuild}/seqvars/dbsnp/{version}/dbsnp.normalized.{chrom}.vcf.gz",
        tbi="work/download/annos/{genomebuild}/seqvars/dbsnp/{version}/dbsnp.normalized.{chrom}.vcf.gz.tbi",
    params:
        chrom=lambda wildcards: f"chr{wildcards.chrom}"
        if wildcards.genomebuild == "grch38"
        else wildcards.chrom,
    threads: THREADS
    shell:
        r"""
        bcftools norm \
            --check-ref s \
            --regions "{params.chrom}" \
            --threads {threads} \
            --multiallelics - \
            --fasta-ref {input.reference} \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        tabix -f {output.vcf}
        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        """


def input_result_dbsnp_tsv(wildcards):
    return {
        "header": "rules/pre-mehari/header/dbsnp.txt",
        "vcf": f"work/download/annos/{wildcards.genomebuild.lower()}/seqvars/dbsnp/{wildcards.version}/dbsnp.normalized.{wildcards.chrom}.vcf.gz",
        "tbi": f"work/download/annos/{wildcards.genomebuild.lower()}/seqvars/dbsnp/{wildcards.version}/dbsnp.normalized.{wildcards.chrom}.vcf.gz.tbi",
    }


rule result_grch3x_dbsnp_tsv:
    input:
        unpack(input_result_dbsnp_tsv),
    output:
        release_info="output/pre-mehari/{genomebuild}/dbSNP/{version}/Dbsnp.{chrom}.release_info",
        tsv="output/pre-mehari/{genomebuild}/dbSNP/{version}/Dbsnp.{chrom}.tsv",
    threads: THREADS
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query {input.vcf} \
                -f '{wildcards.genomebuild}\t%CHROM\t%POS\t%END\t\t%REF\t%ALT\t%ID\n' \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} ((length($6) <= 512) && (length($7) <= 512)) {{ print }}' \
            | sort -S 1G --parallel={threads} -u -t $'\t' -k 2,2 -k 3,3 -k 6,6 -k 7,7
        ) \
        | python rules/pre-mehari/tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nDbsnp\t{wildcards.version}\t{wildcards.genomebuild}\t" > {output.release_info}
        """
