# Download dbSNP, map chr names, normalize, and convert into TSV for import.


rule grch37_dbsnp_b155_download:
    output:
        vcf="GRCh37/dbSNP/b155/download/GCF_000001405.25.gz",
        tbi="GRCh37/dbSNP/b155/download/GCF_000001405.25.gz.tbi",
    log:
        "GRCh37/dbSNP/b155/download/GCF_000001405.25.gz.log",
    shell:
        r"""
        wget --no-check-certificate \
            -o {log} \
            -O {output.vcf} \
            http://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
        wget --no-check-certificate \
            -a {log} \
            -O {output.tbi} \
            http://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch37_dbsnp_b155_map_chr:
    input:
        vcf="GRCh37/dbSNP/b155/download/GCF_000001405.25.gz",
        tbi="GRCh37/dbSNP/b155/download/GCF_000001405.25.gz.tbi",
    output:
        map=temp("GRCh37/dbSNP/b155/download/GCF_000001405.25.map_chr"),
        vcf="GRCh37/dbSNP/b155/download/GCF_000001405.25.map_chr.gz",
        tbi="GRCh37/dbSNP/b155/download/GCF_000001405.25.map_chr.gz.tbi",
    shell:
        r"""
        awk -v RS="(\r)?\n" 'BEGIN {{ FS="\t" }} !/^#/ {{ if ($10 != "na") print $7,$10; else print $7,$5 }}' \
            tools/data/GCF_000001405.25_GRCh37.p13_assembly_report.txt \
        | sed -e 's/chrM/chrMT/g' \
        | sed -e 's/chr//g' \
        > {output.map}

        bcftools annotate --threads=2 --rename-chrs {output.map} {input.vcf} -O z -o {output.vcf}
        tabix -f {output.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch37_dbsnp_b155_normalize:
    input:
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        vcf="GRCh37/dbSNP/b155/download/GCF_000001405.25.map_chr.gz",
    output:
        vcf="GRCh37/dbSNP/b155/download/GCF_000001405.25.normalized.{chrom}.vcf.gz",
        tbi="GRCh37/dbSNP/b155/download/GCF_000001405.25.normalized.{chrom}.vcf.gz.tbi",
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


rule result_grch37_dbsnp_b155_tsv:
    input:
        header="header/dbsnp.txt",
        vcf="GRCh37/dbSNP/b155/download/GCF_000001405.25.normalized.{chrom}.vcf.gz",
        tbi="GRCh37/dbSNP/b155/download/GCF_000001405.25.normalized.{chrom}.vcf.gz.tbi",
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
            | sort -u -t $'\t' -k 2,2 -k 3,3 -k 6,6 -k 7,7 -S {config[sort_memory]}
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nDbsnp\tb155\tGRCh37\t" > {output.release_info}
        """
