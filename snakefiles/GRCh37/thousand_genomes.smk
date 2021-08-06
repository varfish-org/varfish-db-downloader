# Obtain thousand genomes phase 3 calls and import.


rule grch37_thousand_genomes_download:
    output:
        file="GRCh37/thousand_genomes/phase3/download/{filename}",
        md5="GRCh37/thousand_genomes/phase3/download/{filename}.md5",
    log:
        "GRCh37/thousand_genomes/phase3/download/{filename}.log",
    shell:
        r"""
        wget \
            -O {output.file} \
            -o {log} \
            http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/{wildcards.filename}

        pushd $(dirname {output.file})
        md5sum {wildcards.filename} > {wildcards.filename}.md5
        """


rule grch37_thousand_genomes_normalize:
    input:
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        vcf="GRCh37/thousand_genomes/phase3/download/{filename}",
        tbi="GRCh37/thousand_genomes/phase3/download/{filename}.tbi",
    output:
        vcf="GRCh37/thousand_genomes/phase3/download/normalized.{filename}",
    log:
        "GRCh37/thousand_genomes/phase3/download/{filename}.log",
    shell:
        r"""
        bcftools norm \
            --threads 8 \
            --multiallelics - \
            --fasta-ref {input.reference} \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        """


rule grch37_thousand_genomes_sv_download:
    output:
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz",
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.tbi",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi",
        "GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20130502.ALL.panel",
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.md5",
        "GRCh37/thousand_genomes/phase3/download/ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz.tbi.md5",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.md5",
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi.md5",
        "GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20130502.ALL.panel.md5",
    log:
        "GRCh37/thousand_genomes/phase3/download/wget.log",
    shell:
        r"""
        cd $(dirname {output[0]})

        echo "mget ALL.* README*" \
            | lftp http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel \
            -o $(basename {log})

        for i in {output[0]} {output[1]} {output[2]} {output[3]} {output[4]}
        do
            md5sum $(basename $i) > $(basename $i).md5
        done
        """


rule result_grch37_thousand_genomes_sv_tsv:
    input:
        # ``tbi`` is required for vcf parsing but not passed to the to_tsv function.
        "GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi",
        # ``panel`` is required but not passed to the to_tsv function.
        # It is hard-coded in the function and passed to the Converter constructor.
        "GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20130502.ALL.panel",
        header="header/thousand_genomes_sv.txt",
        vcf="GRCh37/thousand_genomes/phase3/download/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
    output:
        tsv="GRCh37/thousand_genomes/phase3/ThousandGenomesSv.tsv",
        release_info="GRCh37/thousand_genomes/phase3/ThousandGenomesSv.release_info",
    run:
        to_tsv(input.vcf, output.tsv, output.release_info, input.header)


ruleorder: grch37_thousand_genomes_normalize > grch37_thousand_genomes_download


rule grch37_thousand_genomes_gt_to_site:
    input:
        panel=(
            "GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20130502.ALL.panel"
        ),
        ped="GRCh37/thousand_genomes/phase3/download/integrated_call_samples_v3.20200731.ALL.ped",
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        vcf="GRCh37/thousand_genomes/phase3/download/normalized.{filename}.genotypes.vcf.gz",
    output:
        vcf="GRCh37/thousand_genomes/phase3/{filename}.sites.vcf.gz",
        vcf_md5="GRCh37/thousand_genomes/phase3/{filename}.sites.vcf.gz.md5",
        tbi="GRCh37/thousand_genomes/phase3/{filename}.sites.vcf.gz.tbi",
        tbi_md5="GRCh37/thousand_genomes/phase3/{filename}.sites.vcf.gz.tbi.md5",
    shell:
        r"""
        var-agg \
            --io-threads 16 \
            --input-panel {input.panel} \
            --input-ped {input.ped} \
            --input-fasta {input.reference} \
            --output {output.vcf} \
            {input.vcf}

        tabix -f {output.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch37_thousand_genomes_joint:
    input:
        vcf=expand(
            "GRCh37/thousand_genomes/phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
            chr=range(1, 23),
        )
        + [
            "GRCh37/thousand_genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.sites.vcf.gz",
            "GRCh37/thousand_genomes/phase3/ALL.chrY.phase3_integrated_v2b.20130502.sites.vcf.gz",
            "GRCh37/thousand_genomes/phase3/ALL.chrMT.phase3_callmom-v0_4.20130502.sites.vcf.gz",
        ],
    output:
        vcf="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
        vcf_md5="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.md5",
        tbi="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi",
        tbi_md5="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi.md5",
    shell:
        r"""
        bcftools merge -O z -o {output.vcf} {input.vcf}
        tabix -f {output.vcf}
        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule result_grch37_thousand_genomes_tsv:
    input:
        header="header/thousand_genomes.txt",
        vcf="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
        tbi="GRCh37/thousand_genomes/phase3/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi",
    output:
        tsv="GRCh37/thousand_genomes/phase3/ThousandGenomes.tsv",
        release_info="GRCh37/thousand_genomes/phase3/ThousandGenomes.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools view {input.vcf} \
                -e "ALT ~ 'CN'" \
            | bcftools query \
                -f "GRCh37\t%CHROM\t%POS\t%END\t%REF\t%ALT\t%AC\t%AN\t%AF\t%AFR_AF\t%AMR_AF\t%EAS_AF\t%EUR_AF\t%SAS_AF\t[%GT\t]\n" \
            | sort -u -k 2,2 -k 3,3 -k 4,4 -k 5,5 -k 6,6 -S {config[sort_memory]} \
            | awk -F $'\t' \
                'BEGIN {{
                    OFS = FS
                }}
                {{
                    delete a
                    a["0|1"] = 0
                    a["1|0"] = 0
                    a["1|1"] = 0
                    for (i=15; i<NF; ++i) {{
                        a[$i] += 1
                    }}
                    print $1,$2,$3,$4"\t",$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,a["0|1"]+a["1|0"],a["1|1"]
                }}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nThousandGenomes\tphase3\tGRCh37\t." > {output.release_info}
        """
