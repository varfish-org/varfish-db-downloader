# Download ExAC, normalize, and convert into TSV for import.


rule grch37_exac_v1_0_0_download:
    output:
        vcf="GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz",
        tbi="GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz.tbi",
    log:
        "GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz.log",
    shell:
        r"""
        wget \
            -o {log} \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/ExAC.r1.sites.vep.vcf.gz
        wget \
            -o {log} \
            -O {output.tbi} \
            https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/ExAC.r1.sites.vep.vcf.gz.tbi

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch37_exac_v1_0_0_cnv_download:
    output:
        "GRCh37/ExAC/r1/download/exac-final-cnv.gene.scores071316",
        "GRCh37/ExAC/r1/download/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
        "GRCh37/ExAC/r1/download/README.cnv_gene_scores",
        "GRCh37/ExAC/r1/download/README.cnv_bed",
        md5sums="GRCh37/ExAC/r1/download/md5sum.txt",
    shell:
        r"""
        cd $(dirname {output.md5sums})
        # mirror server
        for x in {output}; do
            wget -O $(basename $x) \
                https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exacv1_downloads/release0.3.1/cnv/$(basename $x)
        done
        # check md5 sums
        md5sum -c $(basename {output.md5sums})
        """


rule result_grch37_exac_v1_0_0_cnv_tsv:
    input:
        header="header/exaccnv.txt",
        bed="GRCh37/ExAC/r1/download/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
    output:
        tsv="GRCh37/ExAC/r1/ExacCnv.tsv",
        release_info="GRCh37/ExAC/r1/ExacCnv.release_info",
    run:
        to_tsv(input.bed, output.tsv, output.release_info, input.header)


rule grch37_exac_v1_0_0_normalize:
    input:
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        vcf="GRCh37/ExAC/r1/download/ExAC.r1.sites.vep.vcf.gz",
    output:
        vcf="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz",
        vcf_md5="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz.md5",
        tbi="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz.tbi",
        tbi_md5="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz.tbi.md5",
    shell:
        r"""
        bcftools norm \
            --threads 8 \
            --multiallelics - \
            --fasta-ref {input.reference} \
            -O z \
            -o {output.vcf} \
            {input.vcf}
        tabix -f {output.vcf}

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        popd
        """


rule result_grch37_exac_v1_0_0_tsv:
    input:
        vcf="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz",
        tbi="GRCh37/ExAC/r1/ExAC.r1.sites.vep.vcf.gz.tbi",
        header="header/exac.txt",
    output:
        release_info="GRCh37/ExAC/r1/Exac.release_info",
        tsv="GRCh37/ExAC/r1/Exac.tsv",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools reheader \
                --header <(
                    bcftools view \
                        -h \
                        {input.vcf} \
                    | sed '/##INFO=<ID=AC_Het/s/Number=A/Number=./'
                ) \
                {input.vcf} \
            | bcftools query \
                -f 'GRCh37\t%CHROM\t%POS\t%END\t%REF\t%ALT\t%AC_Adj\t%AC_AFR\t%AC_AMR\t%AC_EAS\t%AC_FIN\t%AC_NFE\t%AC_OTH\t%AC_SAS\t%AN_Adj\t%AN_AFR\t%AN_AMR\t%AN_EAS\t%AN_FIN\t%AN_NFE\t%AN_OTH\t%AN_SAS\t%AC_Hemi\t%Hemi_AFR\t%Hemi_AMR\t%Hemi_EAS\t%Hemi_FIN\t%Hemi_NFE\t%Hemi_OTH\t%Hemi_SAS\t%AC_Hom\t%Hom_AFR\t%Hom_AMR\t%Hom_EAS\t%Hom_FIN\t%Hom_NFE\t%Hom_OTH\t%Hom_SAS\t%POPMAX\n' \
                - \
            | awk -F "\t" '!_[$2$3$4$5$6]++' \
            | awk -F $'\t' \
                'BEGIN {{
                    OFS = FS
                    offset = 6
                    populations["ALL"] = offset + 1
                    populations["AFR"] = offset + 2
                    populations["AMR"] = offset + 3
                    populations["EAS"] = offset + 4
                    populations["FIN"] = offset + 5
                    populations["NFE"] = offset + 6
                    populations["OTH"] = offset + 7
                    populations["SAS"] = offset + 8
                    ac = 0
                    an = 8
                    hemi = 16
                    hom = 24
                    popmax = offset + 33
                }}
                {{
                    if ($popmax == "NA") {{
                        ac_popmax = "."
                        an_popmax = "."
                        af_popmax = "."
                        hemi_popmax = "."
                        hom_popmax = "."
                    }}
                    else {{
                        ac_popmax = $(populations[$popmax] + ac)
                        an_popmax = $(populations[$popmax] + an)
                        af_popmax = (an_popmax > 0) ? ac_popmax / an_popmax : "."
                        hemi_popmax = $(populations[$popmax] + hemi)
                        hom_popmax = $(populations[$popmax] + hom)
                    }}
                    af = ($(populations["ALL"] + an) > 0) ? $(populations["ALL"] + ac) / $(populations["ALL"] + an) : "."
                    af_afr = ($(populations["AFR"] + an) > 0) ? $(populations["AFR"] + ac) / $(populations["AFR"] + an) : "."
                    af_amr = ($(populations["AMR"] + an) > 0) ? $(populations["AMR"] + ac) / $(populations["AMR"] + an) : "."
                    af_eas = ($(populations["EAS"] + an) > 0) ? $(populations["EAS"] + ac) / $(populations["EAS"] + an) : "."
                    af_fin = ($(populations["FIN"] + an) > 0) ? $(populations["FIN"] + ac) / $(populations["FIN"] + an) : "."
                    af_nfe = ($(populations["NFE"] + an) > 0) ? $(populations["NFE"] + ac) / $(populations["NFE"] + an) : "."
                    af_oth = ($(populations["OTH"] + an) > 0) ? $(populations["OTH"] + ac) / $(populations["OTH"] + an) : "."
                    af_sas = ($(populations["SAS"] + an) > 0) ? $(populations["SAS"] + ac) / $(populations["SAS"] + an) : "."
                    $4=$4"\t"
                    print $0,ac_popmax,an_popmax,af_popmax,hemi_popmax,hom_popmax,af,af_afr,af_amr,af_eas,af_fin,af_nfe,af_oth,af_sas
                }}' 
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nExac\tr1\tGRCh37\t." > {output.release_info}
        """


rule grch37_exac_constraints_r0_3_1_download:
    output:
        "GRCh37/ExAC_constraints/r0.3.1/download/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",
    shell:
        r"""
        wget \
            -O {output} \
            ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt            
        cd $(dirname {output})
        md5sum $(basename {output}) > $(basename {output}).md5
        """


rule result_grch37_exac_constraints_r0_3_1_tsv:
    input:
        txt="GRCh37/ExAC_constraints/r0.3.1/download/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",
        header="header/exacconstraints.txt",
    output:
        tsv="GRCh37/ExAC_constraints/r0.3.1/ExacConstraints.tsv",
        release_info="GRCh37/ExAC_constraints/r0.3.1/ExacConstraints.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            tail -n +2 {input.txt} \
            | sort -u -S 80% \
            | awk -F $'\t' '
                BEGIN {{ OFS = FS }}
                {{
                    sub(/\.[0-9]+$/, "", $1)
                    print
                }}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nExacConstraints\tr0.3.1\tGRCh37\t" > {output.release_info}
        """
