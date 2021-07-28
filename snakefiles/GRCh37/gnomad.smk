# Download gnomAD exomes and genomes, normalize, and convert into TSV for import.


rule grch37_gnomad_r2_1_1_download:
    output:
        vcf="GRCh37/gnomAD_{kind}/r2.1.1/download/gnomad.{kind}.r2.1.1.sites.chr{chrom,[^.]+}.vcf.bgz",
        tbi="GRCh37/gnomAD_{kind}/r2.1.1/download/gnomad.{kind}.r2.1.1.sites.chr{chrom,[^.]+}.vcf.bgz.tbi",
    log:
        "GRCh37/gnomAD_{kind}/r2.1.1/download/gnomad.{kind}.r2.1.1.sites.chr{chrom}.vcf.bgz.log",
    shell:
        r"""
        wget \
            -o {log} \
            -O {output.vcf} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r2.1.1.sites.{wildcards.chrom}.vcf.bgz                     
        wget \
            -o {log} \
            -O {output.tbi} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/{wildcards.kind}/gnomad.{wildcards.kind}.r2.1.1.sites.{wildcards.chrom}.vcf.bgz.tbi

        pushd $(dirname {output.vcf})
        md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
        md5sum $(basename {output.tbi}) >$(basename {output.tbi}).md5
        """


rule grch37_gnomad_r2_1_1_normalize:
    input:
        reference="GRCh37/reference/hs37d5/hs37d5.fa",
        vcf="GRCh37/gnomAD_{kind}/r2.1.1/download/gnomad.{kind}.r2.1.1.sites.chr{chrom}.vcf.bgz",
    output:
        vcf="GRCh37/gnomAD_{kind}/r2.1.1/download/gnomad.{kind}.r2.1.1.sites.chr{chrom}.normalized.vcf.bgz",
        tbi="GRCh37/gnomAD_{kind}/r2.1.1/download/gnomad.{kind}.r2.1.1.sites.chr{chrom}.normalized.vcf.bgz.tbi",
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
        """


# Rule to generate TSV files from chromosomes
## Note: From the release notes (https://macarthurlab.org/2018/10/17/gnomad-v2-1/)
## Hemizygote counts
## For this release, there are no more hemizygote count annotations. Instead, our
## male counts correctly consider males as hemizygous in non-pseudo-autosomal (non-PAR)
## regions of the sex chromosomes and we provide a nonpar flag to identify all variants that are in the these regions.
rule grch37_gnomad_exomes_r2_1_1_to_tsv:
    input:
        vcf="GRCh37/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr{chrom}.normalized.vcf.bgz",
        tbi="GRCh37/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr{chrom}.normalized.vcf.bgz.tbi",
        header="header/gnomad_exomes.txt",
    output:
        tsv="GRCh37/gnomAD_exomes/r2.1.1/GnomadExomes.{chrom}.tsv",
        release_info="GRCh37/gnomAD_exomes/r2.1.1/GnomadExomes.{chrom}.release_info",
    shell:
        r"""
                (
                    cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
                    bcftools query {input.vcf} \
                        -f "\
        GRCh37\t\
        %CHROM\t\
        %POS\t\
        %END\t\
        \t\
        %REF\t\
        %ALT\t\
        %AC\t\
        %AC_afr\t\
        %AC_amr\t\
        %AC_asj\t\
        %AC_eas\t\
        %AC_fin\t\
        %AC_nfe\t\
        %AC_oth\t\
        %AC_sas\t\
        %AN\t\
        %AN_afr\t\
        %AN_amr\t\
        %AN_asj\t\
        %AN_eas\t\
        %AN_fin\t\
        %AN_nfe\t\
        %AN_oth\t\
        %AN_sas\t\
        %AF\t\
        %AF_afr\t\
        %AF_amr\t\
        %AF_asj\t\
        %AF_eas\t\
        %AF_fin\t\
        %AF_nfe\t\
        %AF_oth\t\
        %AF_sas\t\
        %AC_male\t\
        %AC_afr_male\t\
        %AC_amr_male\t\
        %AC_asj_male\t\
        %AC_eas_male\t\
        %AC_fin_male\t\
        %AC_nfe_male\t\
        %AC_oth_male\t\
        %AC_sas_male\t\
        %nhomalt\t\
        %nhomalt_afr\t\
        %nhomalt_amr\t\
        %nhomalt_asj\t\
        %nhomalt_eas\t\
        %nhomalt_fin\t\
        %nhomalt_nfe\t\
        %nhomalt_oth\t\
        %nhomalt_sas\t\
        %popmax\t\
        %AC_popmax\t\
        %AN_popmax\t\
        %AF_popmax\t\
        .\t\
        %nhomalt_popmax\t\
        %controls_AC\t\
        %controls_AC_afr\t\
        %controls_AC_amr\t\
        %controls_AC_asj\t\
        %controls_AC_eas\t\
        %controls_AC_fin\t\
        %controls_AC_nfe\t\
        %controls_AC_oth\t\
        %controls_AC_sas\t\
        %controls_AF\t\
        %controls_AF_afr\t\
        %controls_AF_amr\t\
        %controls_AF_asj\t\
        %controls_AF_eas\t\
        %controls_AF_fin\t\
        %controls_AF_nfe\t\
        %controls_AF_oth\t\
        %controls_AF_sas\t\
        %controls_AN\t\
        %controls_AN_afr\t\
        %controls_AN_amr\t\
        %controls_AN_asj\t\
        %controls_AN_eas\t\
        %controls_AN_fin\t\
        %controls_AN_nfe\t\
        %controls_AN_oth\t\
        %controls_AN_sas\t\
        %controls_AC_male\t\
        %controls_AC_afr_male\t\
        %controls_AC_amr_male\t\
        %controls_AC_asj_male\t\
        %controls_AC_eas_male\t\
        %controls_AC_fin_male\t\
        %controls_AC_nfe_male\t\
        %controls_AC_oth_male\t\
        %controls_AC_sas_male\t\
        %controls_nhomalt\t\
        %controls_nhomalt_afr\t\
        %controls_nhomalt_amr\t\
        %controls_nhomalt_asj\t\
        %controls_nhomalt_eas\t\
        %controls_nhomalt_fin\t\
        %controls_nhomalt_nfe\t\
        %controls_nhomalt_oth\t\
        %controls_nhomalt_sas\t\
        %controls_popmax\t\
        %controls_AC_popmax\t\
        %controls_AN_popmax\t\
        %controls_AF_popmax\t\
        .\t\
        %controls_nhomalt_popmax\t\
        %nonpar\n" \
                    | sort -t $'\t' -u -k 2,2 -k 3,3 -k 4,4 -k 6,6 -k 7,7 -S 80% \
                    | awk -F $'\t' \
                        'BEGIN {{
                            OFS = FS
                            nonpar = 110
                        }}
                        {{
                            # nonpar is a flag. it will take either "." if not set or "1" if set.
                            if ($nonpar == ".") {{
                                # hemi values
                                for (i=35; i<44; ++i) {{
                                    $i = "."
                                }}
                                # control hemi values
                                for (i=86; i<95; ++i) {{
                                    $i = "."
                                }}
                            }}
                            # Remove nonpar column
                            NF--
                            print $0
                        }}'
                ) \
                | python tools/ucsc_binning.py \
                > {output.tsv}

                echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadExomes\tr2.1.1\tGRCh37\t." > {output.release_info}
        """


# Rule to generate TSV files from chromosomes
rule grch37_gnomad_genomes_r2_1_1_to_tsv:
    input:
        vcf="GRCh37/gnomAD_genomes/r2.1.1/download/gnomad.genomes.r2.1.1.sites.chr{chrom}.normalized.vcf.bgz",
        tbi="GRCh37/gnomAD_genomes/r2.1.1/download/gnomad.genomes.r2.1.1.sites.chr{chrom}.normalized.vcf.bgz.tbi",
        header="header/gnomad_genomes.txt",
    output:
        tsv="GRCh37/gnomAD_genomes/r2.1.1/GnomadGenomes.{chrom}.tsv",
        release_info="GRCh37/gnomAD_genomes/r2.1.1/GnomadGenomes.{chrom}.release_info",
    shell:
        r"""
                (
                    cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
                    bcftools query {input.vcf} \
                        -f "\
        GRCh37\t\
        %CHROM\t\
        %POS\t\
        %END\t\
        \t\
        %REF\t\
        %ALT\t\
        %AC\t\
        %AC_afr\t\
        %AC_amr\t\
        %AC_asj\t\
        %AC_eas\t\
        %AC_fin\t\
        %AC_nfe\t\
        %AC_oth\t\
        %AN\t\
        %AN_afr\t\
        %AN_amr\t\
        %AN_asj\t\
        %AN_eas\t\
        %AN_fin\t\
        %AN_nfe\t\
        %AN_oth\t\
        %AF\t\
        %AF_afr\t\
        %AF_amr\t\
        %AF_asj\t\
        %AF_eas\t\
        %AF_fin\t\
        %AF_nfe\t\
        %AF_oth\t\
        %AC_male\t\
        %AC_afr_male\t\
        %AC_amr_male\t\
        %AC_asj_male\t\
        %AC_eas_male\t\
        %AC_fin_male\t\
        %AC_nfe_male\t\
        %AC_oth_male\t\
        %nhomalt\t\
        %nhomalt_afr\t\
        %nhomalt_amr\t\
        %nhomalt_asj\t\
        %nhomalt_eas\t\
        %nhomalt_fin\t\
        %nhomalt_nfe\t\
        %nhomalt_oth\t\
        %popmax\t\
        %AC_popmax\t\
        %AN_popmax\t\
        %AF_popmax\t\
        .\t\
        %nhomalt_popmax\t\
        %controls_AC\t\
        %controls_AC_afr\t\
        %controls_AC_amr\t\
        %controls_AC_asj\t\
        %controls_AC_eas\t\
        %controls_AC_fin\t\
        %controls_AC_nfe\t\
        %controls_AC_oth\t\
        %controls_AF\t\
        %controls_AF_afr\t\
        %controls_AF_amr\t\
        %controls_AF_asj\t\
        %controls_AF_eas\t\
        %controls_AF_fin\t\
        %controls_AF_nfe\t\
        %controls_AF_oth\t\
        %controls_AN\t\
        %controls_AN_afr\t\
        %controls_AN_amr\t\
        %controls_AN_asj\t\
        %controls_AN_eas\t\
        %controls_AN_fin\t\
        %controls_AN_nfe\t\
        %controls_AN_oth\t\
        %controls_AC_male\t\
        %controls_AC_afr_male\t\
        %controls_AC_amr_male\t\
        %controls_AC_asj_male\t\
        %controls_AC_eas_male\t\
        %controls_AC_fin_male\t\
        %controls_AC_nfe_male\t\
        %controls_AC_oth_male\t\
        %controls_nhomalt\t\
        %controls_nhomalt_afr\t\
        %controls_nhomalt_amr\t\
        %controls_nhomalt_asj\t\
        %controls_nhomalt_eas\t\
        %controls_nhomalt_fin\t\
        %controls_nhomalt_nfe\t\
        %controls_nhomalt_oth\t\
        %controls_popmax\t\
        %controls_AC_popmax\t\
        %controls_AN_popmax\t\
        %controls_AF_popmax\t\
        .\t\
        %controls_nhomalt_popmax\t\
        %nonpar\n" \
                    | sort -t $'\t' -u -k 2,2 -k 3,3 -k 4,4 -k 6,6 -k 7,7 -S 80% \
                    | awk -F $'\t' \
                        'BEGIN {{
                            OFS = FS
                            nonpar = 100
                        }}
                        {{
                            # nonpar is a flag. it will take either "." if not set or "1" if set.
                            if ($nonpar == ".") {{
                                # hemi values
                                for (i=32; i<40; ++i) {{
                                    $i = "."
                                }}
                                # control hemi values
                                for (i=78; i<86; ++i) {{
                                    $i = "."
                                }}
                            }}
                            # Remove nonpar column
                            NF--
                            print $0
                        }}'
                ) \
                | python tools/ucsc_binning.py \
                > {output.tsv}

                echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadGenomes\tr2.1.1\tGRCh37\t." > {output.release_info}
        """


rule grch37_gnomad_genomes_sv_v2_1_download:
    output:
        vcf="GRCh37/gnomAD_SV/v2.1/download/gnomad_v2.1_sv.sites.vcf.gz",
        tbi="GRCh37/gnomAD_SV/v2.1/download/gnomad_v2.1_sv.sites.vcf.gz.tbi",
    shell:
        r"""
        wget -O {output.vcf} https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz
        wget -O {output.tbi} https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz.tbi

        cd $(dirname {output.vcf})
        md5sum gnomad_v2.1_sv.sites.vcf.gz >gnomad_v2.1_sv.sites.vcf.gz.md5
        md5sum gnomad_v2.1_sv.sites.vcf.gz.tbi >gnomad_v2.1_sv.sites.vcf.gz.md5
        """


rule grch37_gnomad_genomes_sv_v2_1_tsv:
    input:
        "GRCh37/gnomAD_SV/v2.1/download/gnomad_v2.1_sv.sites.vcf.gz.tbi",
        header="header/gnomadsv.txt",
        vcf="GRCh37/gnomAD_SV/v2.1/download/gnomad_v2.1_sv.sites.vcf.gz",
    output:
        tsv="GRCh37/gnomAD_SV/v2.1/GnomAdSv.tsv",
        release_info="GRCh37/gnomAD_SV/v2.1/GnomAdSv.release_info",
    run:
        to_tsv(input.vcf, output.tsv, output.release_info, input.header)


rule grch37_gnomad_constraints_v2_1_1_download:
    output:
        "GRCh37/gnomAD_constraints/v2.1.1/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
    shell:
        r"""
        wget \
            -O {output} \
            https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

        cd $(dirname {output})
        md5sum $(basename {output}) > $(basename {output}).md5
        """


rule grch37_gnomad_constraints_v2_1_1_tsv:
    input:
        txt="GRCh37/gnomAD_constraints/v2.1.1/download/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        header="header/gnomadconstraints.txt",
    output:
        tsv="GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.tsv",
        release_info="GRCh37/gnomAD_constraints/v2.1.1/GnomadConstraints.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            zcat {input.txt} \
            | tail -n +2 \
            | sort -u -S 80% \
            | awk -F $'\t' '
                BEGIN {{ OFS = FS }}
                {{
                    for (i=1; i<=NF; ++i) {{
                        if ($i == "NA") {{
                            $i = ""
                        }}
                    }}
                    print
                }}'
        ) > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nGnomadConstraints\tv2.1.1\tGRCh37\t" > {output.release_info}
        """
