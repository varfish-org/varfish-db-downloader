#!/bin/bash

HEADER=../header/exac.tsv
INPUT=../databases/exac/download/ExAC.r0.3.1.sites.vep.vcf.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .vcf.gz).tsv
REF=../downloads/hs37d5.fa

# The problem is that the AC_Het field is not properly defined.
# It is set as Number=A which refers to 1 entry for each alternative allele.
# The field contains 1 if there is only 1 alternative allele but 3 if there are
# two. Setting Number=. works, but the field doesn't get split. Number=R which
# is like Number=A but including reference doesn't work for alternative alleles
# > 2. To cut a long story short, we leave out this field as it is not shown
# in the ExAC web interface anyway.

(
    cat $HEADER;
    bcftools reheader \
        --header <(
            bcftools view \
                -h \
                $INPUT \
            | sed '/##INFO=<ID=AC_Het/s/Number=A/Number=./'
        ) \
        $INPUT \
    | bcftools norm \
        --fasta-ref $REF \
        --multiallelics - \
    | bcftools query \
        -f 'GRCh37\t%CHROM\t%POS\t%REF\t%ALT\t%AC_Adj\t%AC_AFR\t%AC_AMR\t%AC_EAS\t%AC_FIN\t%AC_NFE\t%AC_OTH\t%AC_SAS\t%AN_Adj\t%AN_AFR\t%AN_AMR\t%AN_EAS\t%AN_FIN\t%AN_NFE\t%AN_OTH\t%AN_SAS\t%AC_Hemi\t%Hemi_AFR\t%Hemi_AMR\t%Hemi_EAS\t%Hemi_FIN\t%Hemi_NFE\t%Hemi_OTH\t%Hemi_SAS\t%AC_Hom\t%Hom_AFR\t%Hom_AMR\t%Hom_EAS\t%Hom_FIN\t%Hom_NFE\t%Hom_OTH\t%Hom_SAS\t%POPMAX\n' \
        - \
    | awk -F "\t" '!_[$2$3$4$5]++' \
    | awk -F $'\t' \
        'BEGIN {
            OFS = FS
            offset = 5
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
        }
        {
            if ($popmax == "NA") {
                ac_popmax = "."
                an_popmax = "."
                af_popmax = "."
                hemi_popmax = "."
                hom_popmax = "."
            }
            else {
                ac_popmax = $(populations[$popmax] + ac)
                an_popmax = $(populations[$popmax] + an)
                af_popmax = (an_popmax > 0) ? ac_popmax / an_popmax : "."
                hemi_popmax = $(populations[$popmax] + hemi)
                hom_popmax = $(populations[$popmax] + hom)
            }
            af = ($(populations["ALL"] + an) > 0) ? $(populations["ALL"] + ac) / $(populations["ALL"] + an) : "."
            af_afr = ($(populations["AFR"] + an) > 0) ? $(populations["AFR"] + ac) / $(populations["AFR"] + an) : "."
            af_amr = ($(populations["AMR"] + an) > 0) ? $(populations["AMR"] + ac) / $(populations["AMR"] + an) : "."
            af_eas = ($(populations["EAS"] + an) > 0) ? $(populations["EAS"] + ac) / $(populations["EAS"] + an) : "."
            af_fin = ($(populations["FIN"] + an) > 0) ? $(populations["FIN"] + ac) / $(populations["FIN"] + an) : "."
            af_nfe = ($(populations["NFE"] + an) > 0) ? $(populations["NFE"] + ac) / $(populations["NFE"] + an) : "."
            af_oth = ($(populations["OTH"] + an) > 0) ? $(populations["OTH"] + ac) / $(populations["OTH"] + an) : "."
            af_sas = ($(populations["SAS"] + an) > 0) ? $(populations["SAS"] + ac) / $(populations["SAS"] + an) : "."
            print $0,ac_popmax,an_popmax,af_popmax,hemi_popmax,hom_popmax,af,af_afr,af_amr,af_eas,af_fin,af_nfe,af_oth,af_sas
        }'
) > $OUTPUT
