#!/bin/bash

set -exo pipefail

HEADER=../header/gnomad.tsv
INPUT=../databases/gnomad/download/gnomad.exomes.r2.0.2.sites.vcf.bgz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .vcf.bgz).tsv
REF=../downloads/hs37d5.fa

test -e $INPUT.tbi || bcftools index $INPUT

(
    cat $HEADER;
    bcftools norm \
        --fasta-ref $REF \
        --multiallelics - \
        $INPUT \
    | bcftools query \
        -f 'GRCh37\t%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AC_AFR\t%AC_AMR\t%AC_ASJ\t%AC_EAS\t%AC_FIN\t%AC_NFE\t%AC_OTH\t%AC_SAS\t%AN\t%AN_AFR\t%AN_AMR\t%AN_ASJ\t%AN_EAS\t%AN_FIN\t%AN_NFE\t%AN_OTH\t%AN_SAS\t%Hemi\t%Hemi_AFR\t%Hemi_AMR\t%Hemi_ASJ\t%Hemi_EAS\t%Hemi_FIN\t%Hemi_NFE\t%Hemi_OTH\t%Hemi_SAS\t%Hom\t%Hom_AFR\t%Hom_AMR\t%Hom_ASJ\t%Hom_EAS\t%Hom_FIN\t%Hom_NFE\t%Hom_OTH\t%Hom_SAS\t%POPMAX\n' \
        - \
    | sort -u -k 2,5 -S 80% \
    | awk -F $'\t' \
        'BEGIN {
            OFS = FS
            offset = 5
            populations["ALL"] = offset + 1
            populations["AFR"] = offset + 2
            populations["AMR"] = offset + 3
            populations["ASJ"] = offset + 4
            populations["EAS"] = offset + 5
            populations["FIN"] = offset + 6
            populations["NFE"] = offset + 7
            populations["OTH"] = offset + 8
            populations["SAS"] = offset + 9
            ac = 0
            an = 9
            hemi = 18
            hom = 27
            popmax = offset + 37
        }
        {
            if ($popmax == ".") {
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
            af_asj = ($(populations["ASJ"] + an) > 0) ? $(populations["ASJ"] + ac) / $(populations["ASJ"] + an) : "."        
            af_eas = ($(populations["EAS"] + an) > 0) ? $(populations["EAS"] + ac) / $(populations["EAS"] + an) : "."
            af_fin = ($(populations["FIN"] + an) > 0) ? $(populations["FIN"] + ac) / $(populations["FIN"] + an) : "."
            af_nfe = ($(populations["NFE"] + an) > 0) ? $(populations["NFE"] + ac) / $(populations["NFE"] + an) : "."
            af_oth = ($(populations["OTH"] + an) > 0) ? $(populations["OTH"] + ac) / $(populations["OTH"] + an) : "."
            af_sas = ($(populations["SAS"] + an) > 0) ? $(populations["SAS"] + ac) / $(populations["SAS"] + an) : "."
            print $0,ac_popmax,an_popmax,af_popmax,hemi_popmax,hom_popmax,af,af_afr,af_amr,af_asj,af_eas,af_fin,af_nfe,af_oth,af_sas
        }'
) > $OUTPUT