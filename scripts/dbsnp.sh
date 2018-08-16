#!/bin/bash

set -exo pipefail

HEADER=../header/dbsnp.tsv
INPUT=../databases/dbsnp/download/All_20180423.vcf.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .vcf.gz).tsv
REF=../downloads/hs37d5.fa

(
    cat $HEADER;
    bcftools norm \
        --fasta-ref $REF \
        --multiallelics - \
        $INPUT \
    | bcftools query \
        -f 'GRCh37\t%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
        - \
    | awk -F "\t" '!_[$2$3$4$5]++'
) > $OUTPUT
