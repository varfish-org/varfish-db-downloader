#!/bin/bash

set -exo pipefail

HEADER=../header/dbsnp.tsv
INPUT=../databases/dbsnp/download/All_20180423.vcf.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .vcf.gz).tsv
REF=../downloads/hs37d5.fa

test -e $INPUT.tbi || bcftools index $INPUT

(
    cat $HEADER;
    bcftools norm \
        --fasta-ref $REF \
        --multiallelics - \
        $INPUT \
    | bcftools query \
        -f 'GRCh37\t%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
        - \
    | sort -u -k 2,5 -S 80%
) > $OUTPUT
