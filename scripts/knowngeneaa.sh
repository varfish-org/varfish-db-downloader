#!/bin/bash

set -exo pipefail

HEADER=../header/knowngeneaa.tsv
INPUT=../databases/knowngeneaa/download/knownGene.exonAA.vcf.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .vcf.gz).tsv

test -e $INPUT.tbi || bcftools index $INPUT

(
    cat $HEADER;
    bcftools query \
        -f "%CHROM\t%POS\t%END\t%UCSC_GENE\t%ALIGNMENT\n" \
        $INPUT \
) > $OUTPUT
