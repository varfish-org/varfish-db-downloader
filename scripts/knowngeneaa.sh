#!/bin/bash

set -exo pipefail

HEADER=../header/knowngeneaa.tsv
TEMPPUT=../databases/knowngeneaa/download/knownGene.exonAA.fa.gz
INPUT=../databases/knowngeneaa/download/knownGene.exonAA.vcf.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .vcf.gz).tsv

python knowngeneaa.py \
    ../helpers/data/reference/hs37d5.fa \
    $TEMPPUT \
    --output /dev/stdout \
| bcftools sort \
    -O z \
    -o $INPUT

cd $(dirname $INPUT)
md5sum $INPUT > $INPUT.md5

tabix -f $INPUT
md5sum $INPUT.tbi > $INPUT.tbi.md5

(
    cat $HEADER;
    bcftools query \
        -f "%CHROM\t%POS\t%END\t%UCSC_GENE\t%ALIGNMENT\n" \
        $INPUT \
) > $OUTPUT
