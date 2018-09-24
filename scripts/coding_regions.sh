#!/bin/bash

set -euo pipefail
set -x

INPUT=../databases/ensembl/download/Homo_sapiens.GRCh37.75.gtf.gz
OUTPUT=../databases/coding_regions/refseq_ensembl.bed

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

(
    bedtools slop \
        -g ../downloads/hs37d5.fa.fai \
        -i ../databases/ensembl/Homo_sapiens.GRCh37.75.bed \
        -l 5000 \
        -r 5000 \
    > $TMPDIR/refseq.bed

    bedtools slop \
        -g ../downloads/hs37d5.fa.fai \
        -i ../databases/refseq/ref_GRCh37.p13_top_level.bed \
        -l 5000 \
        -r 5000 \
    > $TMPDIR/ensembl.bed

    sort-bed $TMPDIR/refseq.bed $TMPDIR/ensembl.bed \
    > $TMPDIR/genes.bed

    bedtools merge \
        -i $TMPDIR/genes.bed
) > $OUTPUT

