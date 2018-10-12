#!/bin/bash

set -euo pipefail
set -x

INENSEMBL=../data/ensembl/Homo_sapiens.GRCh37.75.bed
INREFSEQ=../data/refseq/ref_GRCh37.p13_top_level.bed
REF=../data/reference/hs37d5.fa.fai
OUTPUT=../data/coding_regions/refseq_ensembl.bed

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

(
    bedtools slop \
        -g $REF \
        -i $INREFSEQ \
        -l 5000 \
        -r 5000 \
    > $TMPDIR/refseq.bed

    bedtools slop \
        -g $REF \
        -i $INENSEMBL \
        -l 5000 \
        -r 5000 \
    > $TMPDIR/ensembl.bed

    sort-bed $TMPDIR/refseq.bed $TMPDIR/ensembl.bed \
    > $TMPDIR/genes.bed

    bedtools merge \
        -i $TMPDIR/genes.bed
) > $OUTPUT

