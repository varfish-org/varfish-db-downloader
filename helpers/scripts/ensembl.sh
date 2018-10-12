#!/bin/bash

set -euo pipefail
set -x

INPUT=../data/ensembl/download/Homo_sapiens.GRCh37.75.gtf.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .gtf.gz).bed

(
    zcat $INPUT \
    | awk -F $'\t' 'BEGIN { OFS=FS } ($3 == "gene") { print $0; }' \
    | gff2bed \
    | cut -f 1-3 \
    | egrep '^[1-9MXY][0-9T]?'
) > $OUTPUT

