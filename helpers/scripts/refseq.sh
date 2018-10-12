#!/bin/bash

set -euo pipefail
set -x

INPUT=../data/refseq/download/ref_GRCh37.p13_top_level.gff3.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .gff3.gz).bed

(
    zcat $INPUT \
    | perl -p -i -e 's/NC_000001.10/1/,s/NC_000002.11/2/,s/NC_000003.11/3/,s/NC_000004.11/4/,s/NC_000005.9/5/,s/NC_000006.11/6/,s/NC_000007.13/7/,s/NC_000008.10/8/,s/NC_000009.11/9/,s/NC_000010.10/10/,s/NC_000011.9/11/,s/NC_000012.11/12/,s/NC_000013.10/13/,s/NC_000014.8/14/,s/NC_000015.9/15/,s/NC_000016.9/16/,s/NC_000017.10/17/,s/NC_000018.9/18/,s/NC_000019.9/19/,s/NC_000020.10/20/,s/NC_000021.8/21/,s/NC_000022.10/22/,s/NC_000023.10/X/,s/NC_000024.9/Y/,s/NC_012920.1/MT/' \
    | awk -F $'\t' 'BEGIN { OFS=FS } ($3 == "gene") { print $0; }' \
    | gff2bed \
    | cut -f 1-3 \
    | egrep '^[1-9MXY][0-9T]?'
) > $OUTPUT

