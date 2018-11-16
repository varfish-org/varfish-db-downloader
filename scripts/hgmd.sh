#!/bin/bash

set -exo pipefail

HEADER=../header/hgmd_public.tsv
INPUT=../databases/hgmd/download/ensembl_variation.txt
OUTPUT=$(dirname $INPUT)/../hgmd_public.bed

(
    cat $HEADER;
	grep -F -w HGMD_MUTATION $INPUT \
    | awk \
        -F $'\t' \
        'BEGIN { OFS=FS; map[100965601] = "MT"; map[27504] = "11"; map[27505] = "21"; map[27506] = "7"; map[27507] = "Y"; map[27508] = "2"; map[27509] = "17"; map[27510] = "22"; map[27511] = "1"; map[27512] = "18"; map[27513] = "13"; map[27514] = "16"; map[27515] = "6"; map[27516] = "X"; map[27517] = "3"; map[27518] = "9"; map[27519] = "12"; map[27520] = "14"; map[27521] = "15"; map[27522] = "20"; map[27523] = "8"; map[27524] = "4"; map[27525] = "10"; map[27526] = "19"; map[27527] = "5"; } ($7 == "HGMD_MUTATION") { chrom = map[$2]; print "GRCh37", chrom, $3 - 1, $4, $8; }' \
        $INPUT
) > $OUTPUT
