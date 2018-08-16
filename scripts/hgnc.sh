#!/bin/bash

set -exo pipefail

HEADER=../header/hgnc.tsv
INPUT=../databases/hgnc/download/hgnc_complete_set.txt.gz
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .txt).tsv

(
    cat $HEADER;
    gunzip -c $INPUT \
    | tail -n +2
) > $OUTPUT
