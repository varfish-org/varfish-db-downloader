#!/bin/bash

set -exo pipefail

HEADER=../header/hgnc.tsv
INPUT=../databases/hgnc/download/hgnc_complete_set.txt
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .txt).tsv

(
    cat $HEADER;
    tail -n +2 $INPUT \
) > $OUTPUT
