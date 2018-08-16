#!/bin/bash

HEADER=../header/mim2gene_medgen.tsv
INPUT=../databases/omim/download/mim2gene_medgen
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT).tsv

(
    cat $HEADER;
    tail -n +2 $INPUT
) > $OUTPUT
