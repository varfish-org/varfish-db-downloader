#!/bin/bash

set -exo pipefail

HEADER=../header/hpo.tsv
INPUT=../databases/hpo/download/phenotype.hpoa
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .hpoa).tsv

(
    cat $HEADER;
    grep -v '^#' $INPUT
) > $OUTPUT
