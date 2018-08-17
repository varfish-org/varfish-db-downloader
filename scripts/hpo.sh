#!/bin/bash

set -exo pipefail

HEADER=../header/hpo.tsv
INPUT=../databases/hpo/download/phenotype.hpoa
OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .hpoa).tsv

(
    cat $HEADER;
    tail -n +6 $INPUT
) > $OUTPUT
