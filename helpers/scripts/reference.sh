#!/bin/bash

set -exo pipefail

test -f ../data/reference/hs37d5.fa.fai || samtools faidx ../data/reference/hs37d5.fa

