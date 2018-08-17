#!/bin/bash

set -exo pipefail

HEADER=../header/case.tsv

for INPUT in ../databases/case/download/*.ped
do
    OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .ped).tsv
    (
        cat $HEADER;
        awk -F $'\t' '
            BEGIN {
                OFS = FS
            }
            {
                pedigrees[$1] = pedigrees[$1]"{\"patient\":\""$2"\",\"father\":\""$3"\",\"mother\":\""$4"\",\"sex\":"$5",\"affected\":"$6"},"
            }
            END {
                for (pedigree in pedigrees) {
                    print pedigree,"["substr(pedigrees[pedigree], 1, length(pedigrees[pedigree])-1)"]"
                }
            }
        ' $INPUT
    ) > $OUTPUT
done

