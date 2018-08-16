#!/bin/bash

HEADERGENETOKEGG=../header/genetokegg.tsv
HEADERKEGGINFO=../header/kegginfo.tsv
INGENETOKEGG=../databases/kegg/download/genetokegg.tsv
INKEGGINFO=../databases/kegg/download/kegginfo.tsv
OUTGENETOKEGG=$(dirname $INGENETOKEGG)/../$(basename $INGENETOKEGG)
OUTKEGGINFO=$(dirname $INKEGGINFO)/../$(basename $INKEGGINFO)

(
    cat $HEADERGENETOKEGG;
    tail -n +2 $INGENETOKEGG \
    | sort -u \
    | awk -F $'\t' 'BEGIN {
            OFS=FS
        }
        {
            if ($2 == "n/a") {
                next
            }
            print $2,$1
        }'
) > $OUTGENETOKEGG

(
    cat $HEADERKEGGINFO;
    tail -n +2 $INKEGGINFO
) > $OUTKEGGINFO