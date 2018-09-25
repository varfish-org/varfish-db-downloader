#!/bin/bash

set -exo pipefail

HEADERENSEMBLTOKEGG=../header/ensembltokegg.tsv
HEADERREFSEQTOKEGG=../header/refseqtokegg.tsv
HEADERKEGGINFO=../header/kegginfo.tsv

INENSEMBLTOKEGG=../databases/kegg/download/ensembltokegg.tsv
INREFSEQTOKEGG=../databases/kegg/download/refseqtokegg.tsv
INKEGGINFO=../databases/kegg/download/kegginfo.tsv

OUTENSEMBLTOKEGG=$(dirname $INENSEMBLTOKEGG)/../$(basename $INENSEMBLTOKEGG)
OUTREFSEQTOKEGG=$(dirname $INREFSEQTOKEGG)/../$(basename $INREFSEQTOKEGG)
OUTKEGGINFO=$(dirname $INKEGGINFO)/../$(basename $INKEGGINFO)

(
    cat $HEADERENSEMBLTOKEGG;
    tail -n +2 $INENSEMBLTOKEGG \
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
) > $OUTENSEMBLTOKEGG

(
    cat $HEADERREFSEQTOKEGG;
    tail -n +2 $INREFSEQTOKEGG \
    | sort -u \
    | awk -F+ 'BEGIN {
            OFS="\t"
        }
        {
            print $2,$1
        }'
) > $OUTREFSEQTOKEGG

(
    cat $HEADERKEGGINFO;
    tail -n +2 $INKEGGINFO
) > $OUTKEGGINFO