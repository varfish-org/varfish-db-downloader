#!/bin/bash

set -exo pipefail

HEADER=../header/thousand_genomes.tsv
INPUT=../databases/thousand_genomes/download/
OUTPUT=$INPUT/../thousand_genomes.tsv
REF=../helpers/data/reference/hs37d5.fa
FILTER=../helpers/data/coding_regions/refseq_ensembl.bed

(
    cat $HEADER;
    bcftools concat \
        -a \
        -R $FILTER \
        $INPUT/*.vcf.gz \
    | bcftools norm \
        --fasta-ref $REF \
        --multiallelics - \
        - \
    | bcftools view \
        -e "ALT ~ 'CN'" \
        - \
    | bcftools query \
        -f "GRCh37\t%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%AFR_AF\t%AMR_AF\t%EAS_AF\t%EUR_AF\t%SAS_AF\t[%GT\t]\n" \
        - \
    | awk -F $'\t' \
        'BEGIN {
            OFS = FS
        }
        {
            delete a
            a["0|1"] = 0
            a["1|0"] = 0
            a["1|1"] = 0
            for (i=14; i<NF; ++i) {
                a[$i] += 1
            }
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,a["0|1"]+a["1|0"],a["1|1"]
        }' \
    | sort -u -k 2,5 -S 80%
) > $OUTPUT

