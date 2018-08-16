#!/bin/bash

set -exo pipefail

HEADER=../header/clinvar.tsv
INSINGLE=../databases/clinvar/download/clinvar_allele_trait_pairs.single.b37.tsv.gz
INMULTI=../databases/clinvar/download/clinvar_allele_trait_pairs.multi.b37.tsv.gz

for INPUT in $INSINGLE $INMULTI
do
    OUTPUT=$(dirname $INPUT)/../$(basename $INPUT .gz)
    (
        cat $HEADER;
        gunzip -c $INPUT \
        | tail -n +2 \
        | awk -F $'\t' 'BEGIN {
            OFS=FS
            split("11 18 25 27 28 29 30 35 36 37", ff, " ")
        }
        {
            gsub("{", ""); gsub("}", "")
            for (f in ff) {
                n = split($ff[f], a, ";")
                s = "\""a[1]"\""
                for (i=2; i<=n; ++i) {
                    s = s",\""a[i]"\""
                }
                $ff[f] = "{"s"}"
            }
            print "GRCh37",$0,"1"
        }' \
        | sed 's/0000-00-00/1970-01-01/g'
    ) > $OUTPUT
done


# 11 scv
# 18 clinical_significance_ordered
# 25 review_status_ordered
# 27 all_submitters
# 28 submitters_ordered
# 29 all_traits
# 30 all_pmids
# 35 origin
# 36 xrefs
# 37 dates_ordered

#  1	chrom
#  2	pos
#  3	ref
#  4	alt
#  5	start
#  6	stop
#  7	strand
#  8	variation_type
#  9	variation_id
# 10	rcv
# 11	scv
# 12	allele_id
# 13	symbol
# 14	hgvs_c
# 15	hgvs_p
# 16	molecular_consequence
# 17	clinical_significance
# 18	clinical_significance_ordered
# 19	pathogenic
# 20	likely_pathogenic
# 21	uncertain_significance
# 22	likely_benign
# 23	benign
# 24	review_status
# 25	review_status_ordered
# 26	last_evaluated
# 27	all_submitters
# 28	submitters_ordered
# 29	all_traits
# 30	all_pmids
# 31	inheritance_modes
# 32	age_of_onset
# 33	prevalence
# 34	disease_mechanism
# 35	origin
# 36	xrefs
# 37	dates_ordered
