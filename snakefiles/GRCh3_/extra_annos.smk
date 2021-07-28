import sys
import tqdm


rule GRChXX_extra_annos_release_info:
    input:
        "{genome_build}/extra-annos/20200704/ExtraAnno.tsv",
        "{genome_build}/extra-annos/20200704/ExtraAnnoField.tsv",
    output:
        anno="{genome_build}/extra-annos/20200704/ExtraAnno.release_info",
        annofield="{genome_build}/extra-annos/20200704/ExtraAnnoField.release_info",
    shell:
        r"""
        echo -e "table\tversion\tgenomebuild\tnull_value\nExtraAnno\t20200704\t{wildcards.genome_build}\t" > {output.anno}
        echo -e "table\tversion\tgenomebuild\tnull_value\nExtraAnnoField\t20200704\t{wildcards.genome_build}\t" > {output.annofield}
        """


rule GRChXX_extra_annos_tsv_step_1:
    input:
        bed="{genome_build}/extra-annos/20200704/download/refseq_ensembl_exons.bed",
        cadd_snvs="{genome_build}/extra-annos/20200704/download/whole_genome_SNVs_inclAnno.tsv.gz",
    output:
        tsv="{genome_build}/extra-annos/20200704/download/ExtraAnno.tsv",
        fields="{genome_build}/extra-annos/20200704/ExtraAnnoField.tsv",
    shell:
        r"""
        cut_expr=94-102,113,114,116
        cols=$({{ zcat {input.cadd_snvs} || true; }} | head -n 2 | tail -n 1| cut -f $cut_expr)

        (
            echo -e "field\tlabel";
            i=0
            for label in $cols; do
                let "i=$i+1";
                echo -e "$i\t$label";
            done
        ) \
        | sed -e 's/PHRED/CADD-PHRED/g' \
        > {output.fields}

        (
            echo -e "release\tchromosome\tstart\tend\tbin\treference\talternative\tanno_data";
            tabix -R {input.bed} {input.cadd_snvs} \
            | cut -f 1-4,$cut_expr \
            | awk -F $'\t' 'BEGIN {{ OFS=FS }} {{
                printf("{wildcards.genome_build}\t%s\t%d\t%d\t-1\t%s\t%s\t[", $1, $2, $2, $3, $4);
                for (i = 5; i <= NF; i++) {{
                    x = ($i == "NA" || $i == ".") ? "null" : $i;
                    if (i != 5) {{
                        printf(",%s", x);
                    }} else {{
                        printf("%s", x);
                    }}
                }}
                printf("]\n");
            }}';
        ) \
        | python tools/ucsc_binning.py \
        | sort -k2,2g -k7,7n \
        | uniq \
        >> {output.tsv}
        """


class DecodeDotAsNull(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        if obj == ".":
            return None
        else:
            return super().object_hook(obj)


rule GRChXX_extra_annos_tsv_step_2:
    input:
        tsv="{genome_build}/extra-annos/20200704/download/ExtraAnno.tsv",
    output:
        tsv="{genome_build}/extra-annos/20200704/ExtraAnno.tsv",
    run:
        import json
        import csv


        def merge_rows(rows):
            result = rows[0][:7] + [[]]
            try:
                values = [
                    json.loads(row[7].replace(",.,", ",null,"), cls=DecodeDotAsNull) for row in rows
                ]
            except:
                print("\n%s\n" % rows, file=sys.stderr)
                raise
            for i in range(len(values[0])):
                value = None
                for xs in values:
                    if xs[i] is not None:
                        if value is None:
                            value = xs[i]
                        value = max(value, xs[i])
                result[-1].append(value)
            result[-1] = json.dumps(result[-1])
            # print("MERGED\n  %s\nTO\n  %s" % (rows, result), file=sys.stderr)
            return result


        header = ["release", "chromosome", "start", "end", "bin", "ref", "alt", "anno_data"]
        rows = []
        print("Reading from %s" % input.tsv, file=sys.stderr)
        with open(input.tsv, "rt") as inputf:
            reader = csv.reader(inputf, delimiter="\t")
            with open(output.tsv, "wt") as outputf:
                writer = csv.writer(outputf, delimiter="\t")
                writer.writerow(header)

                for row in tqdm.tqdm(reader):
                    if row[0] == header[0]:
                        continue  # skip
                    else:
                        if not rows or row[:7] == rows[0][:7]:
                            rows.append(row)
                        else:
                            writer.writerow(merge_rows(rows))
                            rows = [row]
                if rows:
                    writer.writerow(merge_rows(rows))


rule GRChXX_extra_annos_download:
    output:
        vcf="{genome_build}/extra-annos/20200704/download/whole_genome_SNVs_inclAnno.tsv.gz",
        tbi="{genome_build}/extra-annos/20200704/download/whole_genome_SNVs_inclAnno.tsv.gz.tbi",
    shell:
        r"""
        wget \
            https://kircherlab.bihealth.org/download/CADD/v1.6/{genome_build}/whole_genome_SNVs_inclAnno.tsv.gz \
            -O {output.vcf}
        wget \
            https://kircherlab.bihealth.org/download/CADD/v1.6/{genome_build}/whole_genome_SNVs_inclAnno.tsv.gz.tbi \
            -O {output.tbi}
        """


rule GRChXX_extra_annos_prepare_bed:
    output:
        bed="{genome_build}/extra-annos/20200704/download/refseq_ensembl_exons.bed",
    shell:
        r"""
        bedops -m \
            /fast/projects/cubit/20.05/static_data/exon_list/ENSEMBL/75/{genome_build}/ENSEMBL_75_plus_100bp.bed \
            /fast/projects/cubit/20.05/static_data/exon_list/RefSeq/81/{genome_build}/RefSeq_81_release_105_plus_100bp.bed \
        > {output.bed}
        """
