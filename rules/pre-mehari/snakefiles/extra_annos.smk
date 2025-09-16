import sys
import tqdm


rule grch3x_refseq_exons_download:
    output:
        bed="work/download/pre-mehari/{genomebuild}/exons/refseq_exons.bed",
    params:
        assembly=lambda wildcards: DV.refseq_ref_38_assembly
        if wildcards.genomebuild == "GRCh38"
        else DV.refseq_ref_37_assembly,
        version=lambda wildcards: DV.refseq_38
        if wildcards.genomebuild == "GRCh38"
        else DV.refseq_37,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget --no-check-certificate -O $TMPDIR/chr2acc '{DV.refseq_base_url}/{params.version}/{params.assembly}/{params.assembly}_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc'
        awk 'BEGIN {{ OFS="\t" }} !/^#/ {{ print $2, $1 }}' $TMPDIR/chr2acc \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/names

        wget --no-check-certificate -O $TMPDIR/genomic.gff.gz '{DV.refseq_base_url}/{params.version}/{params.assembly}/{params.assembly}_genomic.gff.gz'
        zgrep -v '^#' $TMPDIR/genomic.gff.gz \
        | LC_ALL=C sort -k1,1 \
        > $TMPDIR/genes

        LC_ALL=C join --check-order -t $'\t' -j 1 $TMPDIR/names $TMPDIR/genes \
        | cut -f 2- \
        | sed -e 's/[=;:,]/ /g' \
        | awk 'BEGIN {{ OFS="\t" }} ($3 == "exon") {{ print $1, $4, $5 }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.bed}
        """


rule grch3x_ensembl_exons_download:
    output:
        bed="work/download/pre-mehari/{genomebuild}/exons/ensembl_exons.bed",
    params:
        ensembl=lambda wildcards: DV.ensembl_38
        if wildcards.genomebuild == "GRCh38"
        else DV.ensembl_37,
        url=lambda wildcards: DV.ensembl_38_archive_ftp
        if wildcards.genomebuild == "GRCh38"
        else DV.ensembl_37_archive_ftp,
    shell:
        r"""
        export TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT ERR

        wget --no-check-certificate -O $TMPDIR/Homo_sapiens.{wildcards.genomebuild}.{params.ensembl}.gtf.gz \
            '{params.url}/release-{params.ensembl}/gtf/homo_sapiens/Homo_sapiens.{wildcards.genomebuild}.{params.ensembl}.gtf.gz'
        zcat $TMPDIR/Homo_sapiens.{wildcards.genomebuild}.{params.ensembl}.gtf.gz \
        | awk '
            BEGIN {{ OFS="\t" }}
            ($3 == "exon") {{ print $1, $4, $5 }}' \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        | sed -e 's/[";]//g' \
        | grep '^[1-9XYM]' \
        > {output.bed}
        """


rule GRChXX_extra_annos_prepare_bed:
    input:
        "work/download/pre-mehari/{genomebuild}/exons/refseq_exons.bed",
        "work/download/pre-mehari/{genomebuild}/exons/ensembl_exons.bed",
    output:
        bed="work/download/pre-mehari/{genomebuild}/exons/refseq_ensembl_exons.bed",
    shell:
        r"""
        bedops --range 100 -m {input} > {output.bed}
        """


rule result_GRChXX_extra_annos_release_info:
    input:
        "output/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnno.tsv",
        "output/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnnoField.tsv",
    output:
        anno="output/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnno.release_info",
        annofield="output/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnnoField.release_info",
    wildcard_constraints:
        release_name="[^/]+",
    shell:
        r"""
        echo -e "table\tversion\tgenomebuild\tnull_value\nExtraAnno\t{wildcards.release_name}\t{wildcards.genomebuild}\t" > {output.anno}
        echo -e "table\tversion\tgenomebuild\tnull_value\nExtraAnnoField\t{wildcards.release_name}\t{wildcards.genomebuild}\t" > {output.annofield}
        """


rule result_GRChXX_extra_annos_tsv_step_1:
    input:
        bed="work/download/pre-mehari/{genomebuild}/exons/refseq_ensembl_exons.bed",
        cadd_snvs="work/download/annos/{genomebuild}/seqvars/cadd/{release_name}/whole_genome_SNVs_inclAnno.tsv.gz",
    output:
        tsv="work/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnno.tsv",  # this is intentionally in work
        fields="output/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnnoField.tsv",
    wildcard_constraints:
        release_name="[^/]+",
    threads: THREADS
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
    shell:
        r"""
        if [[ "{wildcards.genomebuild}" == "GRCh38" ]]; then
          cut_expr=109-112,113-117,129,130,134
        else
          cut_expr=94-102,113,114,116
        fi
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

        if [[ "{wildcards.genomebuild}" == GRCh37 ]]; then
          prefix=""
        else
          prefix="chr"
        fi

        export TMPDIR=$(mktemp -d)
        #trap "rm -rf $TMPDIR" EXIT ERR

        mkdir -p $TMPDIR/split.d
        split -n l/1028 {input.bed} $TMPDIR/split.d/exons

        export prefix
        export cut_expr
        export cols

        write-chunk()
        {{
            set -x
            set -euo pipefail
            input=$1
            output=$TMPDIR/out.d/$(basename $input)

            mkdir -p $TMPDIR/out.d

            # for the $2=$2+1: (bed file is 0-based, tbi probably 1-based) https://www.biostars.org/p/84686/
            (
                echo -e "release\tchromosome\tstart\tend\tbin\treference\talternative\tanno_data";
                awk 'BEGIN{{OFS="\t"}} {{$2=$2+1; print}}' $input \
                | tabix -R - {input.cadd_snvs} \
                | cut -f 1-4,$cut_expr \
                | awk -v "prefix=$prefix" -F $'\t' 'BEGIN {{ OFS=FS }} {{
                    printf("{wildcards.genomebuild}\t%s\t%d\t%d\t-1\t%s\t%s\t[", prefix $1, $2, $2, $3, $4);
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
            | python rules/pre-mehari/tools/ucsc_binning.py \
            | tail -n +2 \
            | sort -S 1G -k2,2g -k7,7n -o $output
        }}
        export -f write-chunk

        echo -e "release\tchromosome\tstart\tend\tbin\treference\talternative\tanno_data" \
        > {output.tsv}

        parallel -j {threads} -t 'write-chunk {{}}' ::: $TMPDIR/split.d/*
        sort -m -S 1G --parallel={threads} -k2,2g -k7,7n $TMPDIR/out.d/* \
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


rule result_GRChXX_extra_annos_tsv_step_2:
    input:
        tsv="work/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnno.tsv",  # this is intentionally in work
    output:
        tsv="output/pre-mehari/{genomebuild}/extra_annos/{release_name}/ExtraAnno.tsv",
    resources:
        runtime=os.environ.get("RUNTIME_ANNONARS_IMPORT", "48h"),
        mem_mb=MEMORY,
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


        header = [
            "release",
            "chromosome",
            "start",
            "end",
            "bin",
            "reference",
            "alternative",
            "anno_data",
        ]
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
