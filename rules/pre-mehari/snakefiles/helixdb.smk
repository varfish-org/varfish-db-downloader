rule grchXX_helixmtdb_vcf:
    input:
        "work/download/annos/{genome_build}/seqvars/helixmtdb/{helix_v}/helixmtdb.tsv",
    output:
        "work/download/annos/{genome_build}/seqvars/helixmtdb/{helix_v}/helixmtdb.splitted.vcf",
    run:
        import vcfpy
        import csv
        import json
        import math
        from collections import defaultdict

        if wildcards.genome_build == "grch37":
            chrom = "MT"
        else:
            chrom = "chrM"

            # As stated in HelixMTdb paper by Bolze et al, 2019,
            # "Selective constraints and pathogenicity of mitochondrial DNA variants
            #  inferred from a novel database of 196,554 unrelated individuals."
        AN = 196_554

        #: Mitochondrial length
        MT_LENGTH = 16_569

        with open(input[0], "r") as fh:
            header = next(fh).strip().split("\t")
            # construct VCF header
            vcf_header = vcfpy.Header(
                lines=[
                    vcfpy.HeaderLine("fileformat", "VCFv4.2"),
                    vcfpy.ContigHeaderLine.from_mapping({"ID": chrom, "length": MT_LENGTH}),
                ],
                samples=vcfpy.SamplesInfos(sample_names=[]),
            )
            vcf_header.add_info_line(
                {
                    "ID": "AC_HOM",
                    "Number": "A",
                    "Type": "Integer",
                    "Description": "Homoplasmy allele count",
                }
            )
            vcf_header.add_info_line(
                {
                    "ID": "AN",
                    "Number": "A",
                    "Type": "Integer",
                    "Description": "Allele number",
                }
            )
            vcf_header.add_info_line(
                {
                    "ID": "AF",
                    "Number": "A",
                    "Type": "Float",
                    "Description": "Allele frequency",
                }
            )
            vcf_header.add_info_line(
                {
                    "ID": "AC_HET",
                    "Number": "A",
                    "Type": "Integer",
                    "Description": "Heteroplasmy allele count",
                }
            )
            vcf_header.add_info_line(
                {
                    "ID": "TRIALLELIC",
                    "Number": "0",
                    "Type": "Flag",
                    "Description": "Site is triallelic",
                }
            )
            vcf_writer = vcfpy.Writer.from_path(
                output[0],
                vcf_header,
            )
            variants = defaultdict(
                lambda: {
                    "reference": "",
                    "alternative": "",
                    "ac_het": 0,
                    "ac_hom": 0,
                    "triallelic": False,
                    "updates": -1,
                }
            )
            for line in fh:
                record = dict(zip(header, line.strip("\n").split("\t")))
                locus, pos = record["locus"].split(":")
                alleles = json.loads(record["alleles"])
                ref = alleles[0]
                alts = alleles[1:]
                for alt in alts:
                    if len(ref) == 1 and len(alt) == 1:
                        var_type = vcfpy.SNV
                    elif len(ref) == len(alt):
                        var_type = vcfpy.MNV
                    else:
                        var_type = vcfpy.INDEL
                    signature = "{:05d}-{}-{}".format(int(pos), ref, alt)
                    variants[signature]["reference"] = ref
                    variants[signature]["alternative"] = vcfpy.Substitution(var_type, alt)
                    variants[signature]["position"] = int(pos)
                    variants[signature]["ac_het"] += math.ceil(
                        int(record["counts_het"]) / len(alts)
                    )
                    variants[signature]["ac_hom"] += int(record["counts_hom"])
                    variants[signature]["updates"] += 1
                    variants[signature]["triallelic"] = len(alts) > 1
            for _, record in sorted(variants.items()):
                info = {
                    "AN": [AN],
                    "AF": [(record["ac_het"] + record["ac_hom"]) / AN],
                    "AC_HOM": [record["ac_hom"]],
                    "AC_HET": [record["ac_het"]],
                }
                if record["updates"] > 0 or record["triallelic"]:
                    info["TRIALLELIC"] = True
                vcf_writer.write_record(
                    vcfpy.Record(
                        chrom,
                        int(record["position"]),
                        [],
                        record["reference"],
                        [record["alternative"]],
                        None,
                        [],
                        info,
                    )
                )
            vcf_writer.close()


rule GRChXX_helixmtdb_normalize:
    input:
        vcf="work/download/annos/{genome_release}/seqvars/helixmtdb/{helix_v}/helixmtdb.splitted.vcf",
        ref="work/reference/{genome_release}/reference.fa",
    output:
        "work/download/annos/{genome_release}/seqvars/helixmtdb/{helix_v}/helixmtdb.splitted.normalized.vcf",
    shell:
        r"""
        if [ "$CI" = "true" ]; then
            echo "Skipping normalization in CI environment."
            touch {output}
            exit 0
        fi
        bcftools norm \
            -m -any \
            -c w \
            -f {input.ref} \
            -o {output} \
            {input.vcf}
        """


def input_helixmtdb_tsv(wildcards):
    return {
        "txt": f"work/download/annos/{wildcards.genome_build.lower()}/seqvars/helixmtdb/{wildcards.helix_v}/helixmtdb.splitted.normalized.vcf",
        "header": "rules/pre-mehari/header/helixmtdb.txt",
    }


rule result_GRChXX_helixmtdb_tsv:
    input:
        unpack(input_helixmtdb_tsv),
    output:
        tsv="output/pre-mehari/{genome_build}/HelixMTdb/{helix_v}/HelixMtDb.tsv",
        release_info="output/pre-mehari/{genome_build}/HelixMTdb/{helix_v}/HelixMtDb.release_info",
    shell:
        r"""
         if [ "$CI" = "true" ]; then
            echo "Skipping normalization in CI environment."
            touch {output.tsv} {output.release_info}
            exit 0
        fi
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query \
                -f "{wildcards.genome_build}\t%CHROM\t%POS\t%END\t\t%REF\t%ALT\t%AC_HOM\t%AN\t%AF\t%AC_HET\t0\t0\t0\t0\t%TRIALLELIC\n" \
                {input.txt} \
            | awk -F$'\t' 'BEGIN{{OFS=FS}}{{if($16=="."){{$16=0}}print}}'
        ) \
        | python rules/pre-mehari/tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHelixMtDb\t{wildcards.helix_v}\t{wildcards.genome_build}\t." > {output.release_info}
        """
