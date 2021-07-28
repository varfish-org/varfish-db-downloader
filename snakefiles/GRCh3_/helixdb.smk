rule grchXX_helixmtdb_download:
    output:
        "{genome_build}/HelixMTdb/20200327/download/HelixMTdb_20200327.tsv",
    shell:
        r"""
        wget https://helix-research-public.s3.amazonaws.com/mito/HelixMTdb_20200327.tsv \
            -O {output}
        """


rule grchXX_helixmtdb_vcf:
    input:
        "{genome_build}/HelixMTdb/20200327/download/HelixMTdb_20200327.tsv",
    output:
        "{genome_build}/HelixMTdb/20200327/download/HelixMTdb_20200327.splitted.vcf",
    run:
        import vcfpy
        import csv
        import json
        import math
        from collections import defaultdict

        if wildcards.genome_build == "GRCh37":
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


rule GRCh37_helixmtdb_normalize:
    input:
        vcf="GRCh37/HelixMTdb/20200327/download/HelixMTdb_20200327.splitted.vcf",
        ref="GRCh37/reference/hs37d5/hs37d5.fa",
    output:
        "GRCh37/HelixMTdb/20200327/download/HelixMTdb_20200327.splitted.normalized.vcf",
    shell:
        r"""
        bcftools norm \
            -m -any \
            -c w \
            -f {input.ref} \
            -o {output} \
            {input.vcf}
        """


rule GRCh38_helixmtdb_normalize:
    input:
        vcf="GRCh38/HelixMTdb/20200327/download/HelixMTdb_20200327.splitted.vcf",
        ref="GRCh38/reference/hs38/hs38.fa",
    output:
        "GRCh38/HelixMTdb/20200327/download/HelixMTdb_20200327.splitted.normalized.vcf",
    shell:
        r"""
        bcftools norm \
            -m -any \
            -c w \
            -f {input.ref} \
            -o {output} \
            {input.vcf}
        """


rule GRChXX_helixmtdb_tsv:
    input:
        txt=(
            "{genome_build}/HelixMTdb/20200327/download/HelixMTdb_20200327.splitted.normalized.vcf"
        ),
        header="header/helixmtdb.txt",
    output:
        tsv="{genome_build}/HelixMTdb/20200327/HelixMtDb.tsv",
        release_info="{genome_build}/HelixMTdb/20200327/HelixMtDb.release_info",
    shell:
        r"""
        (
            cat {input.header} | tr '\n' '\t' | sed -e 's/\t*$/\n/g';
            bcftools query \
                -f "{wildcards.genome_build}\t%CHROM\t%POS\t%END\t\t%REF\t%ALT\t%AC_HOM\t%AN\t%AF\t%AC_HET\t0\t0\t0\t0\t%TRIALLELIC\n" \
                {input.txt} \
            | awk -F$'\t' 'BEGIN{{OFS=FS}}{{if($16=="."){{$16=0}}print}}'
        ) \
        | python tools/ucsc_binning.py \
        > {output.tsv}

        echo -e "table\tversion\tgenomebuild\tnull_value\nHelixMtDb\t20200327\t{wildcards.genome_build}\t." > {output.release_info}
        """
