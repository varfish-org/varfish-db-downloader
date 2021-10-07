"""Import of a database of structural variants for the ``svdbs`` app.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

import gzip
import os
import binning
import vcfpy


#: Headers used in GFF3.
GFF3_HEADERS = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


def list_to_str(liste):
    return "{" + ",".join(['"""%s"""' % x for x in liste if liste]) + "}"


class DgvGoldStandardConverter:
    """Class for the import of the DGV SV GFF3 file."""

    def __init__(self, path, name, database, output_tsv, header):
        #: Path to file to import.
        self.path = path
        #: Name of the database stored in the file.
        self.name = name
        if name not in ("dgv-gs-GRCh37", "dgv-gs-GRCh38"):
            raise Exception("Invalid database name {}".format(name))
        #: DB configuration
        self.database = database
        #: The genome release to import for
        self.genome_release = name.split("-")[-1]
        #: Seen IDs to prevent duplicate import
        self.seen_ids = set()
        #: Output file to write to
        self.fh_tsv = open(output_tsv, "w")
        #: Header for tsv
        self.header = header

    def __del__(self):
        self.fh_tsv.close()

    def convert(self):
        self.fh_tsv.write(self.header + "\n")
        with open(self.path, "rt") as inputf:
            chrom = None
            # Read file and insert into database.
            for line in inputf:
                if line:
                    line = line[:-1]
                if not line:
                    continue  # skip empty lines
                arr = line.split("\t")
                if len(arr) != len(GFF3_HEADERS):
                    raise Exception(
                        "Too few entries (%d vs %d) in line %s"
                        % (len(arr), len(GFF3_HEADERS), repr(line))
                    )
                values = dict(zip(GFF3_HEADERS, arr))
                values["attributes"] = {
                    key: value
                    for key, value in [
                        entry.split("=", 1) for entry in values["attributes"].split(";")
                    ]
                }
                if values["attributes"]["ID"] in self.seen_ids:
                    continue
                else:
                    self.seen_ids.add(values["attributes"]["ID"])
                if ("_" in values["seqid"]) or (values["seqid"] in ("", "chr")):
                    continue
                if values["seqid"].startswith("chr") and self.genome_release == "GRCh37":
                    values["seqid"] = values["seqid"][3:]
                self.write_to_tsv(values)
                # read next line
                if chrom != values["seqid"]:
                    print("Starting contig {}".format(values["seqid"]))
                chrom = values["seqid"]

    def write_to_tsv(self, values):
        """Insert record into database."""
        attributes = values["attributes"]
        pop_sum = {
            key: value
            for key, value in [x.split(" ") for x in attributes["PopulationSummary"].split(":")]
        }
        self.fh_tsv.write(
            "\t".join(
                [
                    self.genome_release,
                    values["seqid"],
                    attributes["outer_start"],
                    attributes["inner_start"],
                    attributes["inner_end"],
                    attributes["outer_end"],
                    str(
                        binning.assign_bin(
                            int(attributes["outer_start"]) - 1, int(attributes["outer_end"])
                        )
                    ),
                    attributes["ID"],
                    attributes["variant_type"],
                    attributes["variant_sub_type"],
                    attributes["num_studies"],
                    list_to_str(attributes["Studies"].split(",")),
                    attributes["num_platforms"],
                    list_to_str(attributes["Platforms"].split(",")),
                    attributes["number_of_algorithms"],
                    list_to_str(attributes["algorithms"].split(",")),
                    attributes["num_variants"],
                    attributes["num_samples"],
                    attributes.get(
                        "num_unique_samples_tested",
                        attributes.get("Number_of_unique_samples_tested"),
                    ),
                    pop_sum["African"],
                    pop_sum["Asian"],
                    pop_sum["European"],
                    pop_sum["Mexican"],
                    pop_sum["MiddleEast"],
                    pop_sum["NativeAmerican"],
                    pop_sum["NorthAmerican"],
                    pop_sum["Oceania"],
                    pop_sum["SouthAmerican"],
                    pop_sum["Admixed"],
                    pop_sum["Unknown"],
                ]
            )
            + "\n"
        )


class DgvConverter:
    """Class for the import of DGV text files."""

    def __init__(self, path, name, database, output_tsv, header):
        #: Path to file to import.
        self.path = path
        #: Name of the database stored in the file.
        self.name = name
        if name not in ("dgv-GRCh37", "dgv-GRCh38"):
            raise Exception("Invalid database name {}".format(name))
        #: DB configuration
        self.database = database
        #: The genome release to import for
        self.genome_release = self.name.split("-")[-1]
        #: Output file to write to
        self.fh_tsv = open(output_tsv, "w")
        #: Header for tsv
        self.header = header

    def __del__(self):
        self.fh_tsv.close()

    def convert(self):
        self.fh_tsv.write(self.header + "\n")
        with open(self.path, "rt") as inputf:
            chrom = None
            header = inputf.readline()[:-1].split("\t")
            # Read file and insert into database.
            linebuf = inputf.readline()
            while linebuf:
                if linebuf:
                    linebuf = linebuf[:-1]
                if not linebuf:
                    continue  # skip empty lines
                arr = linebuf.split("\t")
                if len(arr) != len(header):
                    raise Exception(
                        "Too few entries (%d vs %d) in line %s (%s)"
                        % (len(header), len(arr), repr(linebuf), repr(header))
                    )
                values = dict(zip(header, arr))
                if values["chr"].startswith("chr") and self.genome_release == "GRCh37":
                    values["chr"] = values["chr"][3:]
                elif not values["chr"].startswith("chr") and self.genome_release != "GRCh37":
                    values["chr"] = "chr%s" % values["chr"]
                if ("_" not in values["chr"]) and (values["chr"] not in ("", "chr", "N", "chrN")):
                    self.write_to_tsv(values)
                # read next line
                if chrom != values["chr"]:
                    print("Starting contig {}".format(values["chr"]))
                chrom = values["chr"]
                linebuf = inputf.readline()

    def write_to_tsv(self, values):
        """Insert record into database."""
        self.fh_tsv.write(
            "\t".join(
                [
                    self.genome_release,
                    values["chr"],
                    values["start"],
                    values["end"],
                    str(binning.assign_bin(int(values["start"]) - 1, int(values["end"]))),
                    values["variantaccession"],
                    values["varianttype"],
                    values["variantsubtype"],
                    values["reference"],
                    list_to_str(values["platform"].split(",")),
                    values["samplesize"] or "0",
                    values["observedgains"] or "0",
                    values["observedlosses"] or "0",
                ]
            )
            + "\n"
        )


class ExacCnvConverter:
    """Class for the import of ExAC CNV files."""

    def __init__(self, path, name, database, output_tsv, header):
        #: Path to file to import.
        self.path = path
        #: Name of the database stored in the file.
        self.name = name
        if name != "exac-GRCh37":
            raise Exception("Invalid database name {}".format(name))
        #: DB configuration
        self.database = database
        #: The genome release to import for
        self.genome_release = self.name.split("-")[1]
        #: Output file to write to
        self.fh_tsv = open(output_tsv, "w")
        #: Header for tsv
        self.header = header

    def __del__(self):
        self.fh_tsv.close()

    def convert(self):
        self.fh_tsv.write(self.header + "\n")
        with open(self.path, "rt") as inputf:
            chrom = None
            var_type = None
            for line in inputf:
                if line and line[-1] == "\n":
                    line = line[:-1]
                if line.startswith("track"):
                    if line.split()[1].startswith("name=delControls"):
                        var_type = "deletion"
                    else:
                        if not line.split()[1].startswith("name=dupControls"):
                            raise Exception("Unexpected track line: {}".format(line))
                        var_type = "duplication"
                else:
                    arr = line.split()
                    self.fh_tsv.write(
                        "\t".join(
                            [
                                self.genome_release,
                                arr[0][len("chr") :],
                                arr[1],
                                arr[2],
                                str(binning.assign_bin(int(arr[1]) - 1, int(arr[2]))),
                                var_type,
                                arr[3].split("-")[1],
                                arr[4],
                            ]
                        )
                        + "\n"
                    )
                    # read next line
                    if chrom != arr[0]:
                        print("Starting sv type {} on contig {}".format(var_type, arr[0]))
                    chrom = arr[0]


class ThousandGenomesConverter:
    """Class for the import of thousand genomes SVs."""

    def __init__(self, path, name, panel_path, database, output_tsv, header):
        #: Path to file to import.
        self.path = path
        #: Name of the database stored in the file.
        self.name = name
        if name != "thousand-genomes-svs-GRCh37":
            raise Exception("Invalid database name {}".format(name))
        #: DB configuration
        self.database = database
        #: The genome release to import for
        self.genome_release = self.name.split("-")[-1]
        #: Path to panel mapping
        self.panel_path = panel_path
        #: Output file to write to
        self.fh_tsv = open(output_tsv, "w")
        #: Header for tsv
        self.header = header

    def __del__(self):
        self.fh_tsv.close()

    def convert(self):
        self.fh_tsv.write(self.header + "\n")
        print("Loading panel map...")
        panel_map = self.load_panel_map()
        print("Successfully loaded panel map")
        print("Loading SV VCF file...")
        with vcfpy.Reader.from_path(self.path) as vcf_reader:
            chrom = None
            for i, record in enumerate(vcf_reader):
                if record.CHROM != chrom:
                    print("Starting on chrom %s" % record.CHROM)
                    chrom = record.CHROM
                if i % 100 == 0:
                    print("  @ {}:{:,}".format(record.CHROM, record.POS))
                self.import_sv_vcf_record(panel_map, record)

    def load_panel_map(self):
        """Load and return the panel map file."""
        result = {}
        with open(os.path.join(os.path.dirname(self.path), self.panel_path), "rt") as inputf:
            header = None
            for line in inputf:
                if line and line[-1] == "\n":
                    line = line[:-1]
                if not line:
                    continue
                if not header:
                    header = line.strip().split("\t")
                    if header != ["sample", "pop", "super_pop", "gender"]:
                        raise Exception("Unexpected header!")
                    continue
                else:
                    arr = line.split("\t")
                    result[arr[0]] = {"pop": arr[1], "super_pop": arr[2], "sex": arr[3]}
        return result

    def import_sv_vcf_record(self, panel_map, record):
        """Import the SV VCF file into the database."""
        # Counters
        super_pops = ("All", "AFR", "AMR", "EAS", "EUR", "SAS")
        num_samples = 0
        num_alleles = {key: 0 for key in super_pops}
        num_var_alleles = {key: 0 for key in super_pops}

        # Count statistics
        for call in record.calls:
            sample = call.sample
            gt = call.data.get("GT", ".")
            super_pop = panel_map[sample]["super_pop"]
            sex = panel_map[sample]["sex"]
            # Skip if genotype is no-call
            if gt == ".":
                continue
            # Count alleles contributed by this individual
            if record.CHROM == "X":
                this_alleles = 1 if sex == "male" else 2
            elif record.CHROM == "Y":
                this_alleles = 1 if sex == "male" else 0
            else:
                this_alleles = 2
            if this_alleles == 0:
                continue  # no alleles contributed by this individual
            # Increment allele counters
            num_alleles["All"] += this_alleles
            num_alleles[super_pop] += this_alleles
            num_samples += 1
            if gt in ("0|0", "0/0"):
                continue  # non-variant allele
            elif this_alleles == 1:
                num_var_alleles["All"] += 1
                num_alleles[super_pop] += 1
            elif "0" in gt:  # heterozygous, even if multiallelic (-> CNV)
                num_var_alleles["All"] += 1
                num_alleles[super_pop] += 1
            else:  # homozygous non-ref, even if multiallelic (-> CNV)
                num_var_alleles["All"] += 2
                num_alleles[super_pop] += 2

        # Perform the record creation
        self.fh_tsv.write(
            "\t".join(
                [
                    self.genome_release,
                    record.CHROM,
                    str(record.POS),
                    str(record.INFO.get("END", record.POS)),
                    str(binning.assign_bin(record.POS - 1, record.INFO.get("END", record.POS))),
                    str(record.INFO.get("CIPOS", (0, 0))[0]),
                    str(record.INFO.get("CIPOS", (0, 0))[1]),
                    str(record.INFO.get("CIEND", (0, 0))[0]),
                    str(record.INFO.get("CIEND", (0, 0))[1]),
                    record.INFO.get("SVTYPE"),
                    record.INFO.get("CS"),
                    list_to_str(record.INFO.get("MEINFO", [])),
                    str(num_samples),
                    str(num_alleles["All"]),
                    str(num_var_alleles["All"]),
                ]
                + [str(num_alleles[key]) for key in super_pops if key != "All"]
                + [str(num_var_alleles[key]) for key in super_pops if key != "All"]
            )
            + "\n"
        )


class DbVarConverter:
    """Converter for dbVar SVs."""

    def __init__(self, path, name, database, output_tsv, header):
        #: Path to file to import.
        self.path = path
        #: Name of the database stored in the file.
        self.name = name
        arr = self.name.split("-")
        if (
            len(arr) != 3
            or arr[0] != "dbvar"
            or arr[1] not in ("insertions", "duplications", "deletions")
        ):
            raise Exception("Invalid database name {}".format(name))
        #: DB configuration
        self.database = database
        #: The type of SVs to import
        self.sv_type = arr[1]
        #: The genome release to import for
        self.genome_release = arr[-1]
        #: Output file to write to
        self.fh_tsv = open(output_tsv, "w")
        #: Header for tsv
        self.header = header

    def __del__(self):
        self.fh_tsv.close()

    def convert(self):
        self.fh_tsv.write(self.header + "\n")
        with gzip.open(self.path, "rt") as inputf:
            chrom = None
            first_line = None
            second_line = None
            for i, line in enumerate(inputf):
                if line.endswith("\n"):
                    line = line[:-1]
                if not first_line:
                    first_line = line
                    if first_line != "#NR_SVs %s" % self.genome_release:
                        raise Exception("First line must be '#NR_SVs %s'" % self.genome_release)
                elif not second_line:
                    second_line = line
                    if not second_line.startswith("#chr"):
                        raise Exception("Second line must start with '#chr'")
                    header = second_line[1:].split("\t")
                else:
                    values = dict(zip(header, line.split("\t")))
                    if values["chr"] != chrom:
                        print("Starting on chrom %s" % values["chr"])
                        chrom = values["chr"]
                    if i % 10000 == 0:
                        print("  @ {}:{:,}".format(values["chr"], int(values["outermost_start"])))
                    self.create_record(values)

    def create_record(self, values):
        """Create new entry in dbVar SV table."""
        self.fh_tsv.write(
            "\t".join(
                [
                    self.genome_release,
                    values["chr"],
                    values["outermost_start"],
                    values["outermost_stop"],
                    str(
                        binning.assign_bin(
                            int(values["outermost_start"]) - 1, int(values["outermost_stop"])
                        )
                    ),
                    values["variant_count"],
                    values["variant_type"],
                    values["method"],
                    values["analysis"],
                    values["platform"],
                    values["study"],
                    list_to_str(values.get("clinical_assertion", "").split(";")),
                    list_to_str(values.get("clinvar_accessions", "").split(";")),
                    values["bin_size"],
                    values.get("min_insertion_length", ""),
                    values.get("max_insertion_length", ""),
                ]
            )
            + "\n"
        )


class GnomadSvConverter:
    """Converter for gnomAD SVs"""

    def __init__(self, path, name, database, output_tsv, header):
        #: Path to file to import.
        self.path = path
        #: Name of the database stored in the file.
        self.name = name
        arr = self.name.split("-")
        if len(arr) != 2 or arr[0] != "gnomad_sv":
            raise Exception("Invalid database name {}".format(name))
        #: DB configuration
        self.database = database
        #: The genome release to import for
        self.genome_release = arr[-1]
        #: Output file to write to
        self.fh_tsv = open(output_tsv, "w")
        #: Header for tsv
        self.header = header

    def __del__(self):
        self.fh_tsv.close()

    def convert(self):
        self.fh_tsv.write(self.header + "\n")
        prev_chrom = None
        with vcfpy.Reader.from_path(self.path) as reader:
            for i, record in enumerate(reader):
                self._create_record(record)
                if prev_chrom != record.CHROM:
                    print("Now on chrom chr%s" % record.CHROM)
                if i % 1000 == 0:
                    print("  now @ chr{}:{:,}".format(record.CHROM, record.POS))
                prev_chrom = record.CHROM

    def _create_record(self, record):
        """Create new entry in gnomAD SV table."""
        vals = [
            self.genome_release,
            record.CHROM,
            str(record.POS),
            str(record.INFO.get("END")),
            str(binning.assign_bin(record.INFO.get("END") - 1, record.POS)),
            record.REF,
            list_to_str([alt.serialize() for alt in record.ALT]),
            list_to_str(record.ID),
            record.INFO.get("SVTYPE"),
            str(record.INFO.get("SVLEN")),
            list_to_str(record.FILTER),
            list_to_str(record.INFO.get("EVIDENCE")),
            list_to_str(record.INFO.get("ALGORITHMS")),
            record.INFO.get("CHR2", ""),
            record.INFO.get("CPX_TYPE", ""),
            list_to_str(record.INFO.get("CPX_INTERVALS", [])),
            record.INFO.get("SOURCE", ""),
            record.INFO.get("STRANDS", ""),
            record.INFO.get("UNRESOLVED_TYPE", ""),
            str(record.INFO.get("PCRPLUS_DEPLETED", False)),
            str(record.INFO.get("PESR_GT_OVERDISPERSION", False)),
            list_to_str(record.INFO.get("PROTEIN_CODING_LOF", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__DUP_LOF", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__COPY_GAIN", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__DUP_PARTIAL", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__MSV_EXON_OVR", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__INTRONIC", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__INV_SPAN", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__UTR", [])),
            list_to_str(record.INFO.get("PROTEIN_CODING__NEAREST_TSS", [])),
            str(record.INFO.get("PROTEIN_CODING__INTERGENIC", False)),
            list_to_str(record.INFO.get("PROTEIN_CODING__PROMOTER", [])),
            str(record.INFO.get("AN")),
            list_to_str(record.INFO.get("AC", [])),
            list_to_str(record.INFO.get("AF", [])),
            str(record.INFO.get("N_BI_GENOS", 0)),
            str(record.INFO.get("N_HOMREF", 0)),
            str(record.INFO.get("N_HET", 0)),
            str(record.INFO.get("N_HOMALT", 0)),
            str(record.INFO.get("FREQ_HOMREF", 0.0)),
            str(record.INFO.get("FREQ_HET", 0.0)),
            str(record.INFO.get("FREQ_HOMALT", 0.0)),
            str(record.INFO.get("POPMAX_AF", 0.0)),
            str(record.INFO.get("AFR_AN")),
            list_to_str(record.INFO.get("AFR_AC", [])),
            list_to_str(record.INFO.get("AFR_AF", [])),
            str(record.INFO.get("AFR_N_BI_GENOS", 0)),
            str(record.INFO.get("AFR_N_HOMREF", 0)),
            str(record.INFO.get("AFR_N_HET", 0)),
            str(record.INFO.get("AFR_N_HOMALT", 0)),
            str(record.INFO.get("AFR_FREQ_HOMREF", 0.0)),
            str(record.INFO.get("AFR_FREQ_HET", 0.0)),
            str(record.INFO.get("AFR_FREQ_HOMALT", 0.0)),
            str(record.INFO.get("AMR_AN")),
            list_to_str(record.INFO.get("AMR_AC", [])),
            list_to_str(record.INFO.get("AMR_AF", [])),
            str(record.INFO.get("AMR_N_BI_GENOS", 0)),
            str(record.INFO.get("AMR_N_HOMREF", 0)),
            str(record.INFO.get("AMR_N_HET", 0)),
            str(record.INFO.get("AMR_N_HOMALT", 0)),
            str(record.INFO.get("AMR_FREQ_HOMREF", 0.0)),
            str(record.INFO.get("AMR_FREQ_HET", 0.0)),
            str(record.INFO.get("AMR_FREQ_HOMALT", 0.0)),
            str(record.INFO.get("EAS_AN")),
            list_to_str(record.INFO.get("EAS_AC", [])),
            list_to_str(record.INFO.get("EAS_AF", [])),
            str(record.INFO.get("EAS_N_BI_GENOS", 0)),
            str(record.INFO.get("EAS_N_HOMREF", 0)),
            str(record.INFO.get("EAS_N_HET", 0)),
            str(record.INFO.get("EAS_N_HOMALT", 0)),
            str(record.INFO.get("EAS_FREQ_HOMREF", 0.0)),
            str(record.INFO.get("EAS_FREQ_HET", 0.0)),
            str(record.INFO.get("EAS_FREQ_HOMALT", 0.0)),
            str(record.INFO.get("EUR_AN")),
            list_to_str(record.INFO.get("EUR_AC", [])),
            list_to_str(record.INFO.get("EUR_AF", [])),
            str(record.INFO.get("EUR_N_BI_GENOS", 0)),
            str(record.INFO.get("EUR_N_HOMREF", 0)),
            str(record.INFO.get("EUR_N_HET", 0)),
            str(record.INFO.get("EUR_N_HOMALT", 0)),
            str(record.INFO.get("EUR_FREQ_HOMREF", 0.0)),
            str(record.INFO.get("EUR_FREQ_HET", 0.0)),
            str(record.INFO.get("EUR_FREQ_HOMALT", 0.0)),
            str(record.INFO.get("OTH_AN")),
            list_to_str(record.INFO.get("OTH_AC", [])),
            list_to_str(record.INFO.get("OTH_AF", [])),
            str(record.INFO.get("OTH_N_BI_GENOS", 0)),
            str(record.INFO.get("OTH_N_HOMREF", 0)),
            str(record.INFO.get("OTH_N_HET", 0)),
            str(record.INFO.get("OTH_N_HOMALT", 0)),
            str(record.INFO.get("OTH_FREQ_HOMREF", 0.0)),
            str(record.INFO.get("OTH_FREQ_HET", 0.0)),
            str(record.INFO.get("OTH_FREQ_HOMALT", 0.0)),
        ]
        self.fh_tsv.write("\t".join(vals) + "\n")


DATABASES = {
    "DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3": {
        "genomebuild": "GRCh37",
        "name": "dgv-gs-GRCh37",
        "table": "DgvGoldStandardSvs",
        "release": "20160515",
        "converter": DgvGoldStandardConverter,
    },
    "DGV.GS.hg38.gff3": {
        "genomebuild": "GRCh38",
        "name": "dgv-gs-GRCh38",
        "table": "DgvGoldStandardSvs",
        "release": "20200302",
        "converter": DgvGoldStandardConverter,
    },
    "GRCh37_hg19_variants_2020-02-25.txt": {
        "genomebuild": "GRCh37",
        "name": "dgv-GRCh37",
        "table": "DgvSvs",
        "release": "20200225",
        "converter": DgvConverter,
    },
    "GRCh38_hg38_variants_2020-02-25.txt": {
        "genomebuild": "GRCh38",
        "name": "dgv-GRCh38",
        "table": "DgvSvs",
        "release": "20200225",
        "converter": DgvConverter,
    },
    "exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed": {
        "genomebuild": "GRCh37",
        "name": "exac-GRCh37",
        "table": "ExacCnv",
        "release": "release1",
        "converter": ExacCnvConverter,
    },
    "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz": {
        "genomebuild": "GRCh37",
        "name": "thousand-genomes-svs-GRCh37",
        "table": "ThousandGenomesSv",
        "release": "v8.20130502",
        "converter": ThousandGenomesConverter,
        "panel_path": "integrated_call_samples_v3.20130502.ALL.panel",
    },
    "GRCh37.nr_deletions.tsv.gz": {
        "genomebuild": "GRCh37",
        "name": "dbvar-deletions-GRCh37",
        "table": "DbVarSv",
        "release": "v0.20180824",
        "converter": DbVarConverter,
    },
    "GRCh37.nr_duplications.tsv.gz": {
        "genomebuild": "GRCh37",
        "name": "dbvar-duplications-GRCh37",
        "table": "DbVarSv:duplication",
        "release": "v0.20180824",
        "converter": DbVarConverter,
    },
    "GRCh37.nr_insertions.tsv.gz": {
        "genomebuild": "GRCh37",
        "name": "dbvar-insertions-GRCh37",
        "table": "DbVarSv:insertions",
        "release": "v0.20180824",
        "converter": DbVarConverter,
    },
    "GRCh38.nr_deletions.tsv.gz": {
        "genomebuild": "GRCh38",
        "name": "dbvar-deletions-GRCh38",
        "table": "DbVarSv:deletions",
        "release": "v0.20180824",
        "converter": DbVarConverter,
    },
    "GRCh38.nr_duplications.tsv.gz": {
        "genomebuild": "GRCh38",
        "name": "dbvar-duplications-GRCh38",
        "table": "DbVarSv:duplications",
        "release": "v0.update20180824",
        "converter": DbVarConverter,
    },
    "GRCh38.nr_insertions.tsv.gz": {
        "genomebuild": "GRCh38",
        "name": "dbvar-insertions-GRCh38",
        "table": "DbVarSv:insertions",
        "release": "v0.20180824",
        "converter": DbVarConverter,
    },
    "gnomad_v2.1_sv.sites.vcf.gz": {
        "genomebuild": "GRCh37",
        "name": "gnomad_sv-GRCh37",
        "table": "GnomAdSv",
        "release": "v2.1",
        "converter": GnomadSvConverter,
    },
}


def _write_release_info(release_info, values):
    with open(release_info, "w") as fh:
        fh.write("table\tversion\tgenomebuild\tnull_value\n")
        fh.write(
            "{}\t{}\t{}\t{}".format(values["table"], values["release"], values["genomebuild"], "")
        )


def _read_header(header):
    result = []
    with open(header, "r") as fh:
        for line in fh.readlines():
            line = line.rstrip()
            if not line:
                continue
            result.append(line)
    return "\t".join(result)


def to_tsv(path, output_tsv, release_info, header):
    database = DATABASES.get(os.path.basename(path))
    if not database:
        raise Exception("File name {} is not known".format(os.path.basename(path)))
    kwargs = {}
    if database["name"] == "thousand-genomes-svs-GRCh37":
        kwargs["panel_path"] = database["panel_path"]

    converter = database["converter"](
        path=path,
        name=database["name"],
        database=database,
        output_tsv=output_tsv,
        header=_read_header(header),
        **kwargs,
    )
    converter.convert()
    _write_release_info(
        release_info,
        {
            "genomebuild": database["genomebuild"],
            "table": database["table"],
            "release": database["release"],
            "comment": "",
        },
    )
