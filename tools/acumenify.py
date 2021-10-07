#!/usr/bin/env python

import argparse
import collections
import enum
import json
from pathlib import Path
import re
import sys
import time
import typing

import attr
import cattr
from tqdm import tqdm

from logzero import logger


#: Disable tqdm if no terminal is attached.
DISABLE_TQDM = not sys.stdout.isatty()

#: Human chromosomes.
CHROMS = tuple(list(map(str, range(1, 23))) + ["X", "Y", "MT"])


class ExtractionError(Exception):
    """Raised in the case of extraction problems."""


@attr.s(frozen=True, auto_attribs=True)
class ExtractArgs:
    #: Path to ``import_versions.tsv`` file.
    path_import_versions: str
    #: Path to write statistics files to.
    path_stats: str
    #: Optionally, number of lines to read.
    line_limit: typing.Optional[int] = None
    #: Optionally, name of table group to process.
    only: typing.Tuple[str] = None
    #: Harmonize chromosomes to this prefix.
    harmonize_chroms: str = "chr"
    #: Harmonize chrMT to this name (without prefix).
    harmonize_chrmt: str = "M"


@attr.s(frozen=True, auto_attribs=True)
class ReportArgs:
    #: Path to write report file to.
    path_report: str
    #: Path to ``import_versions.tsv`` file.
    paths_stats: typing.List[str]
    #: Harmonize chromosomes to this prefix.
    harmonize_chroms: str = "chr"
    #: Harmonize chrMT to this name (without prefix).
    harmonize_chrmt: str = "M"


@attr.s(frozen=True, auto_attribs=True)
class ImportVersion:
    build: str
    table_group: str
    version: str

    def get_key(self):
        """Return separated key for import version"""
        return "%s.%s" % (self.build, self.table_group)


@attr.s(frozen=True, auto_attribs=True)
class StatsKey:
    table_group: str
    table: str
    chrom: str


class ColumnType(enum.Enum):
    """Enumeration for colum types."""

    UNKNOWN = 1
    ENUM = 2
    STRING = 3
    INT = 4
    FLOAT = 5
    DNA = 6


#: Header known to hold chromosomes.
HEADER_CHROM = "chromosome"

#: Regular expression for a DNA string.
RE_DNA = r"^[ACGTNacgtn]+$"


@attr.s(frozen=True, auto_attribs=True)
class ChromHarmonizer:
    harmonize_chroms: str
    harmonize_chrmt: str

    def apply(self, value):
        if value == ".":
            return value
        value = self._apply_chr(value)
        m_name = self._apply_chr("M")
        mt_name = self._apply_chr("MT")
        if value in (m_name, mt_name):
            return self._apply_chr(self.harmonize_chrmt)
        else:
            return value

    def _apply_chr(self, value):
        value_chr = value.startswith("chr")
        should_chr = self.harmonize_chroms == "chr"
        if not value_chr and should_chr:
            return "chr" + value
        elif value_chr and not should_chr:
            return value[3:]
        else:
            return value


class ColAggregator:
    """Helper class for aggregating column values."""

    def __init__(self, name, chrom_harmonizer, guess_len=10000, max_enum_size=100):
        #: Column name
        self.name = name
        #: Chromozome harmoniser.
        if self.name == HEADER_CHROM:
            self.fn_harmonize = chrom_harmonizer.apply
        else:
            self.fn_harmonize = lambda x: x  # identity
        #: Number of records to base guess on.
        self.guess_len = guess_len
        #: Maximal enumeration length.
        self.max_enum_size = max_enum_size
        #: The guessed type.
        self.type_ = ColumnType.UNKNOWN
        if name == HEADER_CHROM:
            self.type_ = ColumnType.ENUM
        #: Values seen so far.
        self.values = {}
        #: Number of records read so far.
        self.counter = 0

    def process(self, value):
        self.counter += 1
        if self.counter >= self.guess_len:
            self.finish()
        value = self.fn_harmonize(value)
        if self.values is not None:
            self.values.setdefault(value, 0)
            self.values[value] += 1

    def finish(self):
        """Force guessing of value if not done so far."""
        if self.type_ != ColumnType.UNKNOWN:
            return
        if not self.values:
            return
        for f, t in ((float, ColumnType.INT), (float, ColumnType.FLOAT)):
            for k in self.values.keys():
                try:
                    f(k)
                except ValueError:
                    break
            else:
                self.type_ = t
                break
        else:
            if all(re.match(RE_DNA, v) for v in self.values.keys()):
                self.type_ = ColumnType.DNA
            elif len(self.values) <= self.max_enum_size:
                self.type_ = ColumnType.ENUM
            else:
                self.type_ = ColumnType.STRING
        if self.type_ != ColumnType.ENUM:
            self.values = None

    def to_dict(self):
        """Return dict with results."""
        result = {"type": self.type_.name}
        if self.type_ == ColumnType.ENUM:
            result["values"] = list(sorted(self.values))
        return result


class Aggregator:
    """Helper class for aggregating values"""

    def __init__(self, header, chrom_harmonizer):
        #: Chromozome harmoniser.
        self.chrom_harmonizer = chrom_harmonizer
        #: Count by chromosome, if any, else None.
        self.by_chrom = {}
        #: Aggregate statistics for each column.
        self.by_column = {col: ColAggregator(col, chrom_harmonizer) for col in header}

    def process(self, record: typing.Dict):
        chrom = self.chrom_harmonizer.apply(record.get(HEADER_CHROM, "."))
        self.by_chrom.setdefault(chrom, 0)
        self.by_chrom[chrom] += 1
        for k, v in record.items():
            self.by_column[k].process(v)

    def finish(self):
        for a in self.by_column.values():
            a.finish()
        if "." not in self.by_chrom:
            self.by_chrom["."] = sum(self.by_chrom.values())


def do_extraction(path_base: Path, record: ImportVersion, args: ExtractArgs):
    logger.info("Looking into %s", record)
    # Collect files for each table.
    by_table = {}
    path_dir = path_base / record.build / record.table_group / record.version
    for path_file in path_dir.glob("*.tsv"):
        num_dots = str(path_file.name).count(".")
        if num_dots in (1, 2):
            by_table.setdefault(path_file.name.split(".")[0], []).append(path_file)
        else:
            logger.error("No or more than two dots in filename %s", path_file.name)
            return 1
    # Process for each table.
    result = {}
    for table, paths in by_table.items():
        header = None
        logger.info("processing table %s in %s", table, record.table_group)
        agg = None
        for path in sorted(set(paths)):
            logger.info(
                "  processing file %s%s",
                path,
                " (at most %s lines)" % args.line_limit if args.line_limit else "",
            )
            this_header = None
            with path.open("rt") as inputf:
                prev = time.time()
                modulo = 200_000
                for lineno, line in tqdm(enumerate(inputf), unit="rec", disable=DISABLE_TQDM):
                    if DISABLE_TQDM and lineno and lineno % modulo == 0:
                        curr = time.time()
                        lines_per_sec = modulo / (curr - prev)
                        prev = curr
                        logger.info(
                            "... processing ... at line %s (%.2f per sec)",
                            "{:,}".format(lineno),
                            lines_per_sec,
                        )
                    if args.line_limit and lineno > args.line_limit:
                        logger.debug("Stopping at line limit of %s", args.line_limit)
                        break
                    arr = line.strip().split("\t")
                    if not this_header:
                        this_header = arr
                        if not header:
                            header = this_header
                            agg = Aggregator(
                                header, ChromHarmonizer(args.harmonize_chroms, args.harmonize_chrmt)
                            )
                        else:
                            if header != this_header:
                                raise ExtractionError(
                                    "Incompatible headers, this: %s, first: %s", this_header, header
                                )
                    else:
                        if len(this_header) != len(arr):
                            raise ExtractionError(
                                "Record %s has %d fields but header had %s",
                                lineno,
                                len(arr),
                                len(header),
                            )
                        agg.process(dict(zip(header, arr)))
        agg.finish()
        result[table] = {
            "$schema_version": "0.1.0",
            "import_version": vars(record),
            "by_chrom": agg.by_chrom,
            "by_col": {k: v.to_dict() for k, v in agg.by_column.items()},
        }
    return result


def extract(args: ExtractArgs):
    """Run extraction."""
    logger.info("Running extraction with args\n\n%s", json.dumps(vars(args), indent=2))
    path_import_versions = Path(args.path_import_versions)
    path_stats = Path(args.path_stats)

    if not path_import_versions.exists():
        logger.error("Input file '%s' does not exist", path_import_versions)
        return 1

    path_stats.mkdir(parents=True, exist_ok=True)
    path_base = path_import_versions.parent

    # Load import_versions.tsv file.
    header = None
    stats = {}
    with path_import_versions.open("rt") as inputf:
        for line in inputf:
            arr = line.strip().split("\t")
            if not header:
                expected = ["build", "table_group", "version"]
                header = arr
                if header != expected:
                    logger.error("Header unexpected, WAS: %s, expected: %s", header, expected)
            else:
                record = ImportVersion(**dict(zip(header, arr)))
                if record in stats:
                    logger.debug("Ignoring second occurence of %s", record)
                elif args.only and record.table_group not in args.only:
                    logger.info("Table group %s is not in --only=%s", record.table_group, args.only)
                else:
                    logger.debug("Handling %s", record)
                    stats[record] = True
                    try:
                        res = do_extraction(path_base, record, args)
                    except ExtractionError as e:
                        logger.error("Problem with extraction, giving up: %s", e)
                        return 1
                    for table, values in res.items():
                        path_out = path_stats / ("%s.%s.json" % (record.get_key(), table))
                        logger.info("Writing results to %s", path_out)
                        with path_out.open("wt") as outputf:
                            print(json.dumps(values, indent=2), file=outputf)


class ReportLevel(enum.Enum):
    INFO = 1
    WARNING = 2
    ERROR = 3


@attr.s(frozen=True, auto_attribs=True)
class ReportMessage:
    level: ReportLevel
    stats: str
    table_group: str
    table: str
    msg: str


@attr.s(frozen=True, auto_attribs=True)
class ReportKey:
    stats: typing.Optional[str] = None
    release: typing.Optional[str] = None
    table_group: typing.Optional[str] = None
    table: typing.Optional[str] = None

    def matches(self, other):
        def lower(s):
            if not s:
                return s
            else:
                return s.lower()

        matches = True
        for k in ("stats", "release", "table_group", "table"):
            if getattr(other, k):
                matches = lower(getattr(self, k)) == lower(getattr(other, k))
        return matches


@attr.s(frozen=True, auto_attribs=True)
class SanityCheck:
    report_key: ReportKey
    check: typing.Callable[[ReportKey, typing.Dict, ReportArgs], typing.List[ReportMessage]]


def sc_allchroms(
    key: ReportKey, report: typing.Dict, args: ReportArgs
) -> typing.List[ReportMessage]:
    """Sanity checks that tests that there is data for all chromosomes."""
    h = ChromHarmonizer(args.harmonize_chroms, args.harmonize_chrmt).apply
    result = []

    exceptions = {
        ReportKey(table_group="acmg"): [h(c) for c in CHROMS],
        ReportKey(table_group="DGV", table="DgvGoldStandardSvs"): [h("MT")],
        ReportKey(table_group="DGV", table="DgvSvs"): [h("MT")],
        ReportKey(table_group="ensembl_regulatory", table="EnsemblRegulatoryFeature"): [h("MT")],
        ReportKey(table_group="ensembltogenesymbol", table="EnsemblToGeneSymbol"): [
            h(c) for c in CHROMS
        ],
        ReportKey(table_group="ExAC_constraints", table="ExacConstraints"): [h("MT")],
        ReportKey(table_group="ExAC", table="ExacCnv"): [h("MT"), h("X"), h("Y")],
        ReportKey(table_group="ExAC", table="Exac"): [h("MT")],
        ReportKey(table_group="extra_annos", table="ExtraAnno"): [h("MT")],
        ReportKey(table_group="extra-annos", table="ExtraAnnoField"): [h(c) for c in CHROMS],
        ReportKey(table_group="gnomAD_constraints", table="GnomadConstraints"): [h("MT")],
        ReportKey(table_group="gnomAD_exomes", table="GnomadExomes"): [h("MT")],
        ReportKey(table_group="gnomAD_genomes", table="GnomadGenomes"): [h("MT"), h("Y")],
        ReportKey(table_group="gnomAD_SV", table="GnomAdSv"): [h("MT")],
        ReportKey(table_group="HelixMTdb", table="HelixMtDb"): [h(c) for c in CHROMS if c != "MT"],
        ReportKey(table_group="hgmd_public", table="HgmdPublicLocus"): [h("MT")],
        ReportKey(table_group="kegg"): [h(c) for c in CHROMS],
        ReportKey(table_group="knowngeneaa", table="KnowngeneAA"): [h("MT")],
        ReportKey(table_group="mgi", table="MgiHomMouseHumanSequence"): [h(c) for c in CHROMS],
        ReportKey(table_group="mim2gene", table="Mim2geneMedgen"): [h(c) for c in CHROMS],
        ReportKey(table_group="MITOMAP", table="Mitomap"): [h(c) for c in CHROMS if c != "MT"],
        ReportKey(table_group="mtDB", table="MtDb"): [h(c) for c in CHROMS if c != "MT"],
        ReportKey(table_group="refseqtoensembl"): [h(c) for c in CHROMS],
        ReportKey(table_group="refseqtogenesymbol"): [h(c) for c in CHROMS],
        ReportKey(table_group="tads_hesc", table="TadInterval"): [h("MT"), h("Y")],
        ReportKey(table_group="tads_imr90", table="TadInterval"): [h("MT"), h("Y")],
        ReportKey(table_group="thousand_genomes", table="ThousandGenomesSv"): [h("MT"), h("Y")],
        ReportKey(table_group="vista", table="VistaEnhancer"): [h("MT"), h("Y")],
        ReportKey(table_group="refseq_genes", release="39", table="GeneInterval"): [h("MT")],
    }

    chroms = set(map(h, CHROMS)) | {"."}
    missing = chroms - report["by_chrom"].keys()
    for exc_key, exc_chr in exceptions.items():
        if key.matches(exc_key):
            missing -= set(exc_chr)
    if missing:
        result.append(
            ReportMessage(
                ReportLevel.ERROR,
                key.stats,
                key.table_group,
                key.table,
                "No count for %d chromosomes: %s" % (len(missing), list(sorted(missing))),
            )
        )

    return result


#: Defined simple sanity checks that work on one report at a time."""
SANITY_CHECKS = (SanityCheck(ReportKey(), sc_allchroms),)


@attr.s(frozen=True, auto_attribs=True)
class KnownIssue:
    dataset: str
    table_group: str
    table: str
    message: str

    def matches(self, msg: ReportMessage):
        return all(
            (
                self.dataset == msg.stats,
                self.table_group == msg.table_group,
                self.table == msg.table,
                re.match(self.message, msg.msg),
            )
        )


#: Known issues.
KNOWN_ISSUES = (
    KnownIssue(
        "stats-2020",
        "clinvar",
        "Clinvar",
        re.escape("No count for 3 chromosomes: ['chrM', 'chrX', 'chrY']"),
    ),
    KnownIssue(
        "stats-2020",
        "gnomAD_exomes",
        "GnomadExomes",
        re.escape("No count for 1 chromosomes: ['chrY']"),
    ),
    KnownIssue(
        "stats-2020",
        "gnomAD_exomes",
        "GnomadExomes",
        re.escape("No count for 1 chromosomes: ['chrY']"),
    ),
    KnownIssue(
        "stats-2020",
        "ensembl_regulatory",
        "EnsemblRegulatoryFeature",
        re.escape("No count for 1 chromosomes: ['chrY']"),
    ),
)


def report(args: ReportArgs):
    """Run report generation."""
    logger.info("Running extraction with args\n\n%s", json.dumps(vars(args), indent=2))
    path_report = Path(args.path_report)
    paths_stats = list(map(Path, args.paths_stats))

    logger.info("Loading Reports")
    reports: typing.Dict[ReportKey, typing.Dict] = {}
    for path_stats in paths_stats:
        logger.info("considering %s", path_stats)
        for path_json in path_stats.glob("*.json"):
            key = ReportKey(path_stats.name, *path_json.name.split(".")[:3])
            logger.info("  loading %s", path_json)
            with path_json.open("rt") as inputf:
                reports[key] = json.load(inputf)

    logger.info("Running sanity checks")
    results = []
    for key, report in reports.items():
        for check in SANITY_CHECKS:
            if key.matches(check.report_key):
                results += check.check(key, report, args)

    logger.info("Here is your report")
    fn = {
        ReportLevel.INFO: logger.info,
        ReportLevel.WARNING: logger.warning,
        ReportLevel.ERROR: logger.error,
    }
    counts = {k: 0 for k in ReportLevel}
    known = 0
    for msg in results:
        txt = (
            "%s | %s | %s | %s | %s",
            msg.level.name,
            msg.stats,
            msg.table_group,
            msg.table,
            msg.msg,
        )
        for k in KNOWN_ISSUES:
            if k.matches(msg):
                known += 1
                logger.info("KNOWN " + txt[0], *txt[1:])
                break
        else:
            counts[msg.level] += 1
            fn[msg.level](*txt)
    logger.info("# ___ SUMMARY ___")
    logger.info("# known: %d", known)
    for l in ReportLevel:
        logger.info("# %s: %d", l.name, counts[l])


def main(argv=None):
    parser = argparse.ArgumentParser(description="Sanity checks for varfish-db-downloader results")
    subparsers = parser.add_subparsers(help="sub command help")

    parser_extract = subparsers.add_parser("extract", help="Extract statistics from a TSV file.")
    parser_extract.set_defaults(func=extract, args_class=ExtractArgs)
    parser_extract.add_argument(dest="path_stats", type=str)
    parser_extract.add_argument(dest="path_import_versions", type=str)
    parser_extract.add_argument("--line-limit", type=int)
    parser_extract.add_argument("--only", type=str)

    parser_report = subparsers.add_parser("report", help="Generate report based on statistics TSV.")
    parser_report.set_defaults(func=report, args_class=ReportArgs)
    parser_report.add_argument(dest="path_report", type=str)
    parser_report.add_argument(dest="paths_stats", type=str, nargs="+")

    args = parser.parse_args(argv)

    if hasattr(args, "only"):
        if args.only:
            args.only = tuple(args.only.split(","))
        else:
            args.only = ()

    if not hasattr(args, "func"):
        parser.print_help()
    else:
        return args.func(cattr.structure(vars(args), args.args_class))


if __name__ == "__main__":
    sys.exit(int(main() or 0))
