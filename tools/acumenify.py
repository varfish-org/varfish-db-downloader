#!/usr/bin/env python3
"""Sanity checking for releases as built by ``varfish-db-downloader``.

This tool takes as the input a release built by ``varfish-db-downloader`` which is identified as
a path to the ``import_versions.tsv`` file that is generated in the release directory.  For each
table in the release it will create a JSON file in an outupt directory that contains summary
statistics.


-----
Usage
-----

First point the ``extract`` command at a release directory with ``import_versions.tsv`` file.

::

    python acumenify.py extract \\
        reports/20210728-GRCh37 \\
        releases/20210728/grch37/import_versions.tsv

This will generate report JSON files in the folder ``reports/20210728-GRCh37``.

You can then run sanity checks on one or multiple report folders and also see a tabular report on
the number of lines per chromosome in each file by running the ``report`` command.

::

    python tools/acumenify.py report report.xlsx reports/*

Note that this will also print the report parts to the stderr but large tables will be abbreviated.
Essentially this converts pandas data frames to strings which will leave out columns and rows as
not to overload the user's terminal.
"""

import argparse
import collections
import enum
import functools
import json
import multiprocessing
from pathlib import Path
import re
import sys
import time
import typing

import attr
import pandas as pd
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
    """Arguments for the ``extract`` command."""

    #: Path to ``import_versions.tsv`` file.
    path_import_versions: str
    #: Path to write statistics files to.
    path_stats: str
    #: Degree of parallelism
    processes: int = 0
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
    """Arguments for the ``report`` command."""

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
    """A record from ``import_versions.tsv`` files."""

    #: The genome build, e.g., ``"GRCh37"``.
    build: str
    #: The table group name.
    table_group: str
    #: The table name.
    version: str

    def get_key(self) -> str:
        """Return separated key for import version"""
        return "%s.%s" % (self.build, self.table_group)


@attr.s(frozen=True, auto_attribs=True)
class StatsKey:
    """The key for a statics entry."""

    #: The table group name.
    table_group: str
    #: The table name.
    table: str
    #: The version of data that the statistics were created for.
    chrom: str


class ColumnType(enum.Enum):
    """Enumeration for colum types."""

    #: Type is not known yet.
    UNKNOWN = 1
    #: Type is an enumerate (String-typed with relatively few values).
    ENUM = 2
    #: Arbitrary string.
    STRING = 3
    #: Integer.
    INT = 4
    #: Floating point number.
    FLOAT = 5
    #: DNA String.
    DNA = 6

    def with_generalized(self, other: "ColumnType") -> "ColumnType":
        """Return a ``Column`` that is the generalizatio of ``self`` and ``other``.

        The generalization is a type that is capable of representing values of both ``self`` and
        ``other``.  Note that ``UNKNOWN`` is seen as a wild card and replaced by the other
        respective type.
        """

        return {
            (ColumnType.UNKNOWN, ColumnType.UNKNOWN): ColumnType.UNKNOWN,
            (ColumnType.UNKNOWN, ColumnType.ENUM): ColumnType.ENUM,
            (ColumnType.UNKNOWN, ColumnType.STRING): ColumnType.STRING,
            (ColumnType.UNKNOWN, ColumnType.INT): ColumnType.INT,
            (ColumnType.UNKNOWN, ColumnType.FLOAT): ColumnType.FLOAT,
            (ColumnType.UNKNOWN, ColumnType.DNA): ColumnType.DNA,
            (ColumnType.ENUM, ColumnType.UNKNOWN): ColumnType.ENUM,
            (ColumnType.ENUM, ColumnType.ENUM): ColumnType.ENUM,
            (ColumnType.ENUM, ColumnType.STRING): ColumnType.STRING,
            (ColumnType.ENUM, ColumnType.INT): ColumnType.STRING,
            (ColumnType.ENUM, ColumnType.FLOAT): ColumnType.STRING,
            (ColumnType.ENUM, ColumnType.DNA): ColumnType.STRING,
            (ColumnType.STRING, ColumnType.UNKNOWN): ColumnType.STRING,
            (ColumnType.STRING, ColumnType.ENUM): ColumnType.STRING,
            (ColumnType.STRING, ColumnType.STRING): ColumnType.STRING,
            (ColumnType.STRING, ColumnType.INT): ColumnType.STRING,
            (ColumnType.STRING, ColumnType.FLOAT): ColumnType.STRING,
            (ColumnType.STRING, ColumnType.DNA): ColumnType.STRING,
            (ColumnType.INT, ColumnType.UNKNOWN): ColumnType.INT,
            (ColumnType.INT, ColumnType.ENUM): ColumnType.STRING,
            (ColumnType.INT, ColumnType.STRING): ColumnType.STRING,
            (ColumnType.INT, ColumnType.INT): ColumnType.INT,
            (ColumnType.INT, ColumnType.FLOAT): ColumnType.STRING,
            (ColumnType.INT, ColumnType.DNA): ColumnType.STRING,
            (ColumnType.FLOAT, ColumnType.UNKNOWN): ColumnType.FLOAT,
            (ColumnType.FLOAT, ColumnType.ENUM): ColumnType.STRING,
            (ColumnType.FLOAT, ColumnType.STRING): ColumnType.STRING,
            (ColumnType.FLOAT, ColumnType.INT): ColumnType.FLOAT,
            (ColumnType.FLOAT, ColumnType.FLOAT): ColumnType.FLOAT,
            (ColumnType.FLOAT, ColumnType.DNA): ColumnType.STRING,
            (ColumnType.DNA, ColumnType.UNKNOWN): ColumnType.DNA,
            (ColumnType.DNA, ColumnType.ENUM): ColumnType.STRING,
            (ColumnType.DNA, ColumnType.STRING): ColumnType.STRING,
            (ColumnType.DNA, ColumnType.INT): ColumnType.STRING,
            (ColumnType.DNA, ColumnType.FLOAT): ColumnType.STRING,
            (ColumnType.DNA, ColumnType.DNA): ColumnType.DNA,
        }[(self, other)]


#: Header known to hold chromosomes.
HEADER_CHROM = "chromosome"

#: Regular expression for a DNA string.
RE_DNA = r"^[ACGTNacgtn]+$"


@attr.s(frozen=True, auto_attribs=True)
class ChromHarmonizer:
    """Helper class that harmonizes chromosome name.

    ..code-block:: python

        >>> h = ChromHarmonizer("chr", "M")
        >>> h.apply("chrX")
        "chrX"
        >>> h.apply("X")
        "chrX"
        >>> h.apply("M")
        "chrM"
        >>> h.apply("MT")
        "chrM"
        >>> h = ChromHarmonizer("", "MT")
        >>> h.apply("chrX")
        "X"
        >>> h.apply("X")
        "X"
        >>> h.apply("M")
        "MT"
        >>> h.apply("MT")
        "MT"
    """

    #: Prefix to apply to chromosomes, either ``"chr"`` or ``""``.
    harmonize_chroms: str
    #: The name to normalize the mitochondrial chromosome to, e.g., ``"M"``.
    harmonize_chrmt: str

    def apply(self, value: str) -> str:
        """Apply harmonization to chromosome name ``value``."""
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


def _identity(x):
    """Identity function."""
    return x


class ColAggregator:
    """Helper class for aggregating column values.
    
    The aggregator will consider at most ``guess_len`` values and then guess the values based on
    this.  The distinction between enumerations and strings is made based the value set
    cardinality being not above ``max_enum_size``.  Guessing can be forced by calling ``finished()``.
    """

    def __init__(
        self,
        name,
        chrom_harmonizer: ChromHarmonizer,
        guess_len: int = 10000,
        max_enum_size: int = 100,
    ):
        #: Column name
        self.name = name
        #: Chromozome harmoniser.
        self.chrom_harmonizer = chrom_harmonizer
        if self.name == HEADER_CHROM:
            self.fn_harmonize = chrom_harmonizer.apply
        else:
            self.fn_harmonize = _identity
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

    def with_added(self, other: "ColAggregator") -> "ColAggregator":
        """Return sum of this and other col aggregator.
        
        That is, the count values will be added and the set values will be merged.
        """
        assert self.name == other.name
        assert self.guess_len == other.guess_len
        assert self.max_enum_size == other.max_enum_size
        res = ColAggregator(self.name, self.chrom_harmonizer, self.guess_len, self.max_enum_size)
        res.type_ = self.type_.with_generalized(other.type_)
        if self.values is None and other.values is None:
            res.values = None
        elif self.values is None:
            res.values = dict(other.values)
        elif other.values is None:
            res.values = dict(self.values)
        else:
            res.values = {}
            for k in self.values.keys() | other.values.keys():
                res.values[k] = self.values.get(k, 0) + other.values.get(k, 0)
        res.counter = self.counter + other.counter
        return res

    def process(self, value: str) -> None:
        """Process one value occuring in a column's cell."""
        self.counter += 1
        if self.counter >= self.guess_len:
            self.finish()
        value = self.fn_harmonize(value)
        if self.values is not None:
            self.values.setdefault(value, 0)
            self.values[value] += 1

    def finish(self) -> None:
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

    def to_dict(self) -> typing.Dict:
        """Return dict with results."""
        result = {"type": self.type_.name}
        if self.type_ == ColumnType.ENUM:
            result["values"] = list(sorted(self.values))
        return result


class Aggregator:
    """Helper class for aggregating values column-wise, in multiple columns."""

    def __init__(self, header, chrom_harmonizer):
        #: Header.
        self.header = header
        #: Chromozome harmoniser.
        self.chrom_harmonizer = chrom_harmonizer
        #: Count by chromosome, if any, else None.
        self.by_chrom = {}
        #: Aggregate statistics for each column.
        self.by_column = {col: ColAggregator(col, chrom_harmonizer) for col in header}

    def with_added(self, other):
        """Return ``Aggregator`` object with added values.

        That is, the count values will be added and the set values will be merged.
        """
        assert self.header == other.header
        assert self.by_column.keys() == other.by_column.keys()
        res = Aggregator(self.header, self.chrom_harmonizer)
        res.by_chrom = dict(self.by_chrom)
        for k, v in other.by_chrom.items():
            res.by_chrom[k] = res.by_chrom.get(k, 0) + v
        for k in self.by_column.keys():
            v1 = self.by_column[k]
            v1.finish()
            v2 = other.by_column[k]
            v2.finish()
            res.by_column[k] = v1.with_added(v2)
        return res

    def process(self, record: typing.Dict) -> None:
        """Process a record.

        A record is a ``dict`` that maps header names to column/cell values.
        """
        chrom = self.chrom_harmonizer.apply(record.get(HEADER_CHROM, "."))
        self.by_chrom.setdefault(chrom, 0)
        self.by_chrom[chrom] += 1
        for k, v in record.items():
            self.by_column[k].process(v)

    def finish(self):
        """Force guessing of values if not done so far and compute sum of values for all
        chromosomes with chromosome key ``"."``.
        """
        for a in self.by_column.values():
            a.finish()
        if "." not in self.by_chrom:
            self.by_chrom["."] = sum(self.by_chrom.values())


def do_extraction_job(path: Path, args: ExtractArgs) -> Aggregator:
    """Run one aggregation job for the file at the given ``path``.

    Will display progress using the ``tqdm`` library if only one thread is used and otherwise
    will display progress to logging.
    """
    disable_tqdm = args.processes > 1
    logger.info(
        "  processing file %s%s",
        path,
        " (at most %s lines)" % args.line_limit if args.line_limit else "",
    )

    agg = None
    header = None

    with path.open("rt") as inputf:
        prev = time.time()
        modulo = 200_000  # print progress every so many lines
        for lineno, line in tqdm(enumerate(inputf), unit="rec", disable=disable_tqdm):
            if disable_tqdm and lineno and lineno % modulo == 0:  # progress to logs
                curr = time.time()
                lines_per_sec = modulo / (curr - prev)
                prev = curr
                logger.info(
                    "    ... processing %s ... at line %s (%.2f per sec)",
                    path.name,
                    "{:,}".format(lineno),
                    lines_per_sec,
                )
            if args.line_limit and lineno > args.line_limit:
                logger.debug("    ~> stopping at line limit of %s", args.line_limit)
                break

            while line and line[-1] in "\r\n":  # trim line endings but no other space
                line = line[:-1]

            arr = line.split("\t")
            if not header:  # first line, initialize header and aggregator
                header = arr
                agg = Aggregator(
                    header, ChromHarmonizer(args.harmonize_chroms, args.harmonize_chrmt)
                )
            else:  # every line after the first, process record
                if len(header) != len(arr):
                    raise ExtractionError(
                        "Record %s has %d fields but header had %s"
                        % (lineno, len(arr), len(header))
                    )
                agg.process(dict(zip(header, arr)))

    # Finish up.
    agg.finish()
    logger.info("    => done with %s", path.name)
    return agg


def do_extraction(path_base: Path, record: ImportVersion, args: ExtractArgs) -> typing.Dict:
    """Perform extraction for all tables belonging to the table group defined in ``record``.

    The payload data is assumed to be relative to ``path_base``, parametrize based on ``args``.
    """
    logger.info("Looking into %s", record)

    # Collect files for each table.
    by_table = {}
    path_dir = path_base / record.build / record.table_group / record.version
    logger.debug("  path = %s", path_dir)
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
        logger.info("processing table %s in %s", table, record.table_group)
        if args.processes >= 1:
            pool = multiprocessing.Pool(args.processes)
            aggs = list(
                pool.map(
                    functools.partial(do_extraction_job, args=args), sorted(set(paths)), chunksize=1
                )
            )
        else:
            # args.processes <= 0, do not use multiprocessing.
            aggs = [do_extraction_job(path, args) for path in sorted(set(paths))]
        agg = aggs[0]
        for a in aggs[1:]:
            agg = agg.with_added(a)
        result[table] = {
            "$schema_version": "0.1.0",
            "import_version": vars(record),
            "by_chrom": agg.by_chrom,
            "by_col": {k: v.to_dict() for k, v in agg.by_column.items()},
        }

    return result


def extract(args: ExtractArgs):
    """Run extraction command based on configuration in ``args``."""
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
            if not header:  # first line, read header
                expected = ["build", "table_group", "version"]
                header = arr
                if header != expected:
                    logger.error("Header unexpected, WAS: %s, expected: %s", header, expected)
            else:  # any other line, start extraction for table group
                record = ImportVersion(**dict(zip(header, arr)))
                if record in stats:  # ignore duplicates
                    logger.debug("Ignoring second occurence of %s", record)
                elif args.only and record.table_group not in args.only:  # allow filter on group
                    logger.debug(
                        "Table group %s is not in --only=%s", record.table_group, args.only
                    )
                else:  # the actual handling of the extraction
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
    """Enumeration for report levels."""

    #: Informative message in report.
    INFO = 1
    #: Warning message in report.
    WARNING = 2
    #: Error message in report.
    ERROR = 3


@attr.s(frozen=True, auto_attribs=True)
class ReportMessage:
    """An overall message in the output report."""

    #: The level of the message.
    level: ReportLevel
    #: The name of the input folder.
    stats: str
    #: The table group.
    table_group: str
    #: The table.
    table: str
    #: The message text.
    msg: str

    def is_known(self, known_issues: typing.List["KnownIssue"]) -> bool:
        """Return whether ``self`` is known given the list of ``known_issues``."""
        return any(k.matches(self) for k in known_issues)


@attr.s(frozen=True, auto_attribs=True)
class ReportKey:
    """Key for an entry in a report."""

    #: The input folder.
    stats: typing.Optional[str] = None
    #: The data release.
    release: typing.Optional[str] = None
    #: The Table group.
    table_group: typing.Optional[str] = None
    #: The table.
    table: typing.Optional[str] = None

    def matches(self, other: "ReportKey") -> bool:
        """Whether or not ``self`` matches ``other`` with ``None`` representing wildcars."""

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
    """Wrapper container for sanity check function."""

    #: The sanity check function.
    check: typing.Callable[
        [ReportKey, typing.Dict, ReportArgs, typing.Dict], typing.List[ReportMessage]
    ]


def sc_report_allchroms(
    key: ReportKey, report: typing.Dict, args: ReportArgs, all_reports: typing.Dict
) -> typing.List[ReportMessage]:
    """Sanity check that tests that there is data for all chromosomes."""
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


def sc_report_noextrachroms(
    key: ReportKey, report: typing.Dict, args: ReportArgs, all_reports: typing.Dict
) -> typing.List[ReportMessage]:
    """Sanity check that tests that there are only chr1..chr22, chrY, chrX, chrMT."""
    h = ChromHarmonizer(args.harmonize_chroms, args.harmonize_chrmt).apply
    result = []

    chroms = set(map(h, CHROMS)) | {"."}
    extra = report["by_chrom"].keys() - chroms
    if extra:
        result.append(
            ReportMessage(
                ReportLevel.WARNING,
                key.stats,
                key.table_group,
                key.table,
                "Counts for %d unexpected chromosomes: %s" % (len(extra), list(sorted(extra))),
            )
        )

    return result


def sc_global_missingtable(
    key: ReportKey, report: typing.Dict, args: ReportArgs, all_reports: typing.Dict
) -> typing.List[ReportMessage]:
    """A "global" sanity check that tests that checks for missing tables."""
    result = []

    all_stats = set()

    by_table = {}
    for key in all_reports.keys():
        all_stats.add(key.stats)
        by_table.setdefault((key.table_group, key.table), set()).add(key.stats)

    for (table_group, table), stats in by_table.items():
        if stats != all_stats:
            for s in all_stats - stats:
                result.append(
                    ReportMessage(
                        ReportLevel.ERROR,
                        s,
                        table_group,
                        table,
                        "Missing table %s/%s for dataset %s" % (table_group, table, s),
                    )
                )

    return result


#: Define the sanity checks to run per report.
SANITY_CHECKS_PER_REPORT = (SanityCheck(sc_report_allchroms), SanityCheck(sc_report_noextrachroms))
#: Define sanity checks to be run once globally.
SANITY_CHECKS_GLOBAL = (SanityCheck(sc_global_missingtable),)


@attr.s(frozen=True, auto_attribs=True)
class KnownIssue:
    """A known issue."""

    #: Regular expression for dataset name.
    dataset: typing.Optional[str]
    #: Regular expression for table group.
    table_group: str
    #: Regular expression for table name.
    table: str
    #: Regular expression for message.
    message: str
    #: A comment for the known issue.
    comment: typing.Optional[str] = None

    def matches(self, msg: ReportMessage):
        return all(
            (
                not self.dataset or re.match(self.dataset, msg.stats),
                not self.table_group or re.match(self.table_group, msg.table_group),
                not self.table or re.match(self.table, msg.table),
                not self.message or re.match(self.message, msg.msg),
            )
        )


#: Known issues.
KNOWN_ISSUES = (
    KnownIssue(
        "20201006-GRCh37",
        "clinvar",
        "Clinvar",
        re.escape("No count for 3 chromosomes: ['chrM', 'chrX', 'chrY']"),
        "The first data release only has autosomal chromosomes.",
    ),
    KnownIssue(
        "20201006-GRCh37",
        "gnomAD_exomes",
        "GnomadExomes",
        re.escape("No count for 1 chromosomes: ['chrY']"),
        "The first data release has an old gnomAD version without chrY.",
    ),
    KnownIssue(
        "20201006-GRCh37",
        "ensembl_regulatory",
        "EnsemblRegulatoryFeature",
        re.escape("No count for 1 chromosomes: ['chrY']"),
        "The first data release was accidentally missing chrY for ensembl regulatory features.",
    ),
    KnownIssue(".*-GRCh38", "ExAC", "Exac", None, "ExAC is not available for GRCh38."),
    KnownIssue(".*-GRCh38", "gnomAD_SV", "GnomAdSv", None, "ExAC is not available for GRCh38."),
    KnownIssue(
        ".*-GRCh38",
        r"tads_imr90|tads_hesc",
        None,
        None,
        "These TAD files are not available for GRCh38.",
    ),
    KnownIssue(
        ".*-GRCh38",
        r"thousand_genomes",
        None,
        None,
        "Thousand genomes data not availablef or GRCh38.",
    ),
    KnownIssue(".*-GRCh38", r"vista", None, None, "VISTA is not available for GRCh38."),
)


def report(args: ReportArgs) -> none:
    """Run report generation."""

    def normalize(s):
        """Normalize table (group) name."""
        return s.replace("-", "_")

    logger.info("Running extraction with args\n\n%s", json.dumps(vars(args), indent=2))
    path_report = Path(args.path_report)
    paths_stats = list(map(Path, args.paths_stats))

    logger.info("Loading Reports")
    reports: typing.Dict[ReportKey, typing.Dict] = {}
    for path_stats in paths_stats:
        logger.debug("considering %s", path_stats)
        for path_json in path_stats.glob("*.json"):
            key = ReportKey(path_stats.name, *map(normalize, path_json.name.split(".")[:3]))
            logger.debug("  loading %s", path_json)
            with path_json.open("rt") as inputf:
                reports[key] = json.load(inputf)

    # Build per contig report.
    per_contig_dicts = []
    for key, report in reports.items():
        if report["$schema_version"] != "0.1.0":
            raise RuntimError("Cannot handle schema version %s" % report["$schema_version"])
        for contig, count in report["by_chrom"].items():
            if len(report["by_chrom"]) == 1 or contig != ".":
                per_contig_dicts.append(
                    {
                        "dataset": key.stats,
                        "table_group": normalize(report["import_version"]["table_group"]),
                        "table": normalize(key.table),
                        "version": report["import_version"]["version"],
                        "contig": contig,
                        "count": count,
                    }
                )
    df_per_contig = pd.DataFrame(per_contig_dicts)
    logger.info("Per contig results:\n\n%s\n", df_per_contig)

    df_pivoted = pd.pivot_table(
        data=df_per_contig,
        index=["table_group", "table", "contig"],
        columns=["dataset"],
        values=["count"],
        aggfunc="median",
        margins=True,
        margins_name="median",
    )
    logger.info("Pivoted results:\n\n%s\n", df_pivoted)

    # Build overall sanity check results.
    logger.info("Running sanity checks")
    results = []
    for key, report in reports.items():
        for check in SANITY_CHECKS_PER_REPORT:
            results += check.check(key, report, args, reports)
    for check in SANITY_CHECKS_GLOBAL:
        results += check.check(key, None, args, reports)

    df_result = pd.DataFrame(
        [
            {
                "known_issue": msg.is_known(KNOWN_ISSUES),
                "level": msg.level.name,
                "stats": msg.stats,
                "table_group": msg.table_group,
                "table": msg.table,
                "msg": msg.msg,
            }
            for msg in results
        ]
    )
    logger.info("Overall sanity check messages:\n\n%s\n", df_result)

    summary_counts = {
        "known": 0,
        ReportLevel.INFO.name: 0,
        ReportLevel.WARNING.name: 0,
        ReportLevel.ERROR.name: 0,
    }
    for msg in results:
        if msg.is_known(KNOWN_ISSUES):
            summary_counts["known"] += 1
        else:
            summary_counts[msg.level.name] += 1

    df_summary = pd.DataFrame([{"level": k, "count": v} for k, v in summary_counts.items()])
    logger.info("Overall sanity check stats:\n\n%s\n", df_summary)

    logger.info("Writing output file %s", args.path_report)
    with pd.ExcelWriter(args.path_report) as writer:
        df_summary.to_excel(writer, sheet_name="Overall Stats")
        df_result.to_excel(writer, sheet_name="Sanity Checks")
        df_per_contig.to_excel(writer, sheet_name="Counts")
        df_pivoted.to_excel(writer, sheet_name="Pivoted Counts")

    logger.info("All done. Have a nice day!")


def main(argv=None):
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Sanity checks for varfish-db-downloader results")
    subparsers = parser.add_subparsers(help="sub command help")

    parser_extract = subparsers.add_parser("extract", help="Extract statistics from a TSV file.")
    parser_extract.set_defaults(func=extract, args_class=ExtractArgs)
    parser_extract.add_argument(dest="path_stats", type=str)
    parser_extract.add_argument(dest="path_import_versions", type=str)
    parser_extract.add_argument("--line-limit", type=int)
    parser_extract.add_argument("--only", type=str)
    parser_extract.add_argument("--processes", type=int, default=0, help="Number of processes")

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
