#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import logging
import re
import sys
import typing

import attrs
import cattrs
import vcfpy

@attrs.frozen()
class HgncEntry:
    hgnc_id: str
    ensembl_gene_id: str


@attrs.frozen()
class AlignmentMeta:
    hgnc_id: str
    enst_id: str
    species: str
    exon_idx: int
    exon_count: int
    exon_len: int
    in_frame: int
    out_frame: int
    location: typing.Tuple[str, int, int, str]


@attrs.frozen()
class AlignmentBlock:
    meta: AlignmentMeta
    species: typing.List[str]
    aa_seqs: typing.List[str]


def setup_logging(args):
    """Setup logger."""
    logging.basicConfig(
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
    )
    logger = logging.getLogger("")
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


def read_contigs(fai_path):
    """Read contig information (name, length) from FAI file at ``fai_path``."""
    with open(fai_path, "rt") as fai:
        for line in fai:
            arr = line.split("\t")
            yield arr[0], int(arr[1])


def build_header(contigs, species):
    contig_map = {contig[0]: contig[1] for contig in contigs}

    header = vcfpy.Header()
    header.samples = vcfpy.SamplesInfos([])
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    seen = set()
    for name, length in contigs:
        header.add_contig_line({"ID": name, "length": length})
        seen.add(name)
    for a, b in (("M", "MT"), ("MT", "M"), ("chrM", "chrMT"), ("chrMT", "chrM")):
        if a in seen and not b in seen:
            header.add_contig_line({"ID": b, "length": contig_map[a]})
    header.add_line(vcfpy.HeaderLine("species", ",".join(species)))
    header.add_info_line(
        {
            "ID": "END",
            "Description": "End position of the alignment",
            "Type": "Integer",
            "Number": 1,
        }
    )
    header.add_info_line(
        {"ID": "HGNC_ID", "Description": "HGNC ID of gene", "Type": "String", "Number": 1}
    )
    header.add_info_line(
        {"ID": "ENST_ID", "Description": "ENSEMBL transcript ID", "Type": "String", "Number": 1}
    )
    header.add_info_line(
        {"ID": "EXON", "Description": "Index of exon in transcript", "Type": "Integer", "Number": 1}
    )
    header.add_info_line(
        {
            "ID": "EXON_COUNT",
            "Description": "Number of exons in transcript",
            "Type": "Integer",
            "Number": 1,
        }
    )
    header.add_info_line(
        {
            "ID": "ALIGNMENT",
            "Description": "Amino acid alignment at this location",
            "Type": "String",
            "Number": 1,
        }
    )
    return header


def fasta_header_to_meta(line, have_chr, enst_to_hgnc_id, enst_to_real):
    arr = line[1:].split(" ")
    if len(arr) == 5:
        exon, exon_len, in_frame, out_frame, location = arr
    else:
        (exon, exon_len, in_frame, out_frame), location = arr, ""
    enst_id, species, exon_no, exon_count = exon.split("_")
    enst_nover = enst_id.split(".")[0]
    if enst_nover not in enst_to_hgnc_id:
        logging.debug("ENST %s not found in HGNC ID map", enst_nover)
        return None

    hgnc_id = enst_to_hgnc_id.get(enst_nover, "")
    if location:
        location = location.split(";")[0]
        location, strand = location[:-1], location[-1]
        chrom, range_ = location.split(":")
        if have_chr and not chrom.startswith("chr"):
            chrom = "chr%s" % chrom
        elif not have_chr and chrom.startswith("chr"):
            chrom = chrom[3:]
        start, end = range_.split("-")
        location = (chrom, int(start) - 1, int(end), strand)
    else:
        location = None
    return AlignmentMeta(
        hgnc_id=hgnc_id,
        enst_id=enst_to_real[enst_id.split('.')[0]],
        species=species,
        exon_idx=int(exon_no),
        exon_count=int(exon_count),
        exon_len=int(exon_len),
        in_frame=int(in_frame),
        out_frame=int(out_frame),
        location=location,
    )


def build_block(lines, have_chr, enst_to_hgnc_id, enst_to_real):
    """Build ``AlignmentBlock`` from FASTA lines."""
    hg_meta = fasta_header_to_meta(lines[0], have_chr, enst_to_hgnc_id, enst_to_real)
    if not hg_meta:
        return None
    species = []
    aa_seqs = []
    assert len(lines) % 2 == 0
    for i in range(0, len(lines), 2):
        header, aas = lines[i : i + 2]
        meta = fasta_header_to_meta(header, have_chr, enst_to_hgnc_id, enst_to_real)
        species.append(meta.species)
        aa_seqs.append(aas)
    return AlignmentBlock(meta=hg_meta, species=species, aa_seqs=aa_seqs)


def read_blocks(fa_gz, have_chr, enst_to_hgnc_id, enst_to_real):
    """Read and yield ``AlignmentBlock`` elements from ``fa_gz``."""
    buf = []
    for line in fa_gz:
        line = line.strip()
        if not line:
            if buf:
                val = build_block(buf, have_chr, enst_to_hgnc_id, enst_to_real)
                if val:
                    yield val
                buf = []
        else:
            buf.append(line)
    if buf:
        val = build_block(buf, have_chr, enst_to_hgnc_id)
        if val:
            yield val


def pos_magic(exon_location, rel_start, rel_end):
    _chrom, exon_start, exon_end, strand = exon_location
    assert strand in "+-"
    if strand == "+":
        return exon_start + rel_start, exon_start + rel_end
    else:
        end = exon_end - rel_start
        start = exon_end - rel_end
        return start, end


def block_to_records(block, prev_block):
    """Given an alignment block, yield the VCF records.

    NB: the first/last amino acids are stored for the exon with the major part of the codon.
    """
    logging.debug("Starting new block")
    meta = block.meta
    location = block.meta.location
    in_frame = block.meta.in_frame
    out_frame = block.meta.out_frame

    logging.debug("prev_block: %s", prev_block)
    logging.debug("block: %s", block)
    logging.debug("%s\n%s", meta, "\n".join(block.aa_seqs))

    # Special case handling for first codon.
    assert not in_frame or prev_block
    if in_frame == 1:
        # The alignment for this exon has the amino acid, yield last amino acid for previous exon.
        prev_location = prev_block.meta.location
        start, end = pos_magic(
            prev_location,
            prev_location[2] - prev_location[1] - 1,
            prev_location[2] - prev_location[1],
        )
        logging.debug("in_frame == 1, start, end = %d, %d", start, end)
        yield vcfpy.Record(
            CHROM=prev_block.meta.location[0],
            POS=start + 1,
            ID=[],
            REF="N",
            ALT=[],
            FILTER=[],
            QUAL=None,
            INFO={
                "END": end,
                "HGNC_ID": prev_block.meta.hgnc_id,
                "ENST_ID": prev_block.meta.enst_id,
                "EXON": prev_block.meta.exon_idx,
                "EXON_COUNT": prev_block.meta.exon_count,
                "ALIGNMENT": "".join(seq[0] for seq in block.aa_seqs),
            },
            FORMAT={},
            calls=[],
        )
        # We start at the first codon of this exon
        starts = [0]
        ends = [2]
    elif in_frame == 2:
        # The alignment for the previous exon has the amino acid, yield first amino acid for this exon.
        start, end = pos_magic(location, 0, 1)
        logging.debug("in_frame == 2, start, end = %d, %d", start, end)
        yield vcfpy.Record(
            CHROM=location[0],
            POS=start + 1,
            ID=[],
            REF="N",
            ALT=[],
            FILTER=[],
            QUAL=None,
            INFO={
                "END": end,
                "HGNC_ID": meta.hgnc_id,
                "ENST_ID": prev_block.meta.enst_id,
                "EXON": meta.exon_idx,
                "EXON_COUNT": meta.exon_count,
                "ALIGNMENT": "".join(seq[-1] for seq in prev_block.aa_seqs),
            },
            FORMAT={},
            calls=[],
        )
        # We start at the second codon on this exon
        starts = [1]
        ends = [4]
    else:
        # Start at codon border
        logging.debug("in_frame == 0, start, end = 0, 3")
        starts = [0]
        ends = [3]

    # Handle major part of exon.
    nts = location[2] - location[1]
    assert (nts - (3 - in_frame) - out_frame) % 3 == 0
    starts = starts + list(range(ends[0], nts, 3))
    ends = ends + list(range(ends[0] + 3, nts, 3))
    if ends[-1] != nts:
        ends.append(nts)
    logging.debug("nts=%s", nts)
    logging.debug("starts=%s", starts)
    logging.debug("ends=%s", ends)

    if out_frame == 2:
        # We have the amino acid for the last partial codon.  If out_frame == 1 then the next exon will take care of writing record.
        starts += [ends[-1]]
        ends += [ends[-1] + 2]
    for i, (start, end) in enumerate(zip(starts, ends)):
        if i >= meta.exon_len:
            logging.debug("Too short AA seq found for %s, this happens...", meta.hgnc_id)
            continue
        start, end = pos_magic(location, start, end)
        yield vcfpy.Record(
            CHROM=location[0],
            POS=start + 1,
            ID=[],
            REF="N",
            ALT=[],
            FILTER=[],
            QUAL=None,
            INFO={
                "END": end,
                "HGNC_ID": meta.hgnc_id,
                "ENST_ID": meta.enst_id,
                "EXON": meta.exon_idx,
                "EXON_COUNT": meta.exon_count,
                "ALIGNMENT": "".join(seq[i] for seq in block.aa_seqs),
            },
            FORMAT={},
            calls=[],
        )


def run(args):
    setup_logging(args)
    logging.info("Starting processing")
    logging.info("Arguments are %s", args)

    logging.info("Loading ENST to ENSG map")
    enst_to_real = {}  # GRCh37 has ucsc ids for ENST, so we hack it
    ensg_to_ensts = {}
    with open(args.enst_ensg[0], "rt") as inputf:
        reader = csv.DictReader(inputf, delimiter="\t")
        for row in reader:
            ensg_to_ensts.setdefault(row["ensg"], []).append(row["enst"])
            enst_to_real[row["enst"]] = row.get("real_enst", row["enst"])

    logging.info("Loading HGNC info")
    enst_to_hgnc_id = {}
    with open(args.hgnc_jsonl[0], "rt") as inputf:
        for line in inputf:
            record = json.loads(line)
            if not "ensembl_gene_id" in record:
                continue
            entry: HgncEntry = cattrs.structure(record, HgncEntry)
            for enst in ensg_to_ensts.get(entry.ensembl_gene_id, []):
                enst_to_hgnc_id[enst] = entry.hgnc_id

    logging.info("Reading contigs")
    contigs = list(read_contigs(args.reference[0] + ".fai"))
    have_chr = contigs[0][0].startswith("chr")

    logging.info("Create VCF header and open output file")
    with gzip.open(args.input[0], "rt") as fa_gz:
        for first_block in read_blocks(fa_gz, have_chr, enst_to_hgnc_id, enst_to_real):
            break
    header = build_header(contigs, first_block.species)
    prev_contig = None
    if args.output:
        output = vcfpy.Writer.from_path(args.output, header)
    else:
        output = vcfpy.Writer.from_stream(sys.stdout, header)
    with output as writer:
        with gzip.open(args.input[0], "rt") as fa_gz:
            prev_block = None
            for block in read_blocks(fa_gz, have_chr, enst_to_hgnc_id, enst_to_real):
                records = list(block_to_records(block, prev_block))
                if re.match(args.contig_regexp, records[0].CHROM):
                    if records[0].CHROM != prev_contig:
                        logging.info("Now on contig %s", records[0].CHROM)
                        prev_contig = records[0].CHROM
                    for record in records:
                        writer.write_record(record)
                else:
                    if prev_block != block:
                        logging.info("Will skip records on contig %s", records[0].CHROM)
                prev_block = block

    logging.info("All done. Have a nice day!")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert UCSC amino acid multiz alignment to VCF file"
    )
    parser.add_argument(
        "--contig-regexp",
        type=str,
        default="^(chr)?[1-9XYM][0-9T]?$",
        help="Regular expression to filter contig by",
    )
    parser.add_argument(
        "hgnc_jsonl",
        type=str,
        nargs=1,
        help="Path to HGNC JSONL file"
    )
    parser.add_argument(
        "enst_ensg",
        type=str,
        nargs=1,
        help="Mapping from ENST to ENSG.",
    )
    parser.add_argument(
        "reference",
        metavar="hg19.fa",
        type=str,
        nargs=1,
        help="Path to FAI-indexed FASTA file for taking contig lengths and names from",
    )
    parser.add_argument(
        "input", metavar="ALIGNMENT.fa.gz", type=str, nargs=1, help="The alignment file to process"
    )
    parser.add_argument(
        "--output",
        metavar="OUT.vcf.gz",
        type=str,
        help="Output file to write to, defaults to stdout",
    )
    parser.add_argument("--verbose", action="store_true", default=False, help="Enable verbose mode")

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
