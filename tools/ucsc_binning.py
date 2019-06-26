import sys
import binning
import argparse
from collections import OrderedDict


def tsv_reader(fh, header):
    keys = header.rstrip("\n").split("\t")
    for line in fh:
        if not line.startswith("#"):
            yield OrderedDict(zip(keys, line.rstrip("\n").split("\t")))


def run(args):
    header = next(args.input)
    args.output.write(header)
    for record in tsv_reader(args.input, header):
        if record["end"] == "":
            record["end"] = str(int(record["start"]) + len(record["reference"]) - 1)
        record["bin"] = str(binning.assign_bin(int(record["start"]) - 1, int(record["end"])))
        args.output.write("%s\n" % "\t".join(v for v in record.values()))


def main(argv=None):
    parser = argparse.ArgumentParser(
        "Fill empty ``bin`` column in TSV file with values. "
        "``start`` and ``end`` column must be present in TSV file."
    )
    parser.add_argument(
        '--input',
        type=argparse.FileType('r'),
        nargs='?',
        default=sys.stdin,
        help="Input tsv file to fill the binning column."
    )
    parser.add_argument(
        '--output',
        type=argparse.FileType('w'),
        nargs='?',
        default=sys.stdout,
        help="Output tsv file with filled binning column."
    )
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())

