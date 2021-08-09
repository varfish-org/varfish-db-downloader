import glob
from itertools import product
import os
import re
import sys
import textwrap
import json

from snakemake import shell
from tools.sv_db_to_tsv import to_tsv


# Ensure that the configuration file exists and then load it.
if not os.path.exists("config.yaml"):
    print("No config.yaml exists yet. Try `cp config.yaml.example config.yaml`.", file=sys.stderr)
    sys.exit(1)


configfile: "config.yaml"


# Print configuration.
print("Configuration:", file=sys.stderr)
print("\n---\n%s\n---\n" % json.dumps(config, indent="  "), file=sys.stderr)

#: Use strict mode and also print each command.
shell.prefix("set -x; set -euo pipefail; ")

#: The canonical chromosome names (without "Y").
CHROMS_NO_Y = list(map(str, range(1, 23))) + ["X"]

#: The canonical chromosome names.
CHROMS = CHROMS_NO_Y + ["Y"]

#: List for collecting all result files below.
ALL_RESULT = []


def input_all(wildcards):
    return ALL_RESULT


rule all:
    input:
        input_all,


# Load all snakemake files `snakefiles/*/*.smk`.
snakefiles = list(sorted(glob.glob("snakefiles/*/*.smk")))
print("Loading Snakefiles...", file=sys.stderr)

for path in snakefiles:

    include: path


# Derive output overall output files from rules starting with prefix "output_".
print("Constructing list of output files from all `result_*` rules...\n", file=sys.stderr)

ALL_RESULT = []
for rule in workflow.rules:
    if rule.name.startswith("result_"):
        for genome_build in ("GRCh37", "GRCh38"):
            vals = {**config, "genome_build": genome_build}
            for path in rule.output:
                if "{chrom}" in path:
                    chroms = CHROMS
                else:
                    chroms = ["-"]
                if "{chrom_no_y}" in path:
                    chroms_no_y = CHROMS_NO_Y
                else:
                    chroms_no_y = ["-"]
                for chrom, chrom_no_y in product(chroms, chroms_no_y):
                    path = re.sub(r"{([^,]+)(,.*)?}", r"{\1}", path)
                    ALL_RESULT.append(path.format(chrom=chrom, chrom_no_y=chrom_no_y, **vals))
