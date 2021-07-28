import glob
import os
import re
import sys
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
print("\n---\n%s\n---\n" % json.dumps(config, indent="  "))

#: Use strict mode and also print each command.
shell.prefix("set -x; set -euo pipefail; ")

#: The canonical chromosome names.
CHROMS = list(range(1, 23)) + ["X", "Y"]

#: List for collecting all result files below.
ALL_RESULT = []


def input_all(wildcards):
    return ALL_RESULT


rule all:
    input:
        input_all,


# Load all snakemake files `snakefiles/*/*.smk`.
snakefiles = list(sorted(glob.glob("snakefiles/*/*.smk")))
print("Loading Snakefiles...\n  - %s" % "\n  - ".join(snakefiles), file=sys.stderr)

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
                for chrom in chroms:
                    print(rule.name, chrom, vals, path)
                    path = re.sub(r"{([^,]+)(,.*)?}", r"{\1}", path)
                    ALL_RESULT.append(path.format(chrom=chrom, **vals))

print("  - " + "\n  - ".join(sorted(ALL_RESULT)), file=sys.stderr)
