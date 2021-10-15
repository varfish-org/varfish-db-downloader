#!/usr/bin/bash

export SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

set -euo pipefail

OUT=$1
IN=$2
ECHO=${3-echo}

mkdir -p slurm_log

for table_group in $(tail -n +2 $IN  | sort -u | cut -f 2); do
    if [[ -e $OUT/*.$table_group.*json ]]; then
        >&2 echo "skipping $table_group"
        continue
    fi
    $ECHO sbatch --partition=critical "--job-name=acumenify-$table_group" $SCRIPT_DIR/sbatch_acumenify_job.sh "$OUT" "$IN" "$table_group"
done
