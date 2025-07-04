#!/usr/bin/bash

# Example call to run the workflow on a Slurm cluster.

#SBATCH --job-name=varfish-db-downloader
#SBATCH --output=slurm-varfish-db-downloader-%j.out
#
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2-00:00:00
#SBATCH --memory=2G

set -x
set -euo pipefail

# Number of jobs to run at the same time.
JOBS=${JOBS-500}
# Be relaxed with reruns.
RELAXED_RERUNS=${RELAXED_RERUNS-true}
# Whether to add --keep-going
KEEP_GOING=${KEEP_GOING-false}
# Wait time for files to appear in seconds
LATENCY_WAIT=${LATENCY_WAIT-60}

# Other --default-resource parameters:
#   slurm_partition="$PART"

snakemake \
    --rerun-incomplete \
    $(if [[ "$RELAXED_RERUNS" == true ]]; then \
        echo --rerun-triggers mtime; \
        echo --rerun-triggers params; \
        echo --rerun-triggers input; \
    fi) \
    $(if [[ "$KEEP_GOING" == true ]]; then \
        echo --keep-going; \
    fi) \
    --jobs $JOBS \
    --slurm \
    --latency-wait $LATENCY_WAIT \
    --default-resources \
        'runtime="4h"' \
        mem_mb=4000 \
    -- \
    "${@-all}"
