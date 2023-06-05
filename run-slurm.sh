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

# Default partition.
PART=${PART-critical}

snakemake \
    --jobs 100 \
    --slurm \
    --default-resources \
        slurm_partition=$PART \
        slurm_nodes=1 \
        slurm_ntasks=1 \
        slurm_time=04:00 \
        slurm_mem_mb=2000 \
    -- \
    all
