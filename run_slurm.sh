#!/bin/bash

# Maximal number of jobs to execute at the same time
MAX_JOBS=1000
# Maximal number of jobs per second
MAX_JOBS_PER_SECOND=10
# Number of times to restart jobs
RESTART_TIMES=5

# Check preconditions -------------------------------------------------------

# Ensure slurm_log is a directory
test -d slurm_log || { >&2 echo "${PWD}/slurm_log does not exist"; exit 1; }

# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Create one log directory per Snakemake run --------------------------------

test -z "${SLURM_JOB_ID-}" && SLURM_JOB_ID=$(date +%Y-%m-%d_%H-%M)
LOGDIR=slurm_log/${SLURM_JOB_ID}
mkdir -p ${LOGDIR}

# Kick off Snakemake --------------------------------------------------------

# Using the medium project/queue is a sensible default.
snakemake \
    --printshellcmds \
    --jobs ${MAX_JOBS} \
    --drmaa " --partition=critical --mem=20000 --time=48:00:00 --cpus-per-task=8 --output=${PWD}/${LOGDIR}/slurm-%x-%J.log" \
    --restart-times ${RESTART_TIMES} \
    --rerun-incomplete \
    -- \
    $*

# Print date after finishing, for good measure ------------------------------

>&2 date
>&2 echo "All done. Have a nice day."
