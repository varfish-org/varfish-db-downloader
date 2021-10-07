#!/usr/bin/bash

#SBATCH --output=slurm_log/slurm-%j.log
#
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=1G

OUT=$1
IN=$2
ONLY=$3

python $SCRIPT_DIR/acumenify.py extract "$OUT" "$IN" --only "$ONLY"
