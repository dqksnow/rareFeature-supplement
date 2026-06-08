#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --array=0-79
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --job-name=rareF_CovType
#SBATCH --export=ALL

module purge
module load r/4.4.2

SETTING_ID=${SLURM_ARRAY_TASK_ID}
OUTPUT_DIR="../raw_results/$(date +%m%d%Y)/CoverType"
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

NCORES=$(nproc)
export OMP_NUM_THREADS=$NCORES
export SLURM_CPUS_PER_TASK=$NCORES

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export R_LIBS_USER=~/rlibs

echo "========================================="
echo "Starting CoverType setting $SETTING_ID at $(date)"
echo "Node: $(hostname)"
echo "Detected $NCORES cores available."
echo "Output directory: $OUTPUT_DIR"
echo "========================================="

Rscript rareF_Covertype_HPC_supplement.R \
  --setting-id "$SETTING_ID" \
  --output-dir "$OUTPUT_DIR"

echo "========================================="
echo "Finished CoverType setting $SETTING_ID at $(date)"
echo "========================================="
