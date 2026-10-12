#!/bin/bash
#SBATCH --job-name=recomb
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:00:00
#SBATCH --mem=16G
#SBATCH --partition=shared-cpu

set -eC
source "$SIM2_DIR/setting.sh"
module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0
srun Rscript --verbose "$INFILE_recomb" "$sim_setting" "$n_simu" "$n_array" > "$OUTFILE_recomb" 2>&1
bash "$SIM2_DIR/clean.sh"
