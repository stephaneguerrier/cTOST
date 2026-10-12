#!/bin/bash
#SBATCH --job-name=aqtost_mv
#SBATCH --partition=shared-cpu
#SBATCH --time=0-08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -eC
source "$SIM2_DIR/setting.sh"
module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0
srun Rscript --verbose "$INFILE_sim" "$sim_setting" "$n_simu" "$n_array" > "$OUTFILE_sim" 2>&1
