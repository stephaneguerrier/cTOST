#!/bin/bash
#SBATCH --job-name=dp_prop_tost
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=public-cpu,public-bigmem,shared-cpu,shared-bigmem
#SBATCH --mail-user=stephane.guerrier@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --output ctost/proportion/outfile/outfile_%a.out
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1
INFILE=ctost/proportion/simu.R
OUTFILE=ctost/proportion/report/report_${SLURM_ARRAY_TASK_ID}.Rout
srun R CMD BATCH $INFILE $OUTFILE
