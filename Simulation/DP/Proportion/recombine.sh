#!/bin/bash
#SBATCH --job-name=recombine_dp_prop
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=public-cpu,public-bigmem,shared-cpu,shared-bigmem
#SBATCH --mail-user=stephane.guerrier@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --output ctost/proportion/outfile/outfile_recombine.out
module load GCC/13.2.0 R/4.4.1
INFILE=ctost/proportion/recombine.R
OUTFILE=ctost/proportion/report/recombine.Rout
srun R CMD BATCH $INFILE $OUTFILE
