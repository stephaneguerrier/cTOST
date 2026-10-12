#!/bin/bash
set -e
export SIM2_DIR="${SIM2_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)}"
export sim_setting=${sim_setting:-2}  # setting number: 1 (n=30), 2 (n=10), 3 (n=50)
export n_simu=${n_simu:-50000}        # number of simulations: 50000
export n_array=${n_array:-5000}       # number of arrays for HPC: 5000
[[ "$sim_setting" =~ ^[123]$ ]] || { echo 'Use setting 1, 2 or 3.' >&2; exit 1; }

run_dir="$SIM2_DIR/sim_HPC_out/sim_setting_$sim_setting"
log_dir="$run_dir/outfile"
result_name="results_mvt_sim_setting_${sim_setting}_n_simu_${n_simu}_n_array_${n_array}.rdata"
export INFILE_sim="$SIM2_DIR/sim_HPC_main.R"
export OUTFILE_sim="$run_dir/report/report_${SLURM_ARRAY_TASK_ID:-console}.Rout"
export INFILE_recomb="$SIM2_DIR/out_HPC_combine.R"
export OUTFILE_recomb="$run_dir/report_recomb_HPC.Rout"
