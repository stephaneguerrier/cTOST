#!/bin/bash
set -e
export SIM1_DIR="${SIM1_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)}"
export sim_setting=${sim_setting:-1}  # 1 = size, 2 = power
export n_simu=${n_simu:-50000}        # number of simulations: 50000
export n_array=${n_array:-5000}       # number of array for HPC: 5000
[[ "$sim_setting" =~ ^[12]$ ]] || { echo 'Use setting 1 (size) or 2 (power).' >&2; exit 1; }

run_dir="$SIM1_DIR/sim_HPC_out/sim_setting_$sim_setting"
log_dir="$run_dir/outfile"
if [[ "$sim_setting" == 1 ]]; then
  result_name="results_sim_setting_size_1_n_simu_${n_simu}_n_array_${n_array}.rdata"
else
  result_name="results_sim_setting_2_n_simu_${n_simu}_n_array_${n_array}.rdata"
fi
export INFILE_sim="$SIM1_DIR/sim_HPC_main.R"
export OUTFILE_sim="$run_dir/report/report_${SLURM_ARRAY_TASK_ID:-console}.Rout"
export INFILE_recomb="$SIM1_DIR/out_HPC_combine.R"
export OUTFILE_recomb="$run_dir/report_recomb_HPC.Rout"
