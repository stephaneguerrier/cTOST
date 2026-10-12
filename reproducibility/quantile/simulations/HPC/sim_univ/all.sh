#!/bin/bash
set -e
source "$(dirname "$0")/setting.sh"

if [ -d "$run_dir" ] && [ -n "$(find "$run_dir" -type f -print -quit)" ]; then
  echo "Previous results or logs in $run_dir; keep them separate before" >&2
  exit 1
fi
mkdir -p "$run_dir/tmp" "$run_dir/report" "$log_dir"
unset SIM_RUN_ID
ID=$(sbatch --parsable --export=ALL --array="1-$n_array" --chdir="$SIM1_DIR" \
  --output="$log_dir/outfile_%a.out" "$SIM1_DIR/launch.sh")
ID=${ID%%;*}
[[ "$ID" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid array job ID: $ID" >&2; exit 1; }
export SIM_RUN_ID="$ID"
sbatch --parsable --export=ALL --dependency="afterok:$ID" --kill-on-invalid-dep=yes \
  --chdir="$SIM1_DIR" --output="$log_dir/outfile_recomb_HPC.out" "$SIM1_DIR/recomb.sh"
