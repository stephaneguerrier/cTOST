#!/bin/bash
set -e
source "${SIM1_DIR:-$(dirname "$0")}/setting.sh"

{
  read -r completed_run
  read -r completed_file
} < "$run_dir/.cleanup_ok"
if [ "$completed_run" != "${SIM_RUN_ID:?Missing current array ID}" ] || [ "$completed_file" != "$result_name" ]; then
  echo 'Completion marker does not match' >&2
  exit 1
fi
if [ ! -f "$run_dir/$result_name" ]; then
  echo 'Combined result missing' >&2
  exit 1
fi

find "$run_dir/tmp" "$run_dir/report" "$log_dir" -maxdepth 1 -type f -delete
