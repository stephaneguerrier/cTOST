#!/bin/sh
for epsilon_val in 0.1 0.25 0.5 1
  do
  for n_val in 200 300 400 500 600 700 800 900 1000
    do
    for p1_val in 0.5 0.65 0.8
      do
      for scenario_val in 1 2
        do
        eval "export epsilon_val=$epsilon_val"
        eval "export n_val=$n_val"
        eval "export p1_val=$p1_val"
        eval "export scenario_val=$scenario_val"
        sbatch --array=1-3 ctost/proportion/simu.sh
      done
    done
  done
done
