#!/bin/bash

echo "Running Simulation A"
Rscript --vanilla Simulations/simulations_A.R pilot_sim_abc_clime > Logs/sim_a_clime_log.txt

echo "Running Simulation B"
Rscript --vanilla Simulations/simulations_B.R pilot_sim_abc_clime > Logs/sim_b_clime_log.txt

echo "Running Simulation C"
Rscript --vanilla Simulations/simulations_C.R pilot_sim_abc_clime > Logs/sim_c_clime_log.txt

echo "Running Simulation D"
Rscript --vanilla Simulations/simulations_D_a.R pilot_sim_d_clime > Logs/sim_d_a_clime_log.txt
Rscript --vanilla Simulations/simulations_D_b.R pilot_sim_d_clime > Logs/sim_d_b_clime_log.txt
