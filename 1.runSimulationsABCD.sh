#!/bin/bash

echo "Running Simulation A"
Rscript --vanilla Simulations/simulations_A.R pilot_sim_abc > Logs/sim_a_log.txt

echo "Running Simulation B"
Rscript --vanilla Simulations/simulations_B.R pilot_sim_abc > Logs/sim_b_log.txt

echo "Running Simulation C"
Rscript --vanilla Simulations/simulations_C.R pilot_sim_abc > Logs/sim_c_log.txt

echo "Running Simulation D"
Rscript --vanilla Simulations/simulations_D_a.R pilot_sim_d > Logs/sim_d_a_log.txt
Rscript --vanilla Simulations/simulations_D_b.R pilot_sim_d > Logs/sim_d_b_log.txt
