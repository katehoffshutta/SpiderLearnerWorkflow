source("simulations_master.R")
library(config)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
config_name = args[1]

# read in config parameters from config file
config = config::get(config=config_name, file="Simulations/config.yml")
print(config)

if(config$candidates == "candidates_ld")
  my_candidates = candidates_ld

if(config$candidates == "candidates_hd")
  my_candidates = candidates_hd

if(config$candidates == "candidates_ld_clime")
  my_candidates = candidates_ld_clime

if(config$candidates == "candidates_hd_clime")
  my_candidates = candidates_hd_clime

simCConfigPilot = SimConfig$new(candidates = my_candidates,
                                nPred = 50,
                                nObs = 100,
                                nFolds = config$nFolds,
                                nSim = config$nSim,
                                nCores = config$nCoresInner,
                                seed = config$seed)


# Parallelize over the networks
# 8 processes which need 5 cores each
# Request total of 40 cores
# Memory per core ~400

doParallel::registerDoParallel(cores = 8)
foreach(i=1:8) %dopar% 
{
    print(paste("foreach",i))
    runSimulation(simCConfigPilot, whichNetworks=i)
}

doParallel::stopImplicitCluster()

