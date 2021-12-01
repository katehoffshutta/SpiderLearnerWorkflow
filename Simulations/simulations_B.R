source("simulations_master.R")
library(doParallel)

simBConfigPilot = SimConfig$new(candidates = candidates_ld,
                                nPred = 50,
                                nObs = 1600,
                                nFolds = 10,
                                nSim = 100,
                                nCores = 5,
                                seed = 221)


# Parallelize over the networks
# 8 processes which need 5 cores each
# Request total of 40 cores
# Memory per core ~400

doParallel::registerDoParallel(cores = 8)
foreach(i=1:8) %dopar% 
{
    print(paste("foreach",i))
    runSimulation(simBConfigPilot, whichNetworks=i)
}

doParallel::stopImplicitCluster()
