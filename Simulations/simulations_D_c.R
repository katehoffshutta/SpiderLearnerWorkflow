source("simulations_master.R")
library(doParallel)
source("generateSimGraphs_simD.R")

simDConfigPilot = SimConfig$new(candidates = candidates_hd,
                                nPred = 100,
                                nObs = 60,
                                nFolds = 10,
                                nSim = 100,
                                nCores = 10,
                                seed = 221)


# Parallelize over the networks
# 2 processes which need 10 cores each
# Request total of 20 cores
# Memory per core ~400

doParallel::registerDoParallel(cores = 2)
foreach(i=c(4,8)) %dopar% 
{
    print(paste("foreach",i))
    runSimulation(simDConfigPilot, whichNetworks=i)
}

doParallel::stopImplicitCluster()

