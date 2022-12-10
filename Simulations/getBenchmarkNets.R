# This has to be done on the cluster because the random seed does not reproduce the same results
# with R on different versions of R.

set.seed(42)
source("Simulations/generateSimGraphs.R")
ac = makeGoldStandardNets(50)
save(ac, file="Results/eightNetworksAC.rda")

set.seed(42)
source("Simulations/generateSimGraphs_simD.R")
d = makeGoldStandardNets_simD(100)
save(d, file="Results/eightNetworksD.rda")