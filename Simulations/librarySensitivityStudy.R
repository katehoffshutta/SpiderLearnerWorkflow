# Test libraries of three different sizes
library(config)
library(ensembleGGM)
library(MASS)

config = config::get(config="pilot_sim_library", file="Simulations/config.yml")
print(config)

load("Results/eightNetworksAC.rda")
eightNetworks = ac
precMat = eightNetworks[[6]] # scale-free high density
print(nrow(precMat))

# define all the possible candidates
apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new()
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)

librarySimRes = list()
set.seed(1979)
nSim = config$nSim
for(i in 1:nSim)
{
  print(paste("[Library Simulation]",i))

  data = mvrnorm(1000, mu = rep(0,nrow(precMat)), solve(precMat))
  
  # start with a minimal library with just HGL and MLE
  s = SpiderLearner$new()
  
  s$addCandidate(date)
  s$addCandidate(elderberry)
  
  print("Minimal library")  
  ensMod1 = s$runSpiderLearner(data,K = config$nFolds,FALSE,nCores=config$nCores)

  # next, add the EBIC
  
  s$addCandidate(apple)
  s$addCandidate(banana)
  
  print("Medium library")
  ensMod2 = s$runSpiderLearner(data,K = config$nFolds,FALSE,nCores=config$nCores)
  
  # next, add RIC/StARS/everything else
  
  s$addCandidate(clementine)
  s$addCandidate(fraise)
  s$addCandidate(grape)
  s$addCandidate(honeydew)
  s$addCandidate(icewine)
  
  print("Maximal library")
  ensMod3 = s$runSpiderLearner(data,K = config$nFolds,FALSE,nCores=config$nCores)
  
  ensMods = list(ensMod1, ensMod2, ensMod3)
  
  librarySimRes[[i]] = ensMods
  save(librarySimRes,file="Results/Pilot/librarySimRes.rda")
}

