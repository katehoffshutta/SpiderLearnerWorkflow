# Test libraries of three different sizes
source("../spiderLearnerOOP.R")
source("../spiderLearnerCandidates.R")
library(MASS)

load("../transfer/eightNetworksAC.rda")
precMat = eightNetworks[[6]] # scale-free high density

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
nSim = 100
for(i in 1:nSim)
{
  print(paste("[Library Simulation]",i))

  data = mvrnorm(1000, mu = rep(0,nrow(precMat)), solve(precMat))
  
  # start with a minimal library with just HGL and MLE
  s = SpiderLearner$new()
  
  s$addCandidate(date)
  s$addCandidate(elderberry)
  
  ensMod1 = s$runSpiderLearner(data,K = 10,FALSE,nCores=10)

  # next, add the EBIC
  
  s$addCandidate(apple)
  s$addCandidate(banana)
  
  ensMod2 = s$runSpiderLearner(data,K = 10,FALSE,nCores=10)
  
  # next, add RIC/StARS/everything else
  
  s$addCandidate(clementine)
  s$addCandidate(fraise)
  s$addCandidate(grape)
  s$addCandidate(honeydew)
  s$addCandidate(icewine)
  
  ensMod3 = s$runSpiderLearner(data,K = 10,FALSE,nCores=10)
  
  ensMods = list(ensMod1, ensMod2, ensMod3)
  
  librarySimRes[[i]] = ensMods
}

save(librarySimRes,file="librarySimRes.rda")
