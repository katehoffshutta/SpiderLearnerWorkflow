library(ggplot2)
library(huge)
library(igraph)
library(MASS)
library(moments)
library(pracma)
library(ensembleGGM)


standardize = function(x){return((x-mean(x))/sd(x))}
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  sample = apply(sample,2,standardize)
  return(sample)
}

s = SpiderLearner$new()

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new()
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)

candidates = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine)

for(candidate in candidates)
{
  s$addCandidate(candidate)
}


testRuntime = function(nPred, nObs, nFolds, eightNetworks,sl,index=1)
{  
  X = sampleNetworkData(N = nObs,covMat = solve(eightNetworks[[index]]))
  start_time = Sys.time()
  ensModel = s$runSpiderLearner(X,K=nFolds,standardize=FALSE)
  end_time = Sys.time()
  print(end_time - start_time)
  return(end_time - start_time)
}


source("generateSimGraphs.R")

set.seed(42)
eightNetworks = makeGoldStandardNets(25)

runtimes = list()
print("Runtime 1")
runtimes[[1]] = testRuntime(25,10000,10,eightNetworks,s,1)
print("Runtime 2")
runtimes[[2]] = testRuntime(25,1000,10,eightNetworks,s,1)
print("Runtime 3")
runtimes[[3]] = testRuntime(25,100,10,eightNetworks,s,1)

set.seed(42)
eightNetworks = makeGoldStandardNets(75)

print("Runtime 4")
runtimes[[4]] = testRuntime(75,10000,10,eightNetworks,s,1)
print("Runtime 5")
runtimes[[5]] = testRuntime(75,1000,10,eightNetworks,s,1)
print("Runtime 6")
runtimes[[6]] = testRuntime(75,100,10,eightNetworks,s,1)

set.seed(42)
eightNetworks = makeGoldStandardNets(100)

runtimes[[7]] = testRuntime(100,10000,10,eightNetworks,s,1)
print("Runtime 8")
runtimes[[8]] = testRuntime(100,1000,10,eightNetworks,s,1)

s$removeCandidate("mle")

print("Runtime 9")
runtimes[[9]] = testRuntime(100,100,10,eightNetworks,s,1)

set.seed(42)
print("Runtime 10")
runtimes[[10]] = testRuntime(150,100,10,eightNetworks,s,1)
