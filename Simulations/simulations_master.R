#library(devtools)
#install_github("katehoffshutta/ensembleGGM")
library(ensembleGGM)

source("errorMetrics.R")
source("generateSimGraphs.R")
source("simulations_config.R")

sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  return(sample)
}

runSimulation = function(simConfig, whichNetworks = seq(1:8))
{
  M = simConfig$.nModels
  P = simConfig$.nPred
  NS = simConfig$.nSim
  
  set.seed(42) # this seed is ALWAYS 42
  eightNetworks = makeGoldStandardNets(P)
  print(M)
  print(NS)
  relFrobNormsAfter = matrix(rep(NA,(M+2)*NS),nrow=NS)
  matrixRVs = matrix(rep(NA,(M+2)*NS),nrow=NS)
  lls = matrix(rep(NA,(M+2)*NS),nrow=NS)
  llsTest = matrix(rep(NA,(M+2)*NS),nrow=NS)
  ensModels = list()
  
  for(index in whichNetworks)
  {
    set.seed(simConfig$.seed)
    thisNetwork = eightNetworks[[index]]
    
    s = MakeSpiderLearner(simConfig$.candidates)
    
    for(j in 1:NS)
    {
      print(paste("[Running Simulation Number]",j))
      
      thisSample = sampleNetworkData(N=simConfig$.nObs,covMat=solve(thisNetwork)) # sample network data function includes standardizing
      thisTestSample = sampleNetworkData(N=simConfig$.nObs,covMat=solve(thisNetwork))
      
      ensModel = s$runSpiderLearner(thisSample,K=simConfig$.nFolds,standardize=FALSE,nCores=simConfig$.nCores)
      models = ensModel$fullModels
      
      for(i in 1:M)
      {
        relFrobNormsAfter[j,i] = relativeFrobNormAfter(models[[i]], thisNetwork)
        matrixRVs[j,i] = matrixRV(models[[i]], thisNetwork)
        lls[j,i] = loglikLossfunction(models[[i]],thisSample)
        llsTest[j,i] = loglikLossfunction(models[[i]],thisTestSample)
      }
      
      ## Store metrics for ensemble model
      
      relFrobNormsAfter[j,(M+1)] =  relativeFrobNormAfter(ensModel$optTheta,thisNetwork)
      matrixRVs[j,(M+1)] = matrixRV(ensModel$optTheta,thisNetwork)
      lls[j,(M+1)] = loglikLossfunction(ensModel$optTheta,thisSample)
      llsTest[j,(M+1)] = loglikLossfunction(ensModel$optTheta,thisTestSample)
      
      ## Store metrics for simple mean model
      relFrobNormsAfter[j,(M+2)] =  relativeFrobNormAfter(ensModel$simpleMeanNetwork,thisNetwork)
      matrixRVs[j,(M+2)] = matrixRV(ensModel$simpleMeanNetwork,thisNetwork)
      lls[j,(M+2)] = loglikLossfunction(ensModel$simpleMeanNetwork,thisSample)
      llsTest[j,(M+2)] = loglikLossfunction(ensModel$simpleMeanNetwork,thisTestSample)
      
      ## store ensemble model itself
      ensModels[[j]] = ensModel
      
      ## Save at every iteration of the sims.
      ## Not efficient because of the copying,
      ## But helps in case the simulation is interrupted
      theseResults = list("ensModels"=ensModels,
                          "rfnAfter"=relFrobNormsAfter,
                          "mrv"=matrixRVs,
                          "llTrain"=lls,
                          "llTest"=llsTest)
      
      save(theseResults,
           file=paste(c(simConfig$.today,"_", suffixes[index], "_n_", simConfig$.nObs, "_p_",P,"_simStudy.rda"),collapse=""))
    }
  }
}

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new()
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)

candidates_ld = list(apple,
                     banana,
                     clementine,
                     date,
                     elderberry,
                     fraise,
                     grape,
                     honeydew,
                     icewine)

candidates_hd = candidates_ld[-5]

suffixes = list("erLowPrec_RealData",
                "erHighPrec_RealData",
                "wsLowPrec_RealData",
                "wsHighPrec_RealData",
                "sfLowPrec_RealData",
                "sfHighPrec_RealData",
                "hsLowPrec_RealData",
                "hsHighPrec_RealData")

# initialize = function(candidates,nPred,nObs,nFolds,nSim,nCores,seed)

simAConfigPilot = SimConfig$new(candidates = candidates_ld,
                           nPred = 50,
                           nObs = 10000,
                           nFolds = 2,
                           nSim = 2,
                           nCores = 1,
                           seed = 221)

simBConfig = SimConfig$new(candidates = candidates_ld,
                           nPred = 50,
                           nObs = 1600,
                           nFolds = 10,
                           nSim = 1,
                           nCores = 1,
                           seed = 221)

simCConfig = SimConfig$new(candidates = candidates_ld,
                           nPred = 50,
                           nObs = 100,
                           nFolds = 10,
                           nSim = 1,
                           nCores = 1,
                           seed = 221)

simDConfig = SimConfig$new(candidates = candidates_hd,
                           nPred = 100,
                           nObs = 60,
                           nFolds = 10,
                           nSim = 1,
                           nCores = 1,
                           seed = 221)

simAResults = runSimulation(simAConfigPilot, whichNetworks=1)
