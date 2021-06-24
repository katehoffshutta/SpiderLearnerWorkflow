library(hglasso)
library(qgraph)

source("errorMetrics.R")
source("spiderLearnerOOP.R")
source("spiderLearnerCandidates.R")

sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
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

# run Simulation B

nModels = length(candidates)
nPred = 50
nObs = 1600
nFolds = 10
nSim = 100
today = format(Sys.Date(),format="%Y%m%d")

set.seed(42)
source("generateSimGraphs.R")
eightNetworks = makeGoldStandardNets(nPred)

suffixes = list("erLowPrec_RealData",
                "erHighPrec_RealData",
                "wsLowPrec_RealData",
                "wsHighPrec_RealData",
                "sfLowPrec_RealData",
                "sfHighPrec_RealData",
                "hsLowPrec_RealData",
                "hsHighPrec_RealData")

args = commandArgs(trailingOnly=TRUE)

relFrobNormsAfter = matrix(rep(NA,(nModels+2)*nSim),nrow=nSim)
matrixRVs = matrix(rep(NA,(nModels+2)*nSim),nrow=nSim)
lls = matrix(rep(NA,(nModels+2)*nSim),nrow=nSim)
llsTest = matrix(rep(NA,(nModels+2)*nSim),nrow=nSim)
ensModels = list()

for(index in as.numeric(args))
{	
  mySeed = as.numeric(today)
  set.seed(mySeed)
  print(paste("Seed",mySeed))

  thisNetwork = eightNetworks[[index]]

  for(j in 1:nSim)
  {
    print(paste("[Running Simulation Number]",j))
    
    thisSample = sampleNetworkData(N=nObs,covMat=solve(thisNetwork)) # sample network data function includes standardizing
    thisTestSample = sampleNetworkData(N=nObs,covMat=solve(thisNetwork))
  
    # run one fold per core
    # be careful with this!
    ensModel = s$runSpiderLearner(thisSample,K=nFolds,standardize=FALSE,nCores=nFolds)
    models = ensModel$fullModels
    
    for(i in 1:nModels)
    {
      relFrobNormsAfter[j,i] = relativeFrobNormAfter(models[[i]], thisNetwork)
      matrixRVs[j,i] = matrixRV(models[[i]], thisNetwork)
      lls[j,i] = loglikLossfunction(models[[i]],thisSample)
      llsTest[j,i] = loglikLossfunction(models[[i]],thisTestSample)
    }
    
    ## Store metrics for ensemble model

    relFrobNormsAfter[j,(nModels+1)] =  relativeFrobNormAfter(ensModel$optTheta,thisNetwork)
    matrixRVs[j,(nModels+1)] = matrixRV(ensModel$optTheta,thisNetwork)
    lls[j,(nModels+1)] = loglikLossfunction(ensModel$optTheta,thisSample)
    llsTest[j,(nModels+1)] = loglikLossfunction(ensModel$optTheta,thisTestSample)

    ## Store metrics for simple mean model
    relFrobNormsAfter[j,(nModels+2)] =  relativeFrobNormAfter(ensModel$simpleMeanNetwork,thisNetwork)
    matrixRVs[j,(nModels+2)] = matrixRV(ensModel$simpleMeanNetwork,thisNetwork)
    lls[j,(nModels+2)] = loglikLossfunction(ensModel$simpleMeanNetwork,thisSample)
    llsTest[j,(nModels+2)] = loglikLossfunction(ensModel$simpleMeanNetwork,thisTestSample)

    ## store ensemble model itself
    ensModels[[j]] = ensModel
  }
  
  theseResults = list("ensModels"=ensModels,
                           "rfnAfter"=relFrobNormsAfter,
                           "mrv"=matrixRVs,
                           "llTrain"=lls,
                           "llTest"=llsTest)
  
  save(theseResults,
       file=paste(c(today,"_", suffixes[index], "_n_", nObs, "_p_",nPred,"_simStudy.rda"),collapse=""))

}
