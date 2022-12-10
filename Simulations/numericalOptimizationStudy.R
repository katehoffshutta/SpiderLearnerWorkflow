library(config)
library(ensembleGGM)
library(MASS)

config = config::get(config="pilot_sim_numerical_optimization", file="Simulations/config.yml")
print(config)

load("Results/eightNetworksAC.rda")
eightNetworks = ac
precMat = eightNetworks[[6]] # scale-free high density
print(nrow(precMat))

standardize = function(x){return((x-mean(x))/sd(x))}
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  sample = apply(sample,2,standardize)
  return(sample)
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

candidates = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine)

s = MakeSpiderLearner(candidates)

# for the same input data, are the estimated coefficients the same every time?
# sample the input data
set.seed(0926)

nObs = 1000
thisSample = sampleNetworkData(N=nObs,covMat=solve(precMat))

# make matrix to store results
nIt = config$nSim
ensModelWeights = matrix(rep(NA,nIt*length(candidates)),
                         ncol=length(candidates))

# also store the precision matrices
ensModels = list()

for(i in 1:nIt)
{ 
  print(paste("[numericalOptimizationStudy]: Testing iteration", i))
  ensModel = s$runSpiderLearner(thisSample,config$nFolds,
                              standardize=FALSE,
                              boundedLoss= F, nCores = config$nCores, seedFlag = T, foldseed = 1989)
  ensModelWeights[i,] = ensModel$weights
  write.csv(ensModelWeights, file="Results/Pilot/numericalOptimizationStudyWeights.csv",row.names=F,quote=F)
  ensModels[[i]] = ensModel$optTheta
}

save(ensModels,file="Results/Pilot/numericalOptimizationModels.rda")

