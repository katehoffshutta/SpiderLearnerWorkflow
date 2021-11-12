# code for simulation config object
library(R6)

SimConfig = R6Class("SimConfig", 
                    public = list(
                      .candidates = list(), # list of Candidate objects
                      .nModels = 0,
                      .nPred = 0,
                      .nObs = 0,
                      .nFolds = 1,
                      .nSim = 0,
                      .nCores = 1,
                      .seed = 42,
                      .today = NULL,
                      initialize = function(candidates,nPred,nObs,nFolds,nSim,nCores,seed)
                      {
                        self$.candidates = candidates
                        self$.nModels = length(self$.candidates)
                        self$.nPred = nPred
                        self$.nObs = nObs
                        self$.nFolds= nFolds
                        self$.nSim = nSim
                        self$.nCores = nCores
                        self$.seed = seed
                        self$.today = format(Sys.Date(),format="%Y%m%d")
                      }
                      
                    )
)