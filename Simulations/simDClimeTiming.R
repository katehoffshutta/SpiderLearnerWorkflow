# Show adding CLIME to the ensemble is unreasonably long for
# simulation d

library(clime)
library(ensembleGGM)
library(MASS)

load("../Results/eightNetworksD.rda")
eightNetworks = d
precMat = eightNetworks[[2]] # random graph high density

print(nrow(precMat))
# plot(graph_from_adjacency_matrix(-cov2cor(as.matrix(precMat)),weighted=T,diag=F,mode="undirected"),layout=layout_with_lgl)

standardize = function(x){return((x-mean(x))/sd(x))}
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  sample = apply(sample,2,standardize)
  return(sample)
}

set.seed(42)
mysamp = sampleNetworkData(100,solve(precMat))

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new()
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)
jicama = CLIMECandidate$new()

candidates_clime = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  # elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine,
                  jicama)

candidates_no_clime = list(apple, 
                             banana, 
                             clementine, 
                             date, 
                             # elderberry,
                             fraise,
                             grape,
                             honeydew,
                             icewine)

s_no_clime = MakeSpiderLearner(candidates_no_clime)
s_clime = MakeSpiderLearner(candidates_clime)


start_no_clime = proc.time()
myres_no_clime = s_no_clime$runSpiderLearner(mysamp,10,
                           standardize = F, 
                           boundedLoss = F,
                           nCores = 5)
stop_no_clime = proc.time()

start_clime = proc.time()
myres_clime = s_clime$runSpiderLearner(mysamp,10,
                                       standardize = F, 
                                      boundedLoss = F,
                                      nCores = 5)
stop_clime = proc.time()

stop_clime - start_clime
#user    system   elapsed 
#20909.751    39.081 12128.336 

stop_no_clime - start_no_clime
# user   system  elapsed 
# 1153.180    5.843  313.125 

round(myres_clime$weights,3)
# [1] 0.000 0.000 0.000 0.169 0.000 0.355 0.477 0.000 0.000
round(myres_no_clime$weights,3)
# [1] 0.000 0.000 0.000 0.158 0.000 0.347 0.495 0.000


