# Show adding CLIME to the ensemble

library(clime)
library(ensembleGGM)
library(MASS)

load("../Results/eightNetworksAC.rda")
eightNetworks = ac
precMat = eightNetworks[[2]] # random graph high density
precMat = eightNetworks[[6]] # scale free graph high density

print(nrow(precMat))
# plot(graph_from_adjacency_matrix(-cov2cor(as.matrix(precMat)),weighted=T,diag=F,mode="undirected"),layout=layout_with_lgl)

standardize = function(x){return((x-mean(x))/sd(x))}
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  sample = apply(sample,2,standardize)
  return(sample)
}

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

candidates = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  # elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine,
                  jicama)

candidates = list(icewine, jicama)

s = MakeSpiderLearner(candidates)

myres = s$runSpiderLearner(mysamp,2,
                           standardize = F, 
                           boundedLoss = F)
# climemat = clime(mysamp)
# re.cv= cv.clime(climemat)
# re.clime.opt = clime(mysamp, lambda = re.cv$lambdaopt)

