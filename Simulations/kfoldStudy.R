library(config)
library(ensembleGGM)
library(foreach)
library(doParallel)
library(MASS)

# read in config parameters from config file
config = config::get(config="pilot_sim_kfold", file="Simulations/config.yml")

load("Results/eightNetworksAC.rda")
eightNetworks = ac
precMat = eightNetworks[[2]] # Random graph, high density
print(dim(precMat))
data = mvrnorm(150, mu = rep(0,nrow(precMat)), solve(precMat))

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

s$addCandidate(apple)
s$addCandidate(banana)
s$addCandidate(clementine)
s$addCandidate(date)
s$addCandidate(elderberry)
s$addCandidate(fraise)
s$addCandidate(grape)
s$addCandidate(honeydew)
s$addCandidate(icewine)

s$printLibrary()

R=config$nSim
nSimCores = config$nSimCores

registerDoParallel(nSimCores)
simRes = foreach(r=1:R) %dopar%
{ 
    print(paste("[simulation]",r))

    data = mvrnorm(150, mu = rep(0,nrow(precMat)), solve(precMat))
    
    print(paste(c("[simulation]",r,"k=2"),collapse=" "))
    a = s$runSpiderLearner(data,K=2,standardize = FALSE)
    print(paste(c("[simulation]",r,"k=5"),collapse=" "))
    b = s$runSpiderLearner(data,K=5,standardize = FALSE)
    print(paste(c("[simulation]",r,"k=10"),collapse=" "))
    c = s$runSpiderLearner(data,K=10,standardize = FALSE)
    print(paste(c("[simulation]",r,"k=15"),collapse=" "))
    d = s$runSpiderLearner(data,K=15,standardize = FALSE)
    print(paste(c("[simulation]",r,"k=20"),collapse=" "))
    e = s$runSpiderLearner(data,K=20,standardize = FALSE)
    print(paste(c("[simulation]",r,"k=30"),collapse=" "))
    f = s$runSpiderLearner(data,K=30,standardize = FALSE)
    
    list(a,b,c,d,e,f)
}

stopImplicitCluster()

save(simRes,file="Results/Pilot/k_fold_sim.rda")