library(affy)
library(ensembleGGM)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Uncomment this line to install the curatedOvarianData package
#BiocManager::install("curatedOvarianData")

library(curatedOvarianData)
standardize = function(x){return((x-mean(x))/sd(x))}

data(GSE32062.GPL6480_eset)
lateStage = exprs(GSE32062.GPL6480_eset)
yoshi = read.table("YoshiharaGeneSet.tsv",sep="&")
lateStageSmall = lateStage[which(rownames(lateStage)%in%yoshi[,1]),]
lateStageSmall = t(lateStageSmall)
names(lateStageSmall) = colnames(lateStageSmall)
lateStageSmall = apply(lateStageSmall,2,standardize)

s = SpiderLearner$new()

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new() # don't add MLE forHD case
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

set.seed(1206)
trainIndices = sample(rep(1:10,26))
loglikLossfunction = function(thetaEst,dataTest)
{
  n = nrow(dataTest)
  p = ncol(dataTest)
  firstTerm = -n/2*log(det(solve(thetaEst)))
  secondTerm= 0
  for(i in 1:n)
    secondTerm = secondTerm - 1/2*t(dataTest[i,]) %*% thetaEst %*% dataTest[i,]
  return(firstTerm + secondTerm)
}

ensLoss = rep(NA,10)
candidateLoss = matrix(rep(NA,90),nrow=10,ncol=9)
for(k in 1:10)
{
  print(paste("Working in fold:",k))
  lateStageTrain = lateStageSmall[trainIndices != k,]
  lateStageTest = lateStageSmall[trainIndices == k,]
  slResultTrain = s$runSpiderLearner(lateStageTrain, K = 10, nCores = 10)
  ensLoss[k] = loglikLossfunction(slResultTrain$optTheta,lateStageTest)
  candidateLoss[k,] = sapply(slResultTrain$fullModels,loglikLossfunction,lateStageTest)
  print("ensLoss")
  print(ensLoss)
  print("candidateLoss")
  print(candidateLoss)
}

loglikLossfunction = function(thetaEst,dataTest)
{
  n = nrow(dataTest)
  p = ncol(dataTest)
  firstTerm = -n/2*log(det(solve(thetaEst)))
  secondTerm= 0
  for(i in 1:n)
    secondTerm = secondTerm - 1/2*t(dataTest[i,]) %*% thetaEst %*% dataTest[i,]
  return(firstTerm + secondTerm)
}

pdf("ovarianSmallOOSL.pdf",width=10,height=6)
boxplot(cbind(candidateLoss,ensLoss),
        names=c(paste("Method",1:9),"Ensemble"),
        col = rainbow(10),
        las = 2,
        main = "Out-of-sample Log Likelihood: 10-fold CV")
dev.off()
