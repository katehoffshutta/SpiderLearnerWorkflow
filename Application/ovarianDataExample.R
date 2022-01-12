# This R script runs the ovarian cancer example
# in the application section of the paper

# You will need to install ensembleGGM using install_github:
# library(devtools)
# install_github("katehoffshutta/ensembleGGM")

library(affy)
library(ensembleGGM)
library(tidyverse)

candidateNames = c("glasso-ebic-0", 
                   "glasso-ebic-0.5", 
                   "glasso-ric", 
                   "hglasso",
                   "mle", 
                   "glasso-stars-0.05",
                   "glasso-stars-0.1",
                   "qgraph-ebic-0", 
                   "qgraph-ebic-0.5")


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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Uncomment this line to install the curatedOvarianData package
# BiocManager::install("curatedOvarianData")

library(curatedOvarianData)
standardize = function(x){return((x-mean(x))/sd(x))}

data(GSE32062.GPL6480_eset)
lateStage = exprs(GSE32062.GPL6480_eset)
yoshi = read.table("yoshiharaLimitedGeneSet.tsv",sep="&")
lateStageSubset = lateStage[which(rownames(lateStage)%in%yoshi[,1]),]
lateStageSmall = apply(data.frame(t(lateStageSubset)),2,standardize)
#colnames(lateStageSmall)[which(colnames(lateStageSmall)=="HLA-DPB1")] = "HLA.DPB1"
names(lateStageSmall)=colnames(lateStageSmall)

# Run SpiderLearner

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

#set.seed(1210)
#slResults = s$runSpiderLearner(lateStageSmall, 
#                               K = 10, 
#                               standardize=T, 
#                               nCores = 10)

#save(slResults,file="../Figures/ocSpiderLearnerResults.rda")

load("../Figures/ocSpiderLearnerResults.rda")

# Load validation datasets

data("E.MTAB.386_eset")
data("GSE13876_eset")
data("GSE14764_eset")
data("GSE17260_eset")
data("GSE18520_eset")
data("GSE19829.GPL570_eset")
data("GSE19829.GPL8300_eset")
data("GSE26712_eset")
data("GSE30009_eset")
data("GSE30161_eset")
data("GSE32063_eset")
data("GSE9891_eset")
data("PMID17290060_eset")
data("PMID19318476_eset")
data("TCGA_eset")

identifiers = c("E.MTAB.386",
                "GSE13876",
                "GSE14764",
                "GSE17260",
                "GSE18520",
                "GSE19829.GPL570",
                "GSE19829.GPL8300",
                "GSE26712",
                "GSE30009",
                "GSE30161",
                "GSE32063",
                "GSE9891",
                "PMID17290060",
                "PMID19318476",
                "TCGA")

validData = list(exprs(E.MTAB.386_eset),
                 exprs(GSE13876_eset),
                 exprs(GSE14764_eset),
                 exprs(GSE17260_eset),
                 exprs(GSE18520_eset),
                 exprs(GSE19829.GPL570_eset),
                 exprs(GSE19829.GPL8300_eset),
                 exprs(GSE26712_eset),
                 exprs(GSE30009_eset),
                 exprs(GSE30161_eset),
                 exprs(GSE32063_eset),
                 exprs(GSE9891_eset),
                 exprs(PMID17290060_eset),
                 exprs(PMID19318476_eset),
                 exprs(TCGA_eset))

# Standardize validation datasets

validDataStd = lapply(validData,function(x){t(apply(x,1,standardize))})
validationLikelihoods = matrix(rep(NA,15*10),ncol=10)

genesToGet = colnames(slResults$optTheta)
genesToGet[genesToGet == "HLA.DPB1"] = "HLA-DPB1"
for(d in c(2:6,8,10:15)) # skip datasets 1,7, and 9 as described in selectYoshiharaGeneSubset.R
{
  data = data.frame(validDataStd[[d]])
  data$gene = row.names(data)
  #data$gene[data$gene == "HLA.DPB1"] = "HLA-DPB1"
  matchedData = merge(data.frame("gene"=genesToGet),data,by="gene")
  mismatch = union(setdiff(genesToGet,matchedData$gene),setdiff(matchedData$gene,genesToGet))
  
  if(length(mismatch) == 0)
  {
    print(paste("Dataset",d,": ",loglikLossfunction(slResults$optTheta,t(matchedData[,-1]))),collapse="")
    validationLikelihoods[d,10] = loglikLossfunction(slResults$optTheta,t(matchedData[,-1]))
    for(j in 1:9)
    {
      validationLikelihoods[d,j] = loglikLossfunction(slResults$fullModels[[j]],t(matchedData[,-1]))
    }
  }
  
  else
 {
   print(mismatch)
 }
}

scaledVal = data.frame(t(apply(validationLikelihoods,1,function(x){x<- (x-max(x))/max(abs(x))})))
names(scaledVal) = c(candidateNames,"SpiderLearner")
scaledVal$dataset = 1:15
write.csv(scaledVal,"../Figures/scaledVal.csv")

# We (theoretically) get a better out-of-sample likelihood 
# with the ensemble
# Sample splitting is another way to demo
# the utility of this method

### set.seed(1202)
### trainIndices = sample(rep(1:10,26))

### ensLoss = rep(NA,10)
### candidateLoss = matrix(rep(NA,90),nrow=10,ncol=9)
### for(k in 1:10)
### {
###   print(paste("Working in fold:",k))
###   lateStageTrain = lateStageSmall[trainIndices != k,]
###   lateStageTest = lateStageSmall[trainIndices == k,]
###   slResultTrain = s$runSpiderLearner(lateStageTrain, K = 10, nCores = 10)
###   ensLoss[k] = loglikLossfunction(slResultTrain$optTheta,lateStageTest)
###   candidateLoss[k,] = sapply(slResultTrain$fullModels,loglikLossfunction,lateStageTest)
### }

### df = as.data.frame(cbind(candidateLoss,ensLoss))
### colnames(df)[1:9] = candidateNames
### colnames(df)[10] = "SpiderLearner"
### df$fold = 1:10

### write.csv(df,file="../Figures/internalCVResults.csv")

