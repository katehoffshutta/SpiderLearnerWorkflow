library(devtools)
install_github("katehoffshutta/ensembleGGM")

library(affy)
library(ensembleGGM)
library(tidyverse)

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
# omit COX17 from score for purposes of validation
# most validation datasets don't distinguish between the pseudogene and gene
lateStageSmall = lateStageSmall[,-which(colnames(lateStageSmall) == "COX17")] 

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
#slResults = s$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 2)
#save(slResults,file="ovarianSmall115.rda")
load("ovarianSmall115.rda")

candidateNames = c("glasso-ebic-0", 
                   "glasso-ebic-0.5", 
                   "glasso-ric", 
                   "hglasso",
                   "mle", 
                   "glasso-stars-0.05",
                   "glasso-stars-0.1",
                   "qgraph-ebic-0", 
                   "qgraph-ebic-0.5")

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
standardize = function(x){return((x-mean(x))/sd(x))}
validDataStd = lapply(validData,function(x){t(apply(x,1,standardize))})

validationLikelihoods = matrix(rep(NA,15*10),ncol=10)

genesToGet = row.names(slResults$optTheta)
for(d in 1:length(validDataStd))
{
  data = data.frame(validDataStd[[d]])
  data$gene = row.names(data)
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
}

scaledVal = data.frame(t(apply(validationLikelihoods,1,function(x){x<- (x-max(x))/max(abs(x))})))
names(scaledVal) = c(candidateNames,"SpiderLearner")
scaledVal$dataset = 1:15

scaledVal %>% na.omit() %>% pivot_longer(cols = names(scaledVal)[1:10])  %>%
  ggplot(aes(x=name,y=100*value)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1,size=12)) +
  xlab("Method") +
  ylab("Percent Difference in Log-Likelihood from Best Model") +
  ggtitle("External Validation on 11 Independent Datasets") +
  geom_jitter(width = 0.1,alpha=0.5,size=1.75)

ggsave("ovarianExternalValidation.pdf")



scaledVal %>% 
  na.omit() %>% pivot_longer(cols = names(scaledVal)[1:10])  %>%
  filter(name %in% c("glasso-ebic-0" ,"hglasso","qgraph-ebic-0","SpiderLearner")) %>%
  ggplot(aes(x=name,y=100*value)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1,size=12)) +
  xlab("Method") +
  ylab("Percent Difference in Log-Likelihood from Best Model") +
  ggtitle("External Validation on 11 Independent Datasets") + 
  geom_jitter(width = 0.1,alpha=0.5,size=1.75)

ggsave("ovarianExternalValidationZoom.pdf")

# # We (theoretically) get a better out-of-sample likelihood
# # Let's try sample splitting as another way to demo
# # the utility of this method
# 
# set.seed(1202)
# trainIndices = sample(rep(1:10,26))
# 
# ensLoss = rep(NA,10)
# candidateLoss = matrix(rep(NA,90),nrow=10,ncol=9)
# for(k in 1:10)
# {
#   print(paste("Working in fold:",k))
#   lateStageTrain = lateStageSmall[trainIndices != k,]
#   lateStageTest = lateStageSmall[trainIndices == k,]
#   slResultTrain = s$runSpiderLearner(lateStageTrain, K = 10, nCores = 1)
#   ensLoss[k] = loglikLossfunction(slResultTrain$optTheta,lateStageTest)
#   candidateLoss[k,] = sapply(slResultTrain$fullModels,loglikLossfunction,lateStageTest)
# }
# 

load("candidateLoss.rda")
load("ensLoss.rda")

# pdf("ovarianLargeOOSL.pdf",width=10,height=6)

df = as.data.frame(cbind(candidateLoss,ensLoss))
colnames(df)[1:9] = candidateNames
colnames(df)[10] = "SpiderLearner"
df$fold = 1:10
df %>% as.data.frame  %>%
  pivot_longer(cols = 1:10)  %>%
  ggplot(aes(x=name,y=value)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1,size=12)) +
  xlab("Method") +
  ylab("Out-of-sample Log Likelihood") + 
  ggtitle("10-fold Internal Cross-Validation on Yoshihara Dataset") +
  geom_jitter(width = 0.1,alpha=0.5,size=1.75)

ggsave("ooslValidation.pdf")

# Look at community membership
load("ovarianSmall115.rda")
adjMat = -cov2cor(slResults$optTheta)
colnames(adjMat)=colnames(lateStageSmall)
adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")

absGraph = adjGraph
E(absGraph)$weight = abs(E(absGraph)$weight)
optThetaComm = cluster_fast_greedy(absGraph) #cluster_walktrap(absGraph,steps = optSteps)

clusterList = list()

clusterList$ensemble = optThetaComm

dummyAdjMat = ifelse(-cov2cor(slResults$fullModels[[6]]) == 0, 0, 1)

ensembleBosses = c()
for(i in unique(factor(optThetaComm$membership)))
{
  thisComm = optThetaComm$names[which(optThetaComm$membership==i)]
  thisSubgraph = induced_subgraph(adjGraph,V(adjGraph)[which(colnames(lateStageSmall) %in% thisComm)])
  commHubs = hub_score(thisSubgraph)$vector
  localBoss = thisComm[which.max(commHubs)]
  ensembleBosses = c(ensembleBosses, localBoss)
}

set.seed(46)

myLayout= layout_with_fr(graph_from_adjacency_matrix(dummyAdjMat,diag=F,weighted=T, mode="undirected"))#, weights=NULL)

jpeg("../Figures/ensembleCommunitiesWithHubs_final.jpeg",width=9,height=9,units="in",res=300)
plot(adjGraph,layout=myLayout,
     vertex.size=ifelse(colnames(lateStageSmall)%in%ensembleBosses,22,6),
     vertex.label=ifelse(colnames(lateStageSmall)%in%ensembleBosses,colnames(lateStageSmall),NA),
     edge.width = 5*E(adjGraph)$weight,
     vertex.label.color="black",
     vertex.label.cex=0.8,
     vertex.label.dist=0,
     vertex.label.degree=pi/2,
     vertex.color=optThetaComm$membership, #ifelse(colnames(lateStageSmall)%in%ensembleBosses,"black",optThetaComm$membership),
     edge.color = "gray50")

dev.off()
