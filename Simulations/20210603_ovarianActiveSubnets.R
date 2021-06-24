source("SpiderLearner/SpiderLearner.R")

library(affy)
library(ggplot2)
library(survival)
library(survMisc)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("curatedOvarianData")

library(curatedOvarianData)
standardize = function(x){return((x-mean(x))/sd(x))}

# GSE32062.GPL6480_eset

data(GSE32062.GPL6480_eset)
lateStage = exprs(GSE32062.GPL6480_eset)
lateStageClinical = pData(GSE32062.GPL6480_eset)
lateStageClinical$vital_status = ifelse(lateStageClinical$vital_status == "deceased",1,0)
lateStageClinical$recurrence_status = ifelse(lateStageClinical$recurrence_status == "recurrence",1,0)
yoshi = read.table("YoshiharaGeneSet.tsv",sep="&")
lateStageSmall = lateStage[which(rownames(lateStage)%in%yoshi[,1]),]
lateStageSmall = data.frame(t(lateStageSmall))
names(lateStageSmall) = colnames(lateStageSmall)
#lateStageSmall = lateStageSmall[which(lateStageClinical$summarygrade == "high"),]
#lateStageClinical = lateStageClinical[which(lateStageClinical$summarygrade == "high"),]
lateStageSmall = apply(lateStageSmall,2,standardize)
colnames(lateStageSmall)[which(colnames(lateStageSmall)=="HLA-DPB1")] = "HLA.DPB1"
names(lateStageSmall)=colnames(lateStageSmall)

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


set.seed(543)
slResults = s$runSpiderLearner(lateStageSmall, K = 10, standardize=T,boundedLoss = F)

# load("ovarian.rda") # result from the line above

graphNames = c("glasso - ebic - 0", 
               "glasso - ebic - 0.5",
               "glasso - ric", 
               "hglasso",
               "MLE",
               "glasso - stars - 0.05",
               "glasso - stars - 0.1",
               "qgraph - ebic - 0",
               "qgraph - ebic - 0.5")


adjGraph = graph_from_adjacency_matrix(-cov2cor(slResults$optTheta),weighted=T,mode="undirected",diag=F)
hub_ensemble = hub_score(adjGraph)$vector
hub_candidates = data.frame("gene"=colnames(lateStageSmall),"ensemble"=hub_ensemble)
row.names(hub_candidates)=hub_candidates$gene

for(i in 1:9)
{
  adjGraph = graph_from_adjacency_matrix(-cov2cor(slResults$fullModels[[i]]),weighted=T,mode="undirected",diag=F)
  V(adjGraph)$labels = colnames(lateStageSmall)
  hub_model_1 = hub_score(adjGraph)$vector
  hub_candidates = cbind(hub_candidates,hub_model_1)
  
}

hub_candidates = hub_candidates[order(-hub_candidates$ensemble),]

# What are the hubs in the communities?
## cluster analysis ##

findMaxModStep = function(myGraph,consider=1:10)
{
  # mygraph must have positive edges only
  maxModularity = -1
  besti = 0
  modularities = c()
  
  for(i in consider)
  {
    optThetaComm = cluster_walktrap(myGraph,steps = i)
    thisModularity = modularity(optThetaComm)
    modularities = c(modularities, thisModularity)
    if(thisModularity > maxModularity)
    {
      besti=i
      maxModularity = thisModularity
    }
  }
  return(list("maxMod"=maxModularity,"best_i"=besti,"allMod"=modularities,"all_i"=consider))
}


adjMat = -cov2cor(slResults$optTheta)
# Necessary to correct for the HLA- issue
colnames(adjMat)=colnames(lateStageSmall)
adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")

absGraph = adjGraph
E(absGraph)$weight = abs(E(absGraph)$weight)
optSteps = findMaxModStep(absGraph)$best_i
optThetaComm = cluster_walktrap(absGraph,steps = optSteps)

clusterList = list()

clusterList$ensemble = optThetaComm

dummyAdjMat = ifelse(-cov2cor(slResults$fullModels[[6]]) == 0, 0, 1)

set.seed(60)
myLayout= layout_with_fr(graph_from_adjacency_matrix(dummyAdjMat,diag=F,weighted=T, mode="undirected"))#, weights=NULL)

jpeg("../Figures/ensemble.jpeg",width=8,height=8,units="in",res=300)
plot(optThetaComm,
     adjGraph,layout=myLayout, 
     vertex.size=6,
     vertex.label=NA, 
     edge.width = 5*E(adjGraph)$weight,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     edge.color = "gray50")
dev.off()

ensembleBosses = c()
for(i in unique(factor(optThetaComm$membership)))
{
  thisComm = optThetaComm$names[which(optThetaComm$membership==i)]
  thisSubgraph = induced_subgraph(adjGraph,V(adjGraph)[which(colnames(lateStageSmall) %in% thisComm)])
  commHubs = hub_score(thisSubgraph)$vector
  localBoss = thisComm[which.max(commHubs)]
  ensembleBosses = c(ensembleBosses, localBoss)
}

jpeg("../Figures/ensembleCommunitiesWithHubs.jpeg",width=8,height=8,units="in",res=300)
plot(adjGraph,layout=myLayout, 
     vertex.size=ifelse(colnames(lateStageSmall)%in%ensembleBosses,22,6),
     #mark.groups=communities(optThetaComm),
     vertex.label=ifelse(colnames(lateStageSmall)%in%ensembleBosses,colnames(lateStageSmall),NA), 
     edge.width = 5*E(adjGraph)$weight,
     vertex.label.color="black",
     vertex.label.cex=0.8,
     vertex.label.dist=0,
     vertex.label.degree=pi/2,
     vertex.color=optThetaComm$membership, #ifelse(colnames(lateStageSmall)%in%ensembleBosses,"black",optThetaComm$membership),
     edge.color = "gray50",
     main = "ensemble")
dev.off()

jpeg("../Figures/allComms.jpeg",width=24,height=24,units="in",res=200)
par(mfrow=c(3,3))
localBossList = list()
for(i in 1:9)
{
  adjMat = -cov2cor(slResults$fullModels[[i]])
  # Necessary to correct for the HLA- issue
  colnames(adjMat)=colnames(lateStageSmall)
  adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")
  absGraph = adjGraph
  if(length(E(absGraph))>0) E(absGraph)$weight = abs(E(absGraph)$weight)
  nSteps = findMaxModStep(absGraph)$best_i
  clusters = cluster_walktrap(absGraph,steps = nSteps)
  clusterList[[i+1]] = clusters
  
  localBosses = c()
  for(j in unique(clusters$membership))
  {
    thisComm = colnames(lateStageSmall)[which(clusters$membership==j)]
    thisSubgraph = induced_subgraph(adjGraph,V(adjGraph)[which(colnames(lateStageSmall) %in% thisComm)])
    commHubs = hub_score(thisSubgraph)$vector
    localBoss = thisComm[which.max(commHubs)]
    localBosses = c(localBosses, localBoss)
  }
  
  localBossList[[i]] = localBosses
  plot(clusters,adjGraph,layout=myLayout, 
       vertex.size=6, #ifelse(colnames(lateStageSmall)%in%localBosses,22,6),
       #mark.groups=communities(thisComm),
       vertex.label=NA, #ifelse(colnames(lateStageSmall)%in%localBosses,colnames(lateStageSmall),NA), 
       edge.width = 5*E(adjGraph)$weight,
       vertex.label.color="black",
       vertex.label.cex=1,
       vertex.label.dist=0,
       vertex.label.degree=pi/2,
       vertex.color=clusters$membership, #ifelse(colnames(lateStageSmall)%in%ensembleBosses,"black",optThetaComm$membership),
       edge.color = "gray50")
  title(graphNames[i],cex.main=5)
  
}
dev.off()


# are the community hubs associated with outcomes
lateStageSmallPlus = data.frame(lateStageSmall)
lateStageSmallPlus$days_to_death = c(lateStageClinical$days_to_death)
lateStageSmallPlus$vital_status = lateStageClinical$vital_status
  
myFormula = as.formula(paste("Surv(days_to_death, vital_status) ~ ",paste(ensembleBosses, collapse= "+")))
myMod2 = coxph(myFormula, data = lateStageSmallPlus)

#x = lateStageSmall[,which(colnames(lateStageSmall) %in% ensembleBosses)]
#cvfit <- cv.glmnet(as.matrix(x), Surv(lateStageClinical$days_to_death, lateStageClinical$vital_status), family = "cox", type.measure = "C")
#fit <- glmnet(x, Surv(lateStageClinical$days_to_death, lateStageClinical$vital_status), lambda=cvfit$lambda.min, family = "cox")
#commBetas = as.matrix(fit$beta)
#commBetasList[[j+1]] = myMod2$commBetas[,1]

commBetasList = list()
commBetasList[[1]] = myMod2$coefficients

for(j in c(1:9))
{
  print(paste("Processing method :",j))
  commBetasList[[j+1]] = list()
  
  for(i in sort(unique(clusterList[[(j+1)]]$membership)))
  {
    myFormula = as.formula(paste("Surv(days_to_death, vital_status) ~ ",
                                 paste(localBossList[[j]], collapse= "+")))
    myMod2 = coxph(myFormula, data = lateStageSmallPlus)
    commBetasList[[j+1]] = myMod2$coefficients
  }
}

getGeneralScoreContinuous = function(data,commBetasList)
{
  commScoreTotal = rep(0,nrow(data))
  geneMatches = intersect(names(data),names(commBetasList))
  # gene matches will be alphabetical
  commExtraction = data[,which(names(data) %in% geneMatches)]
  # commExtraction will be alphabetical
  theseBetas = commBetasList[which(names(commBetasList)%in% geneMatches)]
  theseBetas = theseBetas[order(names(theseBetas))]
  commScore = as.matrix(commExtraction) %*% as.matrix(theseBetas)
  return(commScore)
}

getGeneralScore = function(data,commBetasList,cutoff)
{
  commScoreTotal = rep(0,nrow(data))
  geneMatches = intersect(names(data),names(commBetasList))
  # gene matches will be alphabetical
  commExtraction = data[,which(names(data) %in% geneMatches)]
  # commExtraction will be alphabetical
  theseBetas = commBetasList[which(names(commBetasList)%in% geneMatches)]
  theseBetas = theseBetas[order(names(theseBetas))]
  commScore = as.matrix(commExtraction) %*% as.matrix(theseBetas)
  commScoreTotal = commScoreTotal+commScore #Binary
  commScoreBinary = ifelse(commScoreTotal > cutoff,1,0)
  return(commScoreBinary)
}


# Get scores for training data to determine optimal cutoff

overallCutoff = list()
for(j in c(1:10))
{
  scores = getGeneralScoreContinuous(lateStageSmallPlus,commBetasList[[j]])
  seqMin = min(scores)
  seqMax = max(scores)
  survDiffsTrain = c()
  for(cutoff in seq(seqMin,seqMax,by=0.01))
  {
    scTrain = getGeneralScore(lateStageSmallPlus, commBetasList[[j]],cutoff)
    if(sum(scTrain==0)>2 & sum(scTrain==1)>2)
    {
      xyz = survdiff(Surv(lateStageSmallPlus$days_to_death, lateStageSmallPlus$vital_status)~scTrain)
      survDiffsTrain = c(survDiffsTrain,1 - pchisq(xyz$chisq, length(xyz$n) - 1))
    }
    else
      survDiffsTrain = c(survDiffsTrain,NA)
  }

  overallCutoff[[j]] = seq(seqMin,seqMax,by=0.01)[which.min(survDiffsTrain)]

}



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

phenoData = list(pData(E.MTAB.386_eset),
                 pData(GSE13876_eset),
                 pData(GSE14764_eset),
                 pData(GSE17260_eset),
                 pData(GSE18520_eset),
                 pData(GSE19829.GPL570_eset),
                 pData(GSE19829.GPL8300_eset),
                 pData(GSE26712_eset),
                 pData(GSE30009_eset),
                 pData(GSE30161_eset),
                 pData(GSE32063_eset),
                 pData(GSE9891_eset),
                 pData(PMID17290060_eset),
                 pData(PMID19318476_eset),
                 pData(TCGA_eset))

survDiffsTestMat = matrix(rep(NA,length(phenoData)*10),nrow=length(phenoData),byrow=T)

for(id in c(1:(length(phenoData))))
{
  # preprocess
  ocValidation = t(validData[[id]])
  ocPheno = phenoData[[id]]
  ocPheno$vital_status = ifelse(ocPheno$vital_status == "deceased",1,0)
  ocValidationClean = apply(ocValidation,2,standardize)
  colnames(ocValidationClean)[which(colnames(ocValidationClean)=="HLA-DPB1")] = "HLA.DPB1"
  names(ocValidationClean) = colnames(ocValidationClean)
  
  survDiffsTrain = c()
  survDiffsTest = c()
  
  for(j in c(1:10))
  {
    sc = getGeneralScore(ocValidationClean, commBetasList[[j]], overallCutoff[[j]])
    if(sum(sc==0)>2 & sum(sc==1)>2) # exclude cases where everyone gets called the same
    { # need to also exclude cases where there is only one person different
      plot(survfit(Surv(ocPheno$days_to_death, ocPheno$vital_status)~sc), col=1:5,main="Out Of Sample")
      xyz = survdiff(Surv(ocPheno$days_to_death, ocPheno$vital_status)~sc)
      survDiffsTest = c(survDiffsTest,1 - pchisq(xyz$chisq, length(xyz$n) - 1))
    }
    else
      survDiffsTest = c(survDiffsTest,NA)
  }
  
  survDiffsTestMat[id,]=survDiffsTest
  
}

pdf("../Figures/validation_log_rank_tests.pdf",width=10,height=6)
par(mar=c(5,20,2,2))
boxplot(survDiffsTestMat,las=2,horizontal=T,
        main="Validation Performance",
        col=rainbow(10),
        xlab="p-value of log rank test between high and low risk groups",
        names=paste(paste(c("SpiderLearner",graphNames)," - Number of Genes:"),sapply(commBetasList,length)))
dev.off()

apply(survDiffsTestMat,2,function(x){return(sum(x<0.05,na.rm=T))})

# which studies validated?
apply(survDiffsTestMat,2,function(x){return(which(x<0.05))})

tableSource = data.frame("Method"=c("ensemble",graphNames))
tableSource$Cutoff = round(unlist(overallCutoff),3)
tableSource$NValPoss = apply(survDiffsTestMat,2,function(x){sum(!is.na(x))})
tableSource$NVal = round(apply(survDiffsTestMat,2,function(x){sum(x<0.05,na.rm=T)}),2)
tableSource$MeanP = round(apply(survDiffsTestMat,2,function(x){mean(x,na.rm=T)}),2)
tableSource$MedianP = round(apply(survDiffsTestMat,2,function(x){median(x,na.rm=T)}),2)

write.table(tableSource,file="validationTable.asv",quote=F,row.names=F, sep="&")

# What is the median survival time in the high-risk and low-risk groups
# in each of the validated datasets?

validatedMedians = matrix(rep(NA,2*(length(phenoData))),ncol=2)

pdf("../Figures/validatedKMs.pdf",height=8,width=10) # Uncomment this and next line for validated only
par(mfrow=c(2,3))
#pdf("../Figures/allValidationKMs.pdf",height=20,width=10)
#par(mfrow=c(5,3))
for(id in c(1:(length(phenoData))))
{
  # preprocess
  ocValidation = t(validData[[id]])
  ocPheno = phenoData[[id]]
  ocPheno$vital_status = ifelse(ocPheno$vital_status == "deceased",1,0)
  ocValidationClean = apply(ocValidation,2,standardize)
  colnames(ocValidationClean)[which(colnames(ocValidationClean)=="HLA-DPB1")] = "HLA.DPB1"
  names(ocValidationClean) = colnames(ocValidationClean)
  
  for(j in c(1))
  {
    sc = getGeneralScore(ocValidationClean, commBetasList[[j]], overallCutoff[[j]])
    if(sum(sc==0)>2 & sum(sc==1)>2) # exclude cases where everyone gets called the same
    { 
      myFit = survfit(Surv(ocPheno$days_to_death, ocPheno$vital_status)~sc)
      xyz = survdiff(Surv(ocPheno$days_to_death, ocPheno$vital_status)~sc)
      pVal = 1 - pchisq(xyz$chisq, length(xyz$n) - 1)
      if (pVal<0.05) #uncomment this for validated only
      {
        plot(myFit, main=identifiers[id],xlim=c(0,4000), xlab="Days",ylab="Survival Proportion",lwd=3,col=c("orange","deepskyblue"),pch=c(20,20))
        text(3250,0.8,paste("log-rank p=",round(pVal,4)))
        legend(2500,1.0,c("Risk Score=0","Risk Score=1"),col=c("orange","deepskyblue"),pch=c(20,20))
        validatedMedians[id,]= summary(myFit)$table[,'median']
      }
    }
    # uncomment these to make the full plot
    # else
    # {
    #   # just don't do a survdiff
    #   myFit = survfit(Surv(ocPheno$days_to_death, ocPheno$vital_status)~sc)
    #   plot(myFit, col=1:5,main=identifiers[id],xlim=c(0,4000), conf.int = FALSE, xlab="Days",ylab="Survival Proportion",lwd=3)
    #   text(3250,0.8,paste("No Comparison"))
    #   legend(2500,1.0,c("Risk Score=0","Risk Score=1"),col=c("orange","deepskyblue"),pch=c(20,20))
    # }
  
  }

}
dev.off()

validatedMedians = as.data.frame(validatedMedians)
names(validatedMedians)=c("Low","High")
row.names(validatedMedians)=identifiers
validatedShortMedians = validatedMedians[c(2,3,6,12,13,14),]
write.table(validatedShortMedians,file="validatedShortMedians.asv",quote=F,sep="&",row.names=T)

# 20210520 Adding CoxPH analysis for baseline, using available covariates
# from dataset

# Determine what covariates might be good 
for(id in c(1:(length(phenoData))))
{
  print(identifiers[[id]])
  print("age:")
  print(summary(phenoData[[id]]$age_at_initial_pathologic_diagnosis))
  print("summary grade:")
  
  print(summary(factor(phenoData[[id]]$summarygrade)))
  print("summary stage:")
  
  print(summary(factor(phenoData[[id]]$summarystage)))
  print("tumor stage:")
  print(summary(factor(phenoData[[id]]$tumorstage)))
  
  print("grade:")
  print(summary(factor(phenoData[[id]]$grade)))
  print("--------")
}

unadjustedCoefs = matrix(rep(NA,9*15),ncol=9)
adjustedCoefs = matrix(rep(NA,9*15),ncol=9)

for(id in c(1:(length(phenoData))))
{
  # preprocess data (can do this more efficiently later)
  ocValidation = t(validData[[id]])
  ocPheno = phenoData[[id]]
  ocPheno$vital_status = ifelse(ocPheno$vital_status == "deceased",1,0)
  ocValidationClean = apply(ocValidation,2,standardize)
  colnames(ocValidationClean)[which(colnames(ocValidationClean)=="HLA-DPB1")] = "HLA.DPB1"
  names(ocValidationClean) = colnames(ocValidationClean)
  
  # fit the mod with the ensemble risk score, continuous
  sc = getGeneralScoreContinuous(ocValidationClean, commBetasList[[1]])
  thisMod0 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc)
  if(id == 1 | id == 8)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$tumorstage < 4))
  if(id ==2)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$summarygrade))
  if(id ==3)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + factor(ocPheno$summarygrade) + factor(ocPheno$summarystage))
  if(id == 4 | id == 11 | id == 13)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + factor(ocPheno$summarygrade) + factor(ocPheno$tumorstage < 4))
  if(id == 8)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$tumorstage < 4))
  if(id == 9 | id == 10 | id == 14)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$summarygrade) + factor(ocPheno$tumorstage < 4))
  if(id == 12 | id == 15)
    thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$summarygrade) +  factor(ocPheno$summarystage))
  
  unadjusted = summary(thisMod0)$coef
  unadjustedCoefs[id,] = c(unadjusted,summary(thisMod0)$conf.int)
  
  if(!(id %in% c(5,6,7)))
  {
    adjusted = summary(thisMod2)$coef[1,]
    adjustedCoefs[id,] = c(adjusted,summary(thisMod2)$conf.int[1,])
  }
}

# Forest plot of unadjusted and adjusted coefficients in each model
unadjustedCoefs = data.frame(unadjustedCoefs)
adjustedCoefs = data.frame(adjustedCoefs)
names(unadjustedCoefs) = colnames(summary(thisMod0)$coef)
names(adjustedCoefs) = colnames(summary(thisMod0)$coef)
names(unadjustedCoefs)[3] = "SE"
names(adjustedCoefs)[3] = "SE"

fpdf = rbind.data.frame(unadjustedCoefs,adjustedCoefs)
fpdf$model = c(rep("unadjusted",15),rep("adjusted",15))
fpdf$dataset = identifiers
fpdf = na.omit(fpdf)
names(fpdf)[5]="p"
fpdf$dataset = factor(fpdf$dataset, levels=identifiers[order(unadjustedCoefs$coef)])
forest = ggplot(data=fpdf, aes(x=dataset, y=coef, ymin=coef-1.96*SE, 
                                 ymax=coef+1.96*SE,col = model,
                               shape = p<0.05)) + 
  geom_pointrange(lwd=0.8, position = position_jitterdodge(jitter.width=0,dodge.width = 0.7)) +
  geom_hline(aes(yintercept = 0), lty=2) +
  ggtitle("Association of Continuous Risk Score \n with Time to Death") +
  xlab("Model") + ylab("log(Hazard Ratio)") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        title = element_text(face = "bold", size = 12),
        strip.text.x = element_blank(),#element_text(face = "bold", size = 12),
        strip.text.y = element_blank(),
        legend.text = element_text(face = "bold", size = 12),
        legend.position="right",
        legend.box = 'vertical') 

ggsave(forest,filename = "../Figures/coxPHEnsembleScores.jpeg")

unadjustedCoxPHPvals = unadjustedCoefs[,5]
unadjustedKMPvals = survDiffsTestMat[,1]

#pComp = data.frame("KM"= survDiffsTestMat[,1], 
#           "UnadjCoxPH"=unadjustedCoefs[,5], 
#           "AdjCoxPH"=adjustedCoefs[,5])

pComp = data.frame(#"KM"= survDiffsTestMat[,1],
  "UnadjHR"=round(unadjustedCoefs[,6],2),
  "UnadjHRLower"=round(unadjustedCoefs[,8],2),
  "UnadjHRUpper"=round(unadjustedCoefs[,9],2),
  "UnadjCoxPH"=round(unadjustedCoefs[,5],2),
  "AdjHR"=round(adjustedCoefs[,6],2),
  "AdjHRLower"=round(adjustedCoefs[,8],2),
  "AdjHRUpper"=round(adjustedCoefs[,9],2),
  "AdjCoxPH"=round(adjustedCoefs[,5],2))

pComp$Covariates = NA
pComp$Covariates[c(1,8)] = "Age, Tumor Stage < 4"
pComp$Covariates[2] = "Age, Summary Grade"
pComp$Covariates[3] = "Summary Grade, Summary Stage"
pComp$Covariates[c(4,11,13)] = "Summary Grade, Tumor Stage < 4"
pComp$Covariates[c(9,10,14)] = "Age, Summary Grade, Tumor Stage < 4"
pComp$Covariates[c(12,15)] = "Age, Summary Grade, Summary Stage"

row.names(pComp)=identifiers
write.table(pComp,file="pValueComparisions.txt",sep="\t")
write.table(pComp,file="pValueComparisions.asv",sep=" & ",quote=F)

row.names(pComp)=identifiers
write.table(round(pComp,3),file="pValueComparisions.txt",sep="\t")

# Get the coxph vals for risk scores from simple mean and all 9 candidates
# for comparison

megaUnadjustedCoefs = list()
megaAdjustedCoefs = list()

for(j in 1:10)
{
  megaUnadjustedCoefs[[j]] =  matrix(rep(NA,5*15),ncol=5)
  megaAdjustedCoefs[[j]] = matrix(rep(NA,5*15),ncol=5)
    
  for(id in c(1:(length(phenoData))))
  {
    # preprocess data (can do this more efficiently later)
    ocValidation = t(validData[[id]])
    ocPheno = phenoData[[id]]
    ocPheno$vital_status = ifelse(ocPheno$vital_status == "deceased",1,0)
    ocValidationClean = apply(ocValidation,2,standardize)
    colnames(ocValidationClean)[which(colnames(ocValidationClean)=="HLA-DPB1")] = "HLA.DPB1"
    names(ocValidationClean) = colnames(ocValidationClean)
    
    # fit the mod with the ensemble risk score, continuous
    sc = getGeneralScoreContinuous(ocValidationClean, commBetasList[[j]])
    thisMod0 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc)
    if(id == 1 | id == 8)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$tumorstage < 4))
    if(id ==2)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$summarygrade))
    if(id ==3)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + factor(ocPheno$summarygrade) + factor(ocPheno$summarystage))
    if(id == 4 | id == 11 | id == 13)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + factor(ocPheno$summarygrade) + factor(ocPheno$tumorstage < 4))
    if(id == 8)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$tumorstage < 4))
    if(id == 9 | id == 10 | id == 14)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$summarygrade) + factor(ocPheno$tumorstage < 4))
    if(id == 12 | id == 15)
      thisMod2 = coxph(Surv(ocPheno$days_to_death,ocPheno$vital_status) ~ sc + ocPheno$age_at_initial_pathologic_diagnosis + factor(ocPheno$summarygrade) +  factor(ocPheno$summarystage))
    
    unadjusted = summary(thisMod0)$coef
    megaUnadjustedCoefs[[j]][id,] = unadjusted
    
    if(!(id %in% c(5,6,7)))
    {
      adjusted = summary(thisMod2)$coef[1,]
      megaAdjustedCoefs[[j]][id,] = adjusted
    }
  }
}

allUnadjustedP = data.frame(lapply(megaUnadjustedCoefs,function(x){x[,5]}))

pdf("../Figures/validation_coxph_tests.pdf",width=10,height=6)
par(mar=c(5,20,2,2))
boxplot(allUnadjustedP,las=2,horizontal=T,
        main="Validation Performance - CoxPH",
        col=rainbow(10),
        xlab="p-value of score coefficient",
        names=paste(paste(c("SpiderLearner",graphNames)," - Number of Genes:"),sapply(commBetasList,length)))
dev.off()
