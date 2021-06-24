source("spiderLearnerOOP.R")
library(affy)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("curatedOvarianData")

library(curatedOvarianData)
standardize = function(x){return((x-mean(x))/sd(x))}

data(GSE32062.GPL6480_eset)
lateStage = exprs(GSE32062.GPL6480_eset)

# Extract genes included in the HPO ovarian carcinoma pathway
# See https://hpo.jax.org/app/browse/term/HP:0025318
# Reference: Kohler et al. (2021) The Human Phenotype Ontology in 2021, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D1207–D1217, https://doi.org/10.1093/nar/gkaa1043

lateStageSmall = lateStage[c(1680,1681,3027,4564,8930,12243,12245,13694,13695,13701,13979,16082,16875,17980),]
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

#slResults = s$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 2)

save(slResults,file="ovarianSmall.rda")
load("ovarianSmall.rda")
adjMat = -cov2cor(slResults$optTheta)
adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")
absGraph = adjGraph
E(absGraph)$weight = abs(E(absGraph)$weight)
dummyGraph = absGraph
E(dummyGraph)$weight = NA

jpeg("smallOvarian1.jpeg",width=8,height=8,units="in",res=300)
set.seed(60)
myLayout = layout_with_fr(absGraph)
par(mar=c(3,3,3,3))
plot(adjGraph,layout=myLayout, 
     vertex.shape="rectangle",
     vertex.size=50,
     vertex.size2=20,
     vertex.color="white",
     vertex.label=colnames(lateStageSmall),
     edge.width =30*abs(E(adjGraph)$weight),
     vertex.label.color="black",
     vertex.label.cex=1.2,
     vertex.label.font=2,
     edge.color= ifelse(E(adjGraph)$weight>0,"red","blue")) #ifelse(E(adjGraph)$weight>0,"red","blue"))

title("Ensemble Model",line=-1,cex.main=2.5)
dev.off()

jpeg("smallOvarian2.jpeg",width=15,height=8,units="in",res=300)
par(mfrow=c(2,4),mar=c(0,0,0,0))
for(i in 1:9)
{
  adjMat = -cov2cor(slResults$fullModels[[i]])
  adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")
  if(length(E(adjGraph))> 0)
  {
    plot(adjGraph,layout=myLayout,
         vertex.shape="rectangle",
         vertex.size=50,
         vertex.size2=20,
         vertex.color="white",
         vertex.label=colnames(lateStageSmall),
         edge.width =30*abs(E(adjGraph)$weight),
         vertex.label.color="black",
         vertex.label.cex=1.5,
         vertex.label.font=2,
         edge.color= ifelse(E(adjGraph)$weight>0,"red","blue")) #ifelse(E(adjGraph)$weight>0,"red","blue"))
    title(graphNames[i],line=-2,cex.main=2)
  }
}


adjMat = -cov2cor(slResults$simpleMeanNetwork)
adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")
plot(adjGraph,layout=myLayout,
     vertex.shape="rectangle",
     vertex.size=50,
     vertex.size2=20,
     vertex.color="white",
     vertex.label=colnames(lateStageSmall),
     vertex.label.font=2,
     edge.width =30*abs(E(adjGraph)$weight),
     vertex.label.color="black",
     vertex.label.cex=1.5,
     edge.color= ifelse(E(adjGraph)$weight>0,"red","blue")) #ifelse(E(adjGraph)$weight>0,"red","blue"))
title("Simple Mean",line=-2,cex.main=2)

dev.off()


