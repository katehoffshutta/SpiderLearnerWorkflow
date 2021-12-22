library(huge)
library(igraph)
library(moments)
library(pracma)

standardize = function(x) return((x-mean(x))/sd(x))

#### Get discrete uniform distribution from CATHGEN Example ####

cathgenHist = read.csv("mxDist.csv")

#### Make Gold Standard Network Structures #### 

makeGoldStandardNets = function(P)
{
  PminusOne = P-1
  
  # make random graph
  erHigh = sample_gnp(P,0.2)
  erHighDensity = length(E(erHigh))/(P*PminusOne/2)
  erLow = sample_gnp(P,0.11)
  erLowDensity = length(E(erLow))/(P*PminusOne/2)
  print(erLowDensity)

  # make small world graph
  wsHigh = sample_smallworld(dim=1,size=P, nei=10, p=0.5)
  wsHighDensity = length(E(wsHigh))/(P*PminusOne/2)
  wsLow = sample_smallworld(dim=1,size=P, nei=3, p=0.5)
  wsLowDensity = length(E(wsLow))/(P*PminusOne/2)
  print(wsLowDensity)
  
  # make scale-free graph
  sfHigh = sample_pa(n=P, power=1, m=10, directed=F)
  sfHighDensity = length(E(sfHigh))/(P*PminusOne/2)
  sfLow = sample_pa(n=P, power=1, m=3, directed=F)
  sfLowDensity = length(E(sfLow))/(P*PminusOne/2)
  print(sfLowDensity)
  
  # make hub-and-spoke graph
  hsHigh = sample_pa(n=P, power=1.618, m=10, directed=F)
  hsHighDensity = length(E(hsHigh))/(P*PminusOne/2)
  hsLow = sample_pa(n=P, power=1.618, m=3, directed=F)
  hsLowDensity = length(E(sfLow))/(P*PminusOne/2)
  print(hsLowDensity)
  
  #### Assign Edge Weights ####
  
  assignCathgenEdgeWeights = function(graph,ehist)
  {
    P = length(V(graph))
    weights = sample(ehist$mids, replace=T,size=length(E(graph)), prob=ehist$density)  
    E(graph)$weights = weights
    precMat = -1*as_adjacency_matrix(graph, attr = c("weights"), type="both") + diag(P)
    return(list(graph, precMat))
  }
  
  erHighCathgen = assignCathgenEdgeWeights(erHigh,ehist=cathgenHist)
  erLowCathgen = assignCathgenEdgeWeights(erLow,ehist=cathgenHist)
  wsHighCathgen = assignCathgenEdgeWeights(wsHigh,ehist=cathgenHist)
  wsLowCathgen = assignCathgenEdgeWeights(wsLow,ehist=cathgenHist)
  sfHighCathgen = assignCathgenEdgeWeights(sfHigh,ehist=cathgenHist)
  sfLowCathgen = assignCathgenEdgeWeights(sfLow,ehist=cathgenHist)
  hsHighCathgen = assignCathgenEdgeWeights(hsHigh,ehist=cathgenHist)
  hsLowCathgen = assignCathgenEdgeWeights(hsLow,ehist=cathgenHist)
  
  #### Adjust Covariance Matrices to be Positive Definite ####
  
  boost = function(myMatrix)
  {
    minEig = min(eigen(myMatrix)$values)
    if(minEig < 0)
    {
      print("Boosting Needed!")
      return(myMatrix - minEig*1.01*diag(ncol(myMatrix)))
    }
    else
      return(myMatrix)
  }
  
  erLowPrecCathgenBoost = boost(erLowCathgen[[2]])
  erHighPrecCathgenBoost = boost(erHighCathgen[[2]])
  wsLowPrecCathgenBoost = boost(wsLowCathgen[[2]])
  wsHighPrecCathgenBoost = boost(wsHighCathgen[[2]])
  sfLowPrecCathgenBoost = boost(sfLowCathgen[[2]])
  sfHighPrecCathgenBoost = boost(sfHighCathgen[[2]])
  hsLowPrecCathgenBoost = boost(hsLowCathgen[[2]])
  hsHighPrecCathgenBoost = boost(hsHighCathgen[[2]])
  
  adjMatListCathgen = c(erLowPrecCathgenBoost,
                        erHighPrecCathgenBoost,
                        wsLowPrecCathgenBoost,
                        wsHighPrecCathgenBoost,
                        sfLowPrecCathgenBoost,
                        sfHighPrecCathgenBoost,
                        hsLowPrecCathgenBoost,
                        hsHighPrecCathgenBoost)
  
  for(mat in adjMatListCathgen)
  {
    print(isposdef(as.matrix(mat)))
  }
  
  return(adjMatListCathgen)
}

set.seed(42)
adjMatListCathgen = makeGoldStandardNets(50)

graphLabels = c("Random", # - Low Density",
                "Random - High Density",
                "Small World", # - Low Density",
                "Small World - High Density",
                "Scale-Free", # - Low Density",
                "Scale-Free - High Density",
                "Hub-and-Spoke", # - Low Density",
                "Hub-and-Spoke - High Density")

#pdf("graphLayouts.pdf",width=16,height=8)
pdf("graphLayouts_Small.pdf",width=24,height=8)
par(mfrow=c(1,4),mar=c(2,2,2,2))
for(m in c(1,3,5,7)) #,2,4,6,8))
{
  thisGraph = graph_from_adjacency_matrix(-cov2cor(as.matrix(adjMatListCathgen[[m]])),
                                          weighted=T,
                                          mode="undirected",
                                          diag=F)
  dummyGraph = thisGraph
  E(dummyGraph)$weight = NA #ifelse(E(thisGraph)$weight > 0, 1, 0)
  myLayout = layout_with_fr(dummyGraph)
  plot(thisGraph, vertex.label=NA,
       edge.width = ifelse(abs(E(thisGraph)$weight)>0,1.5,0),
       vertex.size = 1.5*(degree(thisGraph)),
       vertex.color = rainbow(40)[degree(thisGraph)],
       layout=myLayout)
       #main = graphLabels[m])
  title(graphLabels[m],line=-3,cex.main=3.5)
}
dev.off()
