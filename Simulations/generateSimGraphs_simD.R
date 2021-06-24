library(ggplot2)
library(huge)
library(igraph)
library(MASS)
library(moments)
library(pracma)

#### Make Gold Standard Network Structures #### 

makeGoldStandardNets = function(P)
{
  PminusOne = P-1
  
  # make random graph
  erHigh = sample_gnp(P,0.20)
  erHighDensity = length(E(erHigh))/(P*PminusOne/2)
  erLow = sample_gnp(P,0.06)
  erLowDensity = length(E(erLow))/(P*PminusOne/2)
  
  print(erHighDensity)
  print(erLowDensity)
  
  # make small world graph
  wsHigh = sample_smallworld(dim=1,size=P, nei=10, p=0.5)
  wsHighDensity = length(E(wsHigh))/(P*PminusOne/2)
  wsLow = sample_smallworld(dim=1,size=P, nei=3, p=0.5)
  wsLowDensity = length(E(wsLow))/(P*PminusOne/2)
  
  print(wsHighDensity)
  print(wsLowDensity)
  
  # make scale-free graph
  sfHigh = sample_pa(n=P, power=1, m=10, directed=F)
  sfHighDensity = length(E(sfHigh))/(P*PminusOne/2)
  sfLow = sample_pa(n=P, power=1, m=3, directed=F)
  sfLowDensity = length(E(sfLow))/(P*PminusOne/2)
  
  print(sfHighDensity)
  print(sfLowDensity)
  
  # make hub-and-spoke graph
  hsHigh = sample_pa(n=P, power=1.618, m=10, directed=F)
  hsHighDensity = length(E(hsHigh))/(P*PminusOne/2)
  hsLow = sample_pa(n=P, power=1.618, m=3, directed=F)
  hsLowDensity = length(E(sfLow))/(P*PminusOne/2)
  
  print(hsHighDensity)
  print(hsLowDensity)
  
  #### Assign Edge Weights ####
    
  assignRealDataEdgeWeights = function(graph,ehist)
  {
    P = length(V(graph))
    weights = sample(ehist$mids, replace=T,size=length(E(graph)), prob=ehist$density)  
    E(graph)$weights = weights
    precMat = -1*as_adjacency_matrix(graph, attr = c("weights"), type="both") + diag(P)
    return(list(graph, precMat))
  }
  
  realDataHist = read.csv("mxDist.csv") 

  erHighRealData = assignRealDataEdgeWeights(erHigh,ehist=realDataHist)
  erLowRealData = assignRealDataEdgeWeights(erLow,ehist=realDataHist)
  wsHighRealData = assignRealDataEdgeWeights(wsHigh,ehist=realDataHist)
  wsLowRealData = assignRealDataEdgeWeights(wsLow,ehist=realDataHist)
  sfHighRealData = assignRealDataEdgeWeights(sfHigh,ehist=realDataHist)
  sfLowRealData = assignRealDataEdgeWeights(sfLow,ehist=realDataHist)
  hsHighRealData = assignRealDataEdgeWeights(hsHigh,ehist=realDataHist)
  hsLowRealData = assignRealDataEdgeWeights(hsLow,ehist=realDataHist)
  
  #### Adjust Covariance Matrices to be Positive Definite ####
  
  boost = function(myMatrix)
  {
    minEig = min(eigen(myMatrix)$values)
    if(minEig < 0)
    {
      print("Boosting!")
      return(myMatrix - minEig*1.01*diag(ncol(myMatrix)))
    }
    else
      return(myMatrix)
  }
  
  erLowPrecRealDataBoost = boost(erLowRealData[[2]])
  erHighPrecRealDataBoost = boost(erHighRealData[[2]])
  wsLowPrecRealDataBoost = boost(wsLowRealData[[2]])
  wsHighPrecRealDataBoost = boost(wsHighRealData[[2]])
  sfLowPrecRealDataBoost = boost(sfLowRealData[[2]])
  sfHighPrecRealDataBoost = boost(sfHighRealData[[2]])
  hsLowPrecRealDataBoost = boost(hsLowRealData[[2]])
  hsHighPrecRealDataBoost = boost(hsHighRealData[[2]])
  
  
  adjMatListRealData = c(erLowPrecRealDataBoost,
                        erHighPrecRealDataBoost,
                        wsLowPrecRealDataBoost,
                        wsHighPrecRealDataBoost,
                        sfLowPrecRealDataBoost,
                        sfHighPrecRealDataBoost,
                        hsLowPrecRealDataBoost,
                        hsHighPrecRealDataBoost)
  
  
  for(mat in adjMatListRealData)
  {
    print(isposdef(as.matrix(mat)))
  }
  
  return(adjMatListRealData)
}

# graphLabels = c("Random - Low Density",
#                 "Random - High Density",
#                 "Small World - Low Density",
#                 "Small World - High Density",
#                 "Scale-Free - Low Density",
#                 "Scale-Free - High Density",
#                 "Hub-and-Spoke - Low Density",
#                 "Hub-and-Spoke - High Density")
# 
# pdf("graphLayouts.pdf",width=8,height=16)
# par(mfrow=c(4,2),mar=c(1,1,1,1))
# for(m in 1:8)
# {
#   thisGraph = graph_from_adjacency_matrix(-cov2cor(as.matrix(adjMatListRealData[[m]])),
#                                           weighted=T,
#                                           mode="undirected",
#                                           diag=F)
#   dummyGraph = thisGraph
#   E(dummyGraph)$weight = NA #ifelse(E(thisGraph)$weight > 0, 1, 0)
#   myLayout = layout_with_fr(dummyGraph)
#   plot(thisGraph, vertex.label=NA,
#        edge.width = ifelse(abs(E(thisGraph)$weight)>0,1.5,0),
#        vertex.size = (degree(thisGraph)+2),
#        vertex.color = rainbow(50)[degree(thisGraph)],
#        layout=myLayout,
#        main = graphLabels[m])
# }
# dev.off()
