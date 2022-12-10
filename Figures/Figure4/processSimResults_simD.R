## process simulation results

library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

setwd("~/research/repos/SpiderLearnerWorkflow/Figures/Figure4")

source("../../Simulations/errorMetrics.R")

load("../../Results/SimD/20211201_erHighPrec_RealData_n_60_p_100_simStudy.rda")
erHighPrecResults = theseResults
load("../../Results/SimD/20211201_erLowPrec_RealData_n_60_p_100_simStudy.rda")
erLowPrecResults = theseResults
load("../../Results/SimD/20211203_wsHighPrec_RealData_n_60_p_100_simStudy.rda")
wsHighPrecResults = theseResults
load("../../Results/SimD/20211201_wsLowPrec_RealData_n_60_p_100_simStudy.rda")
wsLowPrecResults = theseResults
load("../../Results/SimD/20211201_sfHighPrec_RealData_n_60_p_100_simStudy.rda")
sfHighPrecResults = theseResults
load("../../Results/SimD/20211201_sfLowPrec_RealData_n_60_p_100_simStudy.rda")
sfLowPrecResults = theseResults
load("../../Results/SimD/20211203_hsHighPrec_RealData_n_60_p_100_simStudy.rda")
hsHighPrecResults = theseResults
load("../../Results/SimD/20211201_hsLowPrec_RealData_n_60_p_100_simStudy.rda")
hsLowPrecResults = theseResults

allResults = list(erLowPrecResults,
                  erHighPrecResults,
                  wsLowPrecResults,
                  wsHighPrecResults,
                  sfLowPrecResults,
                  sfHighPrecResults,
                  hsLowPrecResults,
                  hsHighPrecResults)

nPred = 100
nTopology = 8
load("../../Results/eightNetworksD.rda")
eightNetworks = d

methods = c("glasso - ebic - 0", 
            "glasso - ebic - 0.5", 
                      "glasso - ric",
                      "hglasso",
                      #"mle",
                      "glasso - stars - 0.05",
                      "glasso - stars - 0.1",
                      "qgraph - ebic - 0",
                      "qgraph - ebic - 0.5",
                      "SpiderLearner",
                      "simple mean")

nMod = length(methods)-2
  
allWeights = list()

for(i in 1:length(allResults))
{
  allWeights[[i]] = t(data.frame(sapply(allResults[[i]][[1]],function(x)x$weights)))
}         

### par(mfrow=c(4,2))

### boxplot(allWeights[[1]],main="Random Graph Low Density Ensemble Weights",names = methods[1:nMod],
###         las=2)
### boxplot(allWeights[[2]],main="Random Graph High Density Ensemble Weights",names = methods[1:nMod],las=2)

### boxplot(allWeights[[3]],main="Small World Graph Low Density Ensemble Weights",names = methods[1:nMod],
###         las=2)
### boxplot(allWeights[[4]],main="Small World Graph High Density Ensemble Weights",names = methods[1:nMod],las=2)

### boxplot(allWeights[[5]],main="Scale Free Graph Low Density Ensemble Weights",names = methods[1:nMod],
###         las=2)
### boxplot(allWeights[[6]],main="Scale Free Graph High Density Ensemble Weights",names = methods[1:nMod],las=2)

### boxplot(allWeights[[7]],main="Hub-and-Spoke Graph Low Density Ensemble Weights",names = methods[1:nMod],
###       las=2)
### boxplot(allWeights[[8]],main="Hub-and-Spoke Graph High Density Ensemble Weights",names = methods[1:nMod],las=2)

weightsTable = data.frame("erLow"=apply(allWeights[[1]],2,mean),
                          "erHigh"=apply(allWeights[[2]],2,mean),
                          "wsLow"=apply(allWeights[[3]],2,mean),
                          "wsHigh"=apply(allWeights[[4]],2,mean),
                          "sfLow"=apply(allWeights[[5]],2,mean),
                          "sfHigh"=apply(allWeights[[6]],2,mean),
                          "hsLow"=apply(allWeights[[7]],2,mean),
                          "hsHigh"=apply(allWeights[[8]],2,mean))

row.names(weightsTable) = methods[1:nMod]
write.table(round(t(weightsTable),2),file="../../Tables//weights_simD.asv",sep="&",row.names=T,quote=F)


nSim=100
rfnLong = data.frame("topology"=c(rep("erdos-renyi",(nMod+2)*2*nSim), 
                                  rep("small world",(nMod+2)*2*nSim),
                                  rep("scale free",(nMod+2)*2*nSim),
                                  rep("hub and spoke",(nMod+2)*2*nSim)))
rfnLong$density = rep(c(rep("low density",(nMod+2)*nSim),rep("high density",(nMod+2)*nSim)),4)
rfnLong$method = rep(c(rep("huge-ebic-0",nSim),
                   rep("huge-ebic-0.5",nSim),
                   rep("glasso-ric",nSim),
                   rep("hub glasso",nSim),
                   #rep("mle",nSim),
                   rep("stars-0.05",nSim),
                   rep("stars-0.1",nSim),
                   rep("qgraph-ebic-0",nSim),
                   rep("qgraph-ebic-0.5",nSim),
                   rep("SpiderLearner",nSim),
                   rep("simple mean",nSim)),8)

rfnLong$density = factor(rfnLong$density, levels=c("low density","high density"))
rfnLong$method = factor(rfnLong$method, levels=c("SpiderLearner",
                                                 "simple mean",
                                                 "huge-ebic-0",
                                                 "huge-ebic-0.5",
                                                 "glasso-ric",
                                                 "hub glasso",
                                                 #"mle",
                                                 "stars-0.05",
                                                 "stars-0.1",
                                                 "qgraph-ebic-0",
                                                 "qgraph-ebic-0.5"))

rfnLong$topology = factor(rfnLong$topology, levels=c("erdos-renyi",
                                                     "small world",
                                                     "scale free",
                                                     "hub and spoke"))


rfnLong$`relative frobenius norm after`= c(sapply(sapply(allResults,function(x){return(x[2])}),as.vector))
rfnLong$`mrv`= c(sapply(sapply(allResults,function(x){return(x[3])}),as.vector))
rfnLong$`llTrain`= c(sapply(sapply(allResults,function(x){return(x[4])}),as.vector))
rfnLong$`llTest`= c(sapply(sapply(allResults,function(x){return(x[5])}),as.vector))

p<-ggplot(rfnLong, aes(x=method,y=`relative frobenius norm after`,color=method)) + 
  geom_boxplot() +
  #geom_jitter(aes(color = method), size=0.5, alpha=0.3) +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -90, hjust = 0,size = 18),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18)) + 
  facet_grid(topology ~ density, scales = "free_y") +
  labs(x="Estimation Method", 
       y="Relative Frobenius Norm",
       title="Simulation D: n=60,p=100,q=5050",
       subtitle="Relative Frobenius Norm")

ggsave("../../Figures/Figure4/simD_rfnAfter.jpeg", plot=p,width=10,height=10,units="in")

p<-ggplot(rfnLong, aes(x=method,y=`mrv`,color=method)) + 
  geom_boxplot() +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -90, hjust = 0,size = 18),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18)) + 
  facet_grid(topology ~ density, scales = "free_y") +
  labs(x="Estimation Method", 
       y="Matrix RV Coefficient",
       title="Simulation D: n=60,p=100,q=5050",
       subtitle="Matrix RV Coefficient")

ggsave("../../Figures/Figure4/simD_mrv.jpeg", plot=p,width=10,height=10,units="in",dpi=150)

p<-ggplot(rfnLong, aes(x=method,y=`llTrain`,color=method)) + 
  geom_boxplot() +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -90, hjust = 0,size = 18),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18)) + 
  facet_grid(topology ~ density, scales = "free_y") +
  labs(x="Estimation Method", 
       y="In-sample Log Likelihood",
       title="Simulation D: n=60,p=100,q=5050",
       subtitle="In-sample Log Likelihood")

ggsave("../../Figures/Figure4/simD_LLTrain.jpeg", plot=p,width=10,height=10,units="in",dpi=150)

p<-ggplot(rfnLong, aes(x=method,y=`llTest`,color=method)) + 
  geom_boxplot() +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -90, hjust = 0,size = 18),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18)) + 
  facet_grid(topology ~ density, scales = "free_y") +
  labs(x="Estimation Method", 
       y="Out-of-sample Log Likelihood",
       title="Simulation D: n=60,p=100,q=5050",
       subtitle="Out-of-sample Log Likelihood")

ggsave("../../Figures/Figure4/simD_LLTest.jpeg", plot=p,width=10,height=10,units="in")

## Boxplot density of each method for ebic. Is it giving us empty networks? usually.

allDensity = list()

for(i in 1:length(allResults))
{
  allDensity[[i]] = t(data.frame(sapply(allResults[[i]][[1]],function(x){return((sum(x$fullModels[[7]]!=0)-nPred)/(nPred^2))})))
}  

boxplot(allDensity)

## Calculate sens and spec

binaryDistance = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
tpArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
fpArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
tnArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
fnArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))

sampleSize = 60
P=100
#zCrit = qnorm(1-0.05/(P^2))
#rCrit = exp(zCrit/(60-3))/(exp(zCrit/(60-3))+1) # reasonable??

for(l in 1:8) # topologies
{
  for(m in 1:(nMod+2)) # methods
  {
    for(k in 1:nSim)
    {
      thres = 0 #rCrit - not fisher threshold for Sim D (see Web Appendix)
      if(m == 9)
        estMat = allResults[[l]]$ensModels[[k]]$optTheta
      if(m == 10)
        estMat = allResults[[l]]$ensModels[[k]]$simpleMeanNetwork
      if(m < 9)
        estMat = allResults[[l]]$ensModels[[k]]$fullModels[[m]]
      diff = binaryMatrixDiff(estMat, as.matrix(eightNetworks[[l]]),thres)
      tp = truePositiveCount(estMat, as.matrix(eightNetworks[[l]]),thres)
      fp = falsePositiveCount(estMat, as.matrix(eightNetworks[[l]]),thres)
      tn = trueNegativeCount(estMat, as.matrix(eightNetworks[[l]]),thres)
      fn = falseNegativeCount(estMat, as.matrix(eightNetworks[[l]]),thres)
      binaryDistance[l,m,k]=diff
      tpArray[l,m,k]=tp
      fpArray[l,m,k]=fp
      tnArray[l,m,k]=tn
      fnArray[l,m,k]=fn
    }
  }
}

sensArray = tpArray/(tpArray + fnArray)
specArray = tnArray/(fpArray + tnArray)
sensMeans = apply(sensArray,c(1,2),mean)
specMeans = apply(specArray,c(1,2),mean)

row.names(sensMeans) <- row.names(specMeans) <- c(
  "erLow",
  "erHigh",
  "wsLow",
  "wsHigh",
  "sfLow",
  "sfHigh",
  "hsLow",
  "hsHigh")

colnames(sensMeans) = methods
colnames(specMeans) = methods


write.table(round(sensMeans[,c(9:10,1:8)],2),"../../Tables//sensMeans_simD.asv",sep="&",row.names=T)
write.table(round(specMeans[,c(9:10,1:8)],4),"../../Tables//specMeans_simD.asv",sep="&",row.names=T)


### Binned Bias and MSE

biasMat = array(rep(NA,nTopology*(nMod+2)*nSim*nPred*nPred),dim=c(nTopology,(nMod+2),nSim,nPred,nPred))
mseMat = array(rep(NA,nTopology*(nMod+2)*nSim*nPred*nPred),dim=c(nTopology,(nMod+2),nSim,nPred,nPred))

for(l in 1:nTopology) # topologies
{
  for(m in 1:(nMod+2)) # methods
  {
    for(k in 1:nSim)
    {
      if(m==9) 
      {
        biasMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$optTheta - as.matrix(eightNetworks[[l]]))
        mseMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$optTheta - as.matrix(eightNetworks[[l]]))^2
      }
      if(m==10)
      {
        biasMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$simpleMeanNetwork - as.matrix(eightNetworks[[l]]))
        mseMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$simpleMeanNetwork - as.matrix(eightNetworks[[l]]))^2
      }
      if(!(m==9|m==10))
      {
        biasMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$fullModels[[m]] - as.matrix(eightNetworks[[l]]))
        mseMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$fullModels[[m]] - as.matrix(eightNetworks[[l]]))^2
      }
    }
  }
}

meanBias = function(l,m)
{
  # mean across the third index, which is simulation
  return(apply(biasMat[l,m,,,],c(2:3),mean))
}

splitBias = function(l,m)
{
  categories = splitToFive(eightNetworks[[l]])
  meanBiasMat = meanBias(l,m)
  category0Bias = meanBiasMat[categories==0]
  category1Bias = meanBiasMat[categories==1]
  category2Bias = meanBiasMat[categories==2]
  category3Bias = meanBiasMat[categories==3]
  category4Bias = meanBiasMat[categories==4]
  return(list(category0Bias,category1Bias,category2Bias,category3Bias,category4Bias))
}

meanBiasArray = array(rep(NA,nTopology*(nMod+2)*5),dim=c(nTopology,(nMod+2),5))
sdBiasArray = array(rep(NA,nTopology*(nMod+2)*5),dim=c(nTopology,(nMod+2),5))

topologies = c("erLow","erHigh","wsLow","wsHigh","sfLow","sfHigh","hsLow","hsHigh")
par(mfrow=c(8,10),mar=c(0,0,0,0))
for(l in 1:nTopology)
{
  for(m in 1:(nMod+2))
  {
    meanBiasArray[l,m,]=sapply(splitBias(l,m),mean)
    sdBiasArray[l,m,]=sapply(splitBias(l,m),sd)
  }
}

cats = c("zeroes","small","med","large","diag")
biasDF=as.data.frame.table(meanBiasArray, responseName = "meanBias",dnn=c(topologies,methods,cats))
library(plyr)
biasDF$Var1 = mapvalues(biasDF$Var1, from = c("A","B","C","D","E","F","G","H"), to = topologies)
biasDF$Var2 = mapvalues(biasDF$Var2, from = c("A","B","C","D","E","F","G","H","I","J"), to = methods)
biasDF$Var3 = mapvalues(biasDF$Var3, from = c("A","B","C","D","E"), to = cats)
sdDF=as.data.frame.table(sdBiasArray, responseName = "sdBias")
biasDF$sdBias = sdDF$sdBias
names(biasDF)=c("topology","method","category","meanBias","sdBias")

biasDF$method = factor(biasDF$method, levels=c("SpiderLearner",
                                               "simple mean",
                                               "glasso - ebic - 0",
                                               "glasso - ebic - 0.5",
                                               "glasso - ric",
                                               "hglasso",
                                               "glasso - stars - 0.05",
                                               "glasso - stars - 0.1",
                                               "qgraph - ebic - 0",
                                               "qgraph - ebic - 0.5"))
p=ggplot(biasDF, aes(x=method, y=meanBias, group=category,color=method)) + 
  facet_grid(category ~ topology,scales="free_y") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),legend.position="bottom") +
  geom_pointrange(size=0.5,aes(ymin=meanBias-sdBias, ymax=meanBias + sdBias)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  ggtitle("Simulation D: Element-wise Bias")

ggsave("../../Figures/Figure4/simD_bias.jpeg", plot=p,width=15,height=8,units="in",dpi=150)


### All the same for MSE

meanRMSE = function(l,m)
{
  # mean across the third index, which is simulation
  return(apply(mseMat[l,m,,,],c(2:3),function(x){return(sqrt(mean(x)))}))
}

splitRMSE = function(l,m)
{
  categories = splitToFive(eightNetworks[[l]])
  meanRMSEMat = meanRMSE(l,m)
  category0RMSE = meanRMSEMat[categories==0]
  category1RMSE = meanRMSEMat[categories==1]
  category2RMSE = meanRMSEMat[categories==2]
  category3RMSE = meanRMSEMat[categories==3]
  category4RMSE = meanRMSEMat[categories==4]
  return(list(category0RMSE,category1RMSE,category2RMSE,category3RMSE,category4RMSE))
}

meanRMSEArray = array(rep(NA,nTopology*(nMod+2)*5),dim=c(nTopology,(nMod+2),5))
sdRMSEArray = array(rep(NA,nTopology*(nMod+2)*5),dim=c(nTopology,(nMod+2),5))

topologies = c("erLow","erHigh","wsLow","wsHigh","sfLow","sfHigh","hsLow","hsHigh")
par(mfrow=c(8,10),mar=c(0,0,0,0))
for(l in 1:nTopology)
{
  for(m in 1:(nMod+2))
  {
    meanRMSEArray[l,m,]=sapply(splitRMSE(l,m),mean)
    sdRMSEArray[l,m,]=sapply(splitRMSE(l,m),sd)
  }
}

cats = c("zeroes","small","med","large","diag")
RMSEDF=as.data.frame.table(meanRMSEArray, responseName = "meanRMSE",dnn=c(topologies,methods,cats))
library(plyr)
RMSEDF$Var1 = mapvalues(RMSEDF$Var1, from = c("A","B","C","D","E","F","G","H"), to = topologies)
RMSEDF$Var2 = mapvalues(RMSEDF$Var2, from = c("A","B","C","D","E","F","G","H","I","J"), to = methods)
RMSEDF$Var3 = mapvalues(RMSEDF$Var3, from = c("A","B","C","D","E"), to = cats)
sdDF=as.data.frame.table(sdRMSEArray, responseName = "sdRMSE")
RMSEDF$sdRMSE = sdDF$sdRMSE
names(RMSEDF)=c("topology","method","category","meanRMSE","sdRMSE")

RMSEDF$method = factor(RMSEDF$method, levels=c("SpiderLearner",
                                               "simple mean",
                                               "glasso - ebic - 0",
                                               "glasso - ebic - 0.5",
                                               "glasso - ric",
                                               "hglasso",
                                               "glasso - stars - 0.05",
                                               "glasso - stars - 0.1",
                                               "qgraph - ebic - 0",
                                               "qgraph - ebic - 0.5"))
p=ggplot(RMSEDF, aes(x=method, y=meanRMSE, group=category,color=method)) + 
  facet_grid(category ~ topology,scales="free_y") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),legend.position="bottom") +
  geom_pointrange(size=0.5,aes(ymin=meanRMSE-sdRMSE, ymax=meanRMSE + sdRMSE)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  ylab("Average MSE over category elements") +
  ggtitle("Simulation D: Element-wise MSE")

ggsave("../../Figures/Figure4/simD_rmse.jpeg", plot=p,width=15,height=8,units="in",dpi=150)


