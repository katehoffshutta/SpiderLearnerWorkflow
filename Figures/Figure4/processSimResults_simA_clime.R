## Process simulation A results to generate table and figure source for paper

library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
setwd("~/research/repos/SpiderLearnerWorkflow/Figures/Figure4/")

source("../../Simulations/errorMetrics.R")

load("../../Results/SimA/clime/20221011_erHighPrec_RealData_n_10000_p_50_simStudy.rda")
erHighPrecResults = theseResults
load("../../Results/SimA/clime/20221011_erLowPrec_RealData_n_10000_p_50_simStudy.rda")
erLowPrecResults = theseResults
load("../../Results/SimA/clime/20221011_wsHighPrec_RealData_n_10000_p_50_simStudy.rda")
wsHighPrecResults = theseResults
load("../../Results/SimA/clime/20221011_wsLowPrec_RealData_n_10000_p_50_simStudy.rda")
wsLowPrecResults = theseResults
load("../../Results/SimA/clime/20221011_sfHighPrec_RealData_n_10000_p_50_simStudy.rda")
sfHighPrecResults = theseResults
load("../../Results/SimA/clime/20221011_sfLowPrec_RealData_n_10000_p_50_simStudy.rda")
sfLowPrecResults = theseResults
load("../../Results/SimA/clime/20221011_hsHighPrec_RealData_n_10000_p_50_simStudy.rda")
hsHighPrecResults = theseResults
load("../../Results/SimA/clime/20221011_hsLowPrec_RealData_n_10000_p_50_simStudy.rda")
hsLowPrecResults = theseResults

allResults = list(erLowPrecResults,
                  erHighPrecResults,
                  wsLowPrecResults,
                  wsHighPrecResults,
                  sfLowPrecResults,
                  sfHighPrecResults,
                  hsLowPrecResults,
                  hsHighPrecResults)

load("../../Results/eightNetworksAC.rda")
eightNetworks = ac

methods = c("glasso - ebic - 0", 
            "glasso - ebic - 0.5", 
                      "glasso - ric",
                      "hglasso",
                      "mle",
                      "glasso - stars - 0.05",
                      "glasso - stars - 0.1",
                      "qgraph - ebic - 0",
                      "qgraph - ebic - 0.5",
                      "clime",
                      "SpiderLearner",
                      "simple mean")
nPred = 50
nTopology = 8
nMod = length(methods)-2 # number of models, not counting SpiderLearner and simple mean
  
allWeights = list()

for(i in 1:length(allResults))
{
  allWeights[[i]] = t(data.frame(sapply(allResults[[i]][[1]],function(x)x$weights)))
}         

# If desired, inspect weight distributions w/boxplots

### par(mfrow=c(4,2))

### boxplot(allWeights[[1]],main="Random Graph Low Density SpiderLearner Weights",names = methods[1:nMod],
###         las=2)
### boxplot(allWeights[[2]],main="Random Graph High Density SpiderLearner Weights",names = methods[1:nMod],las=2)

### boxplot(allWeights[[3]],main="Small World Graph Low Density SpiderLearner Weights",names = methods[1:nMod],
###         las=2)
### boxplot(allWeights[[4]],main="Small World Graph High Density SpiderLearner Weights",names = methods[1:nMod],las=2)

### boxplot(allWeights[[5]],main="Scale Free Graph Low Density SpiderLearner Weights",names = methods[1:nMod],
###         las=2)
### boxplot(allWeights[[6]],main="Scale Free Graph High Density SpiderLearner Weights",names = methods[1:nMod],las=2)

### boxplot(allWeights[[7]],main="Hub-and-Spoke Graph Low Density SpiderLearner Weights",names = methods[1:nMod],
###       las=2)
### boxplot(allWeights[[8]],main="Hub-and-Spoke Graph High Density SpiderLearner Weights",names = methods[1:nMod],las=2)

weightsTable = data.frame("erLow"=apply(allWeights[[1]],2,mean),
                          "erHigh"=apply(allWeights[[2]],2,mean),
                          "wsLow"=apply(allWeights[[3]],2,mean),
                          "wsHigh"=apply(allWeights[[4]],2,mean),
                          "sfLow"=apply(allWeights[[5]],2,mean),
                          "sfHigh"=apply(allWeights[[6]],2,mean),
                          "hsLow"=apply(allWeights[[7]],2,mean),
                          "hsHigh"=apply(allWeights[[8]],2,mean))

row.names(weightsTable) = methods[1:nMod]
write.table(round(t(weightsTable),2),file="../../Tables/weights_simA_clime.asv",sep="&",row.names=T,quote=F)

nSim=100
nMod=10
rfnLong = data.frame("topology"=c(rep("erdos-renyi",(nMod+2)*2*nSim), 
                                  rep("small world",(nMod+2)*2*nSim),
                                  rep("scale free",(nMod+2)*2*nSim),
                                  rep("hub and spoke",(nMod+2)*2*nSim)))
rfnLong$density = rep(c(rep("low density",(nMod+2)*nSim),rep("high density",(nMod+2)*nSim)),4)
rfnLong$method = rep(c(rep("huge-ebic-0",nSim),
                   rep("huge-ebic-0.5",nSim),
                   rep("glasso-ric",nSim),
                   rep("hub glasso",nSim),
                   rep("mle",nSim),
                   rep("stars-0.05",nSim),
                   rep("stars-0.1",nSim),
                   rep("qgraph-ebic-0",nSim),
                   rep("qgraph-ebic-0.5",nSim),
                   rep("clime",nSim),
                   rep("SpiderLearner",nSim),
                   rep("simple mean",nSim)),8)

rfnLong$density = factor(rfnLong$density, levels=c("low density","high density"))
rfnLong$method = factor(rfnLong$method, levels=c("SpiderLearner",
                                                 "simple mean",
                                                 "huge-ebic-0",
                                                 "huge-ebic-0.5",
                                                 "glasso-ric",
                                                 "hub glasso",
                                                 "mle",
                                                 "stars-0.05",
                                                 "stars-0.1",
                                                 "qgraph-ebic-0",
                                                 "qgraph-ebic-0.5",
                                                 "clime"))

rfnLong$topology = factor(rfnLong$topology, levels=c("erdos-renyi",
                                                     "small world",
                                                     "scale free",
                                                     "hub and spoke"))

rfnLong$`relative frobenius norm after`= c(sapply(sapply(allResults,function(x){return(x[2])}),as.vector))
rfnLong$`mrv`= c(sapply(sapply(allResults,function(x){return(x[3])}),as.vector))
rfnLong$`llTrain`= c(sapply(sapply(allResults,function(x){return(x[4])}),as.vector))
rfnLong$`llTest`= c(sapply(sapply(allResults,function(x){return(x[5])}),as.vector))

## Added 20210921: Extra plot of just the 4 low density for dissertation defense

defense_names <- c(
  `erdos-renyi` = "Random",
  `small world` = "Small World",
  `scale free` = "Scale-Free",
  `hub and spoke`="Hub-and-Spoke"
)

p<-ggplot(rfnLong[rfnLong$density == "low density",], aes(x=method,y=`relative frobenius norm after`,fill=method)) + 
  geom_boxplot() +
  theme(legend.position="none",
        axis.text.x=element_text(angle = 45, hjust = 1,size = 18),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=24),
        axis.text.y=element_text(size=24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18)) + 
  facet_wrap(~topology,ncol=4,scales="free",labeller = as_labeller(defense_names)) +
  labs(x="Estimation Method", 
       y="Error (Relative Frobenius Norm)")
       #title="Simulation A: n=10000,p=50,q=1275",
       #subtitle="Relative Frobenius Norm")

ggsave("simResults_defense_clime.jpeg",width=20,height=8,units="in")
## Paper figures

p<-ggplot(rfnLong, aes(x=method,y=`relative frobenius norm after`,color=method)) + 
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
       y="Relative Frobenius Norm",
       title="Simulation A: n=10000,p=50,q=1275",
       subtitle="Relative Frobenius Norm")

ggsave("../../Figures/Figure4/simA_rfnAfter_clime.jpeg", plot=p,width=10,height=10,units="in")

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
       title="Simulation A: n=10000,p=50,q=1275",
       subtitle="Matrix RV Coefficient")

ggsave("../../Figures/Figure4/simA_mrv_clime.jpeg", plot=p,width=10,height=10,units="in",dpi=150)

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
       title="Simulation A: n=10000,p=50,q=1275",
       subtitle="In-sample Log Likelihood")

ggsave("../../Figures/Figure4/simA_LLTrain_clime.jpeg", plot=p,width=10,height=10,units="in",dpi=150)

p<-ggplot(rfnLong, aes(x=method,y=`llTest`,color=method)) + 
  geom_boxplot() +
  facet_grid(topology ~ density, scales = "free_y") +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -90, hjust = 0,size = 18),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=18),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18)) + 
  labs(x="Estimation Method", 
       y="Out-of-sample Log Likelihood",
       title="Simulation A: n=10000,p=50,q=1275",
       subtitle="Out-of-sample Log Likelihood")

ggsave("../../Figures/Figure4/simA_LLTest_clime.jpeg", plot=p,width=10,height=10,units="in")

## Calculate sens and spec

binaryDistance = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
tpArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
fpArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
tnArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))
fnArray = array(rep(NA, 8*(nMod+2)*nSim),dim=c(8,(nMod+2),nSim))

sampleSize = 10000
P=50
zCrit = qnorm(1-0.05/(P^2))
rCrit = (exp(zCrit/sqrt(sampleSize-3-(P-2))*2)-1)/(exp(zCrit/sqrt(sampleSize-3-(P-2))*2)+1)

for(l in 1:8) # topologies
{
  for(m in 1:(nMod+2)) # methods
  {
    for(k in 1:nSim)
    {
      thres = rCrit 
      if(m == 10)
        estMat = allResults[[l]]$ensModels[[k]]$optTheta
      if(m == 11)
        estMat = allResults[[l]]$ensModels[[k]]$simpleMeanNetwork
      if(m < 10)
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


write.table(round(sensMeans[,c(10:11,1:9)],2),"../../Tables/sensMeans_simA.asv",sep="&",row.names=T)
write.table(round(specMeans[,c(10:11,1:9)],4),"../../Tables/specMeans_simA.asv",sep="&",row.names=T)

### Binned Bias and MSE


biasMat = array(rep(NA,nTopology*(nMod+2)*nSim*nPred*nPred),dim=c(nTopology,(nMod+2),nSim,nPred,nPred))
mseMat = array(rep(NA,nTopology*(nMod+2)*nSim*nPred*nPred),dim=c(nTopology,(nMod+2),nSim,nPred,nPred))

for(l in 1:nTopology) # topologies
{
  for(m in 1:(nMod+2)) # methods
  {
    for(k in 1:nSim)
    {
      if(m==10) 
      {
        biasMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$optTheta - as.matrix(eightNetworks[[l]]))
        mseMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$optTheta - as.matrix(eightNetworks[[l]]))^2
      }
      if(m==11)
      {
        biasMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$simpleMeanNetwork - as.matrix(eightNetworks[[l]]))
        mseMat[l,m,k,,] = (allResults[[l]]$ensModels[[k]]$simpleMeanNetwork - as.matrix(eightNetworks[[l]]))^2
      }
      if(!(m==10|m==11))
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
par(mfrow=c(8,11),mar=c(0,0,0,0))
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
biasDF$Var2 = mapvalues(biasDF$Var2, from = c("A","B","C","D","E","F","G","H","I","J","K"), to = methods)
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
                                                 "mle",
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
  ggtitle("Simulation A: Element-wise Bias")
  
ggsave("../../Figures/Figure4/simA_bias_clime.jpeg", plot=p,width=15,height=8,units="in",dpi=150)


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
par(mfrow=c(8,11),mar=c(0,0,0,0))
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
RMSEDF$Var2 = mapvalues(RMSEDF$Var2, from = c("A","B","C","D","E","F","G","H","I","J","K"), to = methods)
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
                                               "mle",
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
  ggtitle("Simulation A: Element-wise MSE")

ggsave("../../Figures/Figure4/simA_rmse_clime.jpeg", plot=p,width=15,height=8,units="in",dpi=150)
