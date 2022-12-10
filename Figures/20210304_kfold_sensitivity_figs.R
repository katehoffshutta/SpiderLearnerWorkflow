library(MASS)
load("Results/eightNetworksAC.rda")

precMat = ac[[2]]

load("bootstrap/simRes_er_50.rda")
simResA= simRes
load("bootstrap/simRes_er_50_s42.rda")
simResB = simRes

simRes = list()
for(i in 1:50)
  simRes[[i]] = simResA[[i]]

for(i in 51:100)
  simRes[[i]] = simResB[[i-50]]

# get standard deviation

R=2
k2 = array(rep(NA,R*50*50),dim=c(R,50,50))
k5 = array(rep(NA,R*50*50),dim=c(R,50,50))
k10 = array(rep(NA,R*50*50),dim=c(R,50,50))
k15 = array(rep(NA,R*50*50),dim=c(R,50,50))
k20 = array(rep(NA,R*50*50),dim=c(R,50,50))
k30 = array(rep(NA,R*50*50),dim=c(R,50,50))

R=2
for(r in 1:R)
{
  k2[r,,] = simRes[[r]][[1]]$optTheta
  k5[r,,] = simRes[[r]][[2]]$optTheta
  k10[r,,] = simRes[[r]][[3]]$optTheta
  k15[r,,] = simRes[[r]][[4]]$optTheta
  k20[r,,] = simRes[[r]][[5]]$optTheta
  k30[r,,] = simRes[[r]][[6]]$optTheta
}

sd2 = apply(k2,2:3,sd)
sd5 = apply(k5,2:3,sd)
sd10 = apply(k10,2:3,sd)
sd15 = apply(k15,2:3,sd)
sd20 = apply(k20,2:3,sd)
sd30 = apply(k30,2:3,sd)

k2Bias = array(rep(NA,R*50*50),dim=c(R,50,50))
k5Bias = array(rep(NA,R*50*50),dim=c(R,50,50))
k10Bias = array(rep(NA,R*50*50),dim=c(R,50,50))
k15Bias = array(rep(NA,R*50*50),dim=c(R,50,50))
k20Bias = array(rep(NA,R*50*50),dim=c(R,50,50))
k30Bias = array(rep(NA,R*50*50),dim=c(R,50,50))

for(r in 1:R)
{
  k2Bias[r,,] = k2[r,,]-as.matrix(precMat)
  k5Bias[r,,] = k5[r,,]-as.matrix(precMat)
  k10Bias[r,,] = k10[r,,]-as.matrix(precMat)
  k15Bias[r,,] = k15[r,,]-as.matrix(precMat)
  k20Bias[r,,] = k20[r,,]-as.matrix(precMat)
  k30Bias[r,,] = k30[r,,]-as.matrix(precMat)
}

mean2 = apply(k2Bias,2:3,mean)
mean5 = apply(k5Bias,2:3,mean)
mean10 = apply(k10Bias,2:3,mean)
mean15 = apply(k15Bias,2:3,mean)
mean20 = apply(k20Bias,2:3,mean)
mean30 = apply(k30Bias,2:3,mean)

## Plot pooled sd
plot(0,0,xlim=c(0,30),ylim=c(-7,-2),xlab="K",ylab="entrywise log SD", main="each profile = one matrix entry")
for(i in 1:50)
{
  for(j in 1:i)
  {
    profile=c(sd2[i,j],
              sd5[i,j],
              sd10[i,j],
              sd15[i,j],
              sd20[i,j],
              sd30[i,j])
    lines(c(2,5,10,15,20,30),log(profile),col="lightgrey",lwd=0.3)
  }
}

## Plot pooled mean bias
plot(0,0,xlim=c(0,30),ylim=c(-0.2,0.2),xlab="K",ylab="entrywise log SD", main="each profile = one matrix entry")
for(i in 1:50)
{
  for(j in 1:i)
  {
    profile=c(mean2[i,j],
              mean5[i,j],
              mean10[i,j],
              mean15[i,j],
              mean20[i,j],
              mean30[i,j])
    lines(c(2,5,10,15,20,30),profile,col="lightgrey",lwd=0.5)
  }
}

nonzeroEntries = as.matrix(precMat)[as.matrix(precMat-diag(50))!=0]
quantile(abs(nonzeroEntries))
lowWeightCutoff = quantile(abs(nonzeroEntries),0.1)
highWeightCutoff = quantile(abs(nonzeroEntries),0.9)

zeroEntries = data.frame(0,0)
smallEntries = data.frame(0,0)
mediumEntries = data.frame(0,0)
largeEntries = data.frame(0,0)

for(i in 2:50)
{
  for(j in 1:(i-1))
  {
    print(j)
    if(abs(precMat[i,j]) == 0) 
      zeroEntries = rbind(zeroEntries,c(i,j))
    if(0 < abs(precMat[i,j]) & abs(precMat[i,j]) <= lowWeightCutoff) 
      smallEntries = rbind(smallEntries,c(i,j))
    if(lowWeightCutoff < abs(precMat[i,j]) & abs(precMat[i,j]) <= highWeightCutoff) 
      mediumEntries = rbind(mediumEntries,c(i,j))
    if(highWeightCutoff < abs(precMat[i,j])) 
      largeEntries = rbind(largeEntries,c(i,j))
  }
}

zeroEntries = zeroEntries[-1,]
smallEntries = smallEntries[-1,]
mediumEntries = mediumEntries[-1,]
largeEntries = largeEntries[-1,]

jpeg("kProfilesBinned_RandomGraph.jpeg",width=8,height=8,units="in",res=400)
par(mfrow=c(2,2))
plot(0,0,xlim=c(0,30),ylim=c(-7,-2),xlab="K",ylab="entrywise log SD", main="zero entries \n each profile = one matrix entry")
zeroEntrySum = rep(0,6)
for(m in 1:nrow(zeroEntries))
{
    i = zeroEntries[m,1]
    j = zeroEntries[m,2]
    profile=c(sd2[i,j],
              sd5[i,j],
              sd10[i,j],
              sd15[i,j],
              sd20[i,j],
              sd30[i,j])
    lines(c(2,5,10,15,20,30),log(profile),col="lightgrey",lwd=0.2)
    zeroEntrySum = zeroEntrySum + log(profile)
}
lines(c(2,5,10,15,20,30), zeroEntrySum/nrow(zeroEntries),col="black",lwd=3)


plot(0,0,xlim=c(0,30),ylim=c(-7,-2),xlab="K",ylab="entrywise log SD", main="small entries \n each profile = one matrix entry")
smallEntrySum = rep(0,6)
for(m in 1:nrow(smallEntries))
{
  i = smallEntries[m,1]
  j = smallEntries[m,2]
  profile=c(sd2[i,j],
            sd5[i,j],
            sd10[i,j],
            sd15[i,j],
            sd20[i,j],
            sd30[i,j])
  lines(c(2,5,10,15,20,30),log(profile),col="lightgrey",lwd=0.2)
  smallEntrySum = smallEntrySum + log(profile)
}

lines(c(2,5,10,15,20,30), smallEntrySum/nrow(smallEntries),col="black",lwd=3)
      
plot(0,0,xlim=c(0,30),ylim=c(-7,-2),xlab="K",ylab="entrywise log SD", main="medium entries \n each profile = one matrix entry")

mediumEntrySum = rep(0,6)
for(m in 1:nrow(mediumEntries))
{
  i = mediumEntries[m,1]
  j = mediumEntries[m,2]
  profile=c(sd2[i,j],
            sd5[i,j],
            sd10[i,j],
            sd15[i,j],
            sd20[i,j],
            sd30[i,j])
  lines(c(2,5,10,15,20,30),log(profile),col="lightgrey",lwd=0.2)
  mediumEntrySum = mediumEntrySum + log(profile)
}

lines(c(2,5,10,15,20,30), mediumEntrySum/nrow(mediumEntries),col="black",lwd=3)


plot(0,0,xlim=c(0,30),ylim=c(-7,-2),xlab="K",ylab="entrywise log SD", main="large entries \n each profile = one matrix entry")

largeEntrySum = c(0,6)
for(m in 1:nrow(largeEntries))
{
  i = largeEntries[m,1]
  j = largeEntries[m,2]
  profile=c(sd2[i,j],
            sd5[i,j],
            sd10[i,j],
            sd15[i,j],
            sd20[i,j],
            sd30[i,j])
  lines(c(2,5,10,15,20,30),log(profile),col="lightgrey",lwd=0.2)
  largeEntrySum = largeEntrySum + log(profile)
  
}

lines(c(2,5,10,15,20,30), largeEntrySum/nrow(largeEntries),col="black",lwd=3)

dev.off()

jpeg("kProfilesBinned_RandomGraph_bias.jpeg",width=8,height=8,units="in",res=400)
par(mfrow=c(2,2))
plot(0,0,xlim=c(0,30),ylim=c(-0.01,0.01),xlab="K",ylab="entrywise bias", main="zero entries \n each profile = one matrix entry")

zeroEntrySum = rep(0,6)
for(m in 1:nrow(zeroEntries))
{
  i = zeroEntries[m,1]
  j = zeroEntries[m,2]
  profile=c(mean2[i,j],
            mean5[i,j],
            mean10[i,j],
            mean15[i,j],
            mean20[i,j],
            mean30[i,j])
  lines(c(2,5,10,15,20,30),profile,col="lightgrey",lwd=0.2)
  zeroEntrySum = zeroEntrySum + profile
}

lines(c(2,5,10,15,20,30),zeroEntrySum/nrow(zeroEntries),col="black",lwd=3)
legend(20,0.2,c("average bias"),col="black",lty=1,lwd=3)

plot(0,0,xlim=c(0,30),ylim=c(-0.05,0.05),xlab="K",ylab="entrywise mean", main="small entries \n each profile = one matrix entry")
smallEntrySum = rep(0,6)
for(m in 1:nrow(smallEntries))
{
  i = smallEntries[m,1]
  j = smallEntries[m,2]
  profile=c(mean2[i,j],
            mean5[i,j],
            mean10[i,j],
            mean15[i,j],
            mean20[i,j],
            mean30[i,j])
  lines(c(2,5,10,15,20,30),profile,col="lightgrey",lwd=0.2)
  smallEntrySum = smallEntrySum + profile
}

lines(c(2,5,10,15,20,30),smallEntrySum/nrow(smallEntries),col="black",lwd=3)
legend(20,0.2,c("average bias"),col="black",lty=1,lwd=3)

plot(0,0,xlim=c(0,30),ylim=c(-0.1,0.1),xlab="K",ylab="entrywise mean", main="medium entries \n each profile = one matrix entry")
mediumEntrySum = rep(0,6)
for(m in 1:nrow(mediumEntries))
{
  i = mediumEntries[m,1]
  j = mediumEntries[m,2]
  profile=c(mean2[i,j],
            mean5[i,j],
            mean10[i,j],
            mean15[i,j],
            mean20[i,j],
            mean30[i,j])
  lines(c(2,5,10,15,20,30),profile,col="lightgrey",lwd=0.2)
  mediumEntrySum = mediumEntrySum + profile
}

lines(c(2,5,10,15,20,30),mediumEntrySum/nrow(mediumEntries),col="black",lwd=3)
legend(20,0.2,c("average bias"),col="black",lty=1,lwd=3)

plot(0,0,xlim=c(0,30),ylim=c(-0.2,0.2),xlab="K",ylab="entrywise mean", main="large entries \n each profile = one matrix entry")
largeEntrySum = rep(0,6)
for(m in 1:nrow(largeEntries))
{
  i = largeEntries[m,1]
  j = largeEntries[m,2]
  profile=c(mean2[i,j],
            mean5[i,j],
            mean10[i,j],
            mean15[i,j],
            mean20[i,j],
            mean30[i,j])
  lines(c(2,5,10,15,20,30),profile,col="lightgrey",lwd=0.2)
  largeEntrySum = largeEntrySum + profile
}
lines(c(2,5,10,15,20,30),largeEntrySum/nrow(largeEntries),col="black",lwd=3)

dev.off()

## Question 2: Out of sample likelihood?

oosls2 = c()
oosls5 = c()
oosls10 = c()
oosls15 = c()
oosls20 = c()
oosls30 = c()

set.seed(900)
newData = mvrnorm(150, mu = rep(0,nrow(precMat)), solve(precMat))
for(i in 1:100)
{
  thisTheta2 = simRes[[i]][[1]]$optTheta
  thisTheta5 = simRes[[i]][[2]]$optTheta
  thisTheta10 = simRes[[i]][[3]]$optTheta
  thisTheta15 = simRes[[i]][[4]]$optTheta
  thisTheta20 = simRes[[i]][[5]]$optTheta
  thisTheta30 = simRes[[i]][[6]]$optTheta
  
  oosls2 = c(oosls2,logLik(thisTheta2,newData))
  oosls5 = c(oosls5,logLik(thisTheta5,newData))
  oosls10 = c(oosls10,logLik(thisTheta10,newData))
  oosls15 = c(oosls15,logLik(thisTheta15,newData))
  oosls20 = c(oosls20,logLik(thisTheta20,newData))
  oosls30 = c(oosls30,logLik(thisTheta30,newData))
}

pdf("kOOSL.pdf")
boxplot(oosls2,oosls5,oosls10,oosls15,oosls20,oosls30,
        names=c(2,5,10,15,20,30),xlab="Number of Folds (K)",
        ylab="Out-of-sample Log Likelihood",
        main="Effect of K on Variability of OOSL")
dev.off()

### How does the choice of k affect the variability of the weights in the model?

weights2 = matrix(rep(NA,9*100),ncol=9)
weights5 = matrix(rep(NA,9*100),ncol=9)
weights10 = matrix(rep(NA,9*100),ncol=9)
weights15 = matrix(rep(NA,9*100),ncol=9)
weights20 = matrix(rep(NA,9*100),ncol=9)
weights30 = matrix(rep(NA,9*100),ncol=9)

for(i in 1:100)
{
  weights2[i,] = simRes[[i]][[1]]$weights
  weights5[i,] = simRes[[i]][[2]]$weights
  weights10[i,] = simRes[[i]][[3]]$weights
  weights15[i,] = simRes[[i]][[4]]$weights
  weights20[i,] = simRes[[i]][[5]]$weights
  weights30[i,] = simRes[[i]][[6]]$weights
}

pdf("kWeights.pdf")
par(mfrow=c(3,2))
boxplot(weights2,ylim=c(0,0.7),main="2-fold")
boxplot(weights5,ylim=c(0,0.7),main="5-fold")
boxplot(weights10,ylim=c(0,0.7),main="10-fold")
boxplot(weights15,ylim=c(0,0.7),main="15-fold")
boxplot(weights20,ylim=c(0,0.7),main="20-fold")
boxplot(weights30,ylim=c(0,0.7),main="30-fold")
dev.off()

