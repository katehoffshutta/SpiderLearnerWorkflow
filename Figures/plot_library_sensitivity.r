library(MASS)
## Process results

rfnSmall = list()
rfnMed = list()
rfnLarge = list()

ooslSmall = list()
ooslMed = list()
ooslLarge = list()

load("Results/Pilot/librarySimRes.rda")
load("Results/eightNetworksAC.rda")
eightNetworks = ac
set.seed(1492)

newData = mvrnorm(1000, mu = rep(0,nrow(precMat)), solve(precMat))


for(i in 1:nSim)
{
  rfnSmall[i] = relativeFrobNormAfter(librarySimRes[[i]][[1]]$optTheta,eightNetworks[[6]])
  rfnMed[i] = relativeFrobNormAfter(librarySimRes[[i]][[2]]$optTheta,eightNetworks[[6]])
  rfnLarge[i] = relativeFrobNormAfter(librarySimRes[[i]][[3]]$optTheta,eightNetworks[[6]])

  ooslSmall[i] = logLik(librarySimRes[[i]][[1]]$optTheta,newData)
  ooslMed[i] = logLik(librarySimRes[[i]][[2]]$optTheta,newData)
  ooslLarge[i] = logLik(librarySimRes[[i]][[3]]$optTheta,newData)
  
}

pdf("libraryBoxplots.pdf")
par(mfrow=c(1,2))
boxplot(unlist(rfnSmall),unlist(rfnMed),unlist(rfnLarge),
        names=c("Small","Medium","Large"), col=c(2,3,4),
        xlab = "Library Size", ylab = "Relative Frobenius Norm",
        main = "Relative Frobenius Norm \n by Library Size")
boxplot(unlist(ooslSmall),unlist(ooslMed),unlist(ooslLarge),
        names=c("Small","Medium","Large"), col=c(2,3,4),
        xlab = "Library Size", ylab = "Out-of-sample Likelihood",
        main = "Out-of-sample Likelihood \n by Library Size")
dev.off()

# Look at the profiles to assess trends
plot(0,0,ylim=c(0.09,0.115),xlim=c(1,3))
for(i in 1:100)
{
  lines(c(1,2,3),c(rfnSmall[[i]],rfnMed[[i]],rfnLarge[[i]]),col=rainbow(100)[i])
} 

plot(0,0,ylim=c(-25640,-25550),xlim=c(1,3))
for(i in 1:100)
{
  lines(c(1,2,3),c(ooslSmall[[i]],ooslMed[[i]],ooslLarge[[i]]),col=rainbow(100)[i])
} 
