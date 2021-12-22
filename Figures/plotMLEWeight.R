load("Results/mleWeights_pkg.rda")
sparsityVec = c(0.05,0.10,0.25,0.5,0.75,1)
plot(0,0,xlim=c(0,1),ylim=c(0,1))

jpeg("Figures/mleWeight.jpeg",width=8,height=8,units="in",res=300)

plot(0,0,xlim=c(0,1),ylim=c(0,0.7),cex=0,xlab="Density = 1-Sparsity",ylab="Weight", main="Weight of MLE in Ensemble")

nIt = 30

for(n in 1:nIt)
{
  points(sparsityVec,ensModelWeights[n,,5],cex=0.2)
}

allMeans = apply(ensModelWeights,2:3,mean)
for(i in 5:5)
{
  points(sparsityVec, allMeans[,i],cex=2)
}

legend(0.55,0.1,c("individual simulation","mean of N=30 simulations"),pch=c(1,1),pt.cex=c(0.2,2))
dev.off()