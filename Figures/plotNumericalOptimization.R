# plot the numerical optimization simulation results

res = read.csv("../Results/numericalOptimizationStudyWeights.csv")

names(res) = c("glasso - ebic - 0", 
                         "glasso - ebic - 0.5", 
                         "glasso - ric",
                         "hglasso",
                         "mle",
                         "glasso - stars - 0.05",
                         "glasso - stars - 0.1",
                         "qgraph - ebic - 0",
                         "qgraph - ebic - 0.5")

pdf("revision_numerical_optimization_seed.pdf",width=10,height=5)
par(cex.axis=0.5)
boxplot(res,xlab="candidate",
        ylab="weight in ensemble",
        ylim=c(0,1),
        main="Variability in Ensemble Coefficients \n n=100 iterations, folds fixed")
dev.off()

round(apply(res,2,sd),5)

# look at the actual variability of the models

mods = load("../Results/numericalOptimizationModels.rda")
library(abind)

modarray = abind(ensModels,along=3)
modmean = apply(modarray,1:2,mean)
modsd = apply(modarray,1:2,sd)

cv = abs(modsd/modmean)
pdf("revision_numerical_optimization_cv.pdf",width=10,height=4)
plot(modmean[lower.tri(modmean,diag=F)],cv[lower.tri(cv,diag=F)],
     pch=20,cex=0.62, xlab="mean edge weight",ylab="cv of edge weight",
     main="Coefficient of Variability of SpiderLearner \n Edge Weights (N=100 Simulations)")
dev.off()
mean(cv)
