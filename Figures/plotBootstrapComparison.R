# plot the results of the bootstrap comparison

res = read.csv("../Results/bootstrapComparisonResults.csv")

hist(res$diffOOSL/abs(res$bootstrapOOSL))
# mean percent increase in OOSL in ensemble vs. bootstrap
mean(res$diffOOSL/abs(res$bootstrapOOSL))
sd(res$diffOOSL/abs(res$bootstrapOOSL))

jpeg("../Figures/chen_bootstrap_example.jpeg",width=8,height=6,units="in",res=200)
boxplot(res$bootstrapOOSL, res$ensOOSL, 
        names=c("Bootstrap Ensemble","SpiderLearner Ensemble"),
        xlab = "Model",
        ylab = "log likelihood on test data",
        main = "Out-of-sample performance of \n SpiderLearner vs. bootstrapping an ensemble")
dev.off()