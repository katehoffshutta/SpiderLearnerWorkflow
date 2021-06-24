unbd = read.csv("../Simulations/ensModelWeightsUnbd.csv")
bd = read.csv("../Simulations/ensModelWeightsBd.csv")

methods = c("glasso-ebic-0","glasso-ebic-0.5","glasso-ric","hglasso","mle")

pdf("boundedLoss2.pdf",width=20,height=8)
par(mfrow=c(1,2))
boxplot(unbd,names=methods,ylab="weights",main="Unbounded Loss Function Weights")
boxplot(bd,names=methods,ylab="weights",main="Bounded Loss Function Weights")
dev.off()
