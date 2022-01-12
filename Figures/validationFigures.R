# Make Figure 4a

df = read.csv("internalCVResults.csv")[-1]

df$fold = 1:10
df %>% as.data.frame  %>%
  pivot_longer(cols = 1:10)  %>%
  ggplot(aes(x=name,y=value)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1,size=12)) +
  xlab("Method") +
  ylab("Out-of-sample Log Likelihood") + 
  ggtitle("10-fold Internal Cross-Validation on Yoshihara Dataset") +
  geom_jitter(width = 0.1,alpha=0.5,size=1.75)

ggsave("ooslValidation.pdf")

# Make Figure 4b and 4c

scaledVal = read.csv("scaledVal.csv")[-1]

scaledVal %>% na.omit() %>% pivot_longer(cols = names(scaledVal)[1:10])  %>%
  ggplot(aes(x=name,y=100*value)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1,size=12)) +
  xlab("Method") +
  ylab("Percent Difference in Log-Likelihood from Best Model") +
  ggtitle("External Validation on 11 Independent Datasets") +
  geom_jitter(width = 0.1,alpha=0.5,size=1.75)

ggsave("ovarianExternalValidation.pdf")

scaledVal %>%
  na.omit() %>% pivot_longer(cols = names(scaledVal)[1:10])  %>%
  filter(name %in% c("glasso.ebic.0" ,"hglasso","qgraph.ebic.0","SpiderLearner")) %>%
  ggplot(aes(x=name,y=100*value)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1,size=12)) +
  xlab("Method") +
  ylab("Percent Difference in Log-Likelihood from Best Model") +
  ggtitle("External Validation on 11 Independent Datasets") +
  geom_jitter(width = 0.1,alpha=0.5,size=1.75)

ggsave("ovarianExternalValidationZoom.pdf")

# Make Figure 4d:community membership

load("ocSpiderLearnerResults.rda")

adjMat = -cov2cor(slResults$optTheta)
colnames(adjMat)=colnames(lateStageSmall)
adjGraph = graph_from_adjacency_matrix(adjMat,weighted=T, diag=F, mode="undirected")

absGraph = adjGraph
E(absGraph)$weight = abs(E(absGraph)$weight)
optThetaComm = cluster_fast_greedy(absGraph)

clusterList = list()

clusterList$ensemble = optThetaComm

ensembleBosses = c()
for(i in unique(factor(optThetaComm$membership)))
{
  thisComm = optThetaComm$names[which(optThetaComm$membership==i)]
  thisSubgraph = induced_subgraph(adjGraph,V(adjGraph)[which(colnames(lateStageSmall) %in% thisComm)])
  commHubs = hub_score(thisSubgraph)$vector
  localBoss = thisComm[which.max(commHubs)]
  ensembleBosses = c(ensembleBosses, localBoss)
}

set.seed(5555)

jpeg("../Figures/ensembleCommunitiesWithHubs_final.jpeg",width=9,height=9,units="in",res=300)
plot(absGraph,layout=layout_with_fr,
     vertex.size=ifelse(colnames(lateStageSmall)%in%ensembleBosses,22,6),
     vertex.label=ifelse(colnames(lateStageSmall)%in%ensembleBosses,colnames(lateStageSmall),NA),
     edge.width = 5*E(absGraph)$weight,
     vertex.label.color="black",
     vertex.label.cex=0.8,
     vertex.label.dist=0,
     vertex.label.degree=pi/2,
     vertex.color=optThetaComm$membership, #ifelse(colnames(lateStageSmall)%in%ensembleBosses,"black",optThetaComm$membership),
     edge.color = "gray50")

dev.off()