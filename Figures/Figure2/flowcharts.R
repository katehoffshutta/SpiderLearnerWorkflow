library(igraph)

miniMat = matrix(c(c(1,1,1,1,0.1,0.1,0.1,0.1,0.1),
                 c(1,1,1,1,0,0,0,0,0.1),
                 c(1,1,1,1,0,0,0,0,0),
                 c(1,1,1,1,0,0,0,0,0),
                 c(0.1,0,0,0,1,1,1,0,0),
                 c(0.1,0,0,0,1,1,1,0,0),
                 c(0.1,0,0,0,1,1,1,0,0),
                 c(0.1,0,0,0,0,0,0,1,1),
                 c(0.1,0.1,0,0,0,0,0,1,1)),ncol=9,byrow=T)

miniGraph = graph_from_adjacency_matrix(miniMat,
                                         mode="undirected",
                                         weighted=T,
                                         diag=F)

set.seed(14)
jpeg("ggm.jpeg")
plot(miniGraph,
     layout = layout_with_gem,
     edge.color="black",
     edge.width=3,
     vertex.border.size=3,
     vertex.color="black",
     vertex.label=NA,
     vertex.size=20)
dev.off()

miniClust=cluster_walktrap(miniGraph,steps = 2) 

new_cols <- c("darkgreen", "red", "blue")[membership(miniClust)]

set.seed(14)
jpeg("ggmClust.jpeg")
plot(miniClust,miniGraph,
     layout = layout_with_gem,
     col = new_cols,
     edge.width=3,
     mark.border="black",
     mark.col=c("lightgreen","pink","lightblue"),
     edge.color="black",
     vertex.label=NA,
     vertex.label.color="black",
     vertex.size=20)
dev.off()


new_cols <- c("darkgreen", "red", "blue")[membership(miniClust)]
new_cols[c(1:3,5:6,8)]=0

set.seed(14)
jpeg("ggmClustHubs.jpeg")
plot(miniClust,miniGraph,
     layout = layout_with_gem,
     col = new_cols,
     edge.width=3,
     mark.border="black",
     mark.col=c("lightgreen","pink","lightblue"),
     edge.color="black",
     vertex.label=NA,
     vertex.label.color="black",
     vertex.size=20)
dev.off()
