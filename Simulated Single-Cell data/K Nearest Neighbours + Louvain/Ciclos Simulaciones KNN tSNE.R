#### Librerias y funciones #####
#Jaccard index (cindex)
#source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/clusterIndex.R") 
library(mclustcomp)
# RNA-seq
library("scater")
library("splatter")
library(scran)
library(bluster)
library(igraph)

### refs ####
# https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rclustering/rclustering/

#### CICLO 1.4 ####
# In this cycle we run the batches simulation with two types of batches, 
# reduce the dimension with UMAP and make clusters with Louvain. 
# We keep the hyper parameters of UMAP, the vector
# of the cluster generated and the original clusters. At the end we will have the
# Jacard index to evaluate how the dimension reduction method works.

# Set the times to run the cycle
k <- 20
#Set the range of the components in the reduced dimension
a1 <- 2
b1 <- 3
n_comp <- seq(a1,b1)

jaccard <- NULL
hyper <- NULL
clusters_ori <- matrix(0, nrow = 800, ncol = k)
clusters_gen <- matrix(0, nrow = 800, ncol = k)

for (i in 1:k) {
  SIM <- splatSimulateGroups(batchCells = c(400, 400),
                             de.facLoc = 0.1,
                             de.facScale = 0.4,
                             dropout.type = "experiment",
                             seed = i*29)
  SIM <- logNormCounts(SIM)
  jaccard[i] <- -50
  g <- colData(SIM)$Batch
  grupos <- NULL
  for (j in 1:length(g)) {
    if(g[j]=="Batch1"){
      grupos[j] <- 1
    } else { if(g[j]=="Batch2"){
      grupos[j] <- 2
    }
    }
  }
  for (a in n_comp) {
        SIM_tsne <- runTSNE(SIM, 
                            ncomponents=a,
                            altexp = NULL, 
                            name = "TSNE")
        clust.louvain <- clusterCells(SIM_tsne, use.dimred = "TSNE", BLUSPARAM=NNGraphParam(cluster.fun="louvain", k=80))
        jaccard_temp <- mclustcomp(as.numeric(as.vector(clust.louvain)), grupos,
                                   types = 'jaccard')
        if(jaccard_temp$scores > jaccard[i]){
          jaccard[i] <- jaccard_temp$scores
          hyper[i] <- a
          clusters_ori[,i] <- grupos
          clusters_gen[,i] <- clust.louvain
    }
  }
}

coords_red <- SIM_tsne$data[,1:2]
colLabels(SIM_tsne) <- clust.louvain
plot_red <- plotReducedDim(SIM_tsne, "TSNE", colour_by="label")
#jaccard
#hyper
#head(clusters_gen)
#head(coords_red)
write.table(jaccard, file = "Jaccard_tsne_13.csv")
write.table(hyper, file =  "Hyper_TSNE_13.csv")
write.table(clusters_ori, file = "Cluster_originales_tsne_13.csv")
write.table(clusters_gen, file = "Cluster_generados_tsne_13.csv")
write.table(coords_red, file =  "Coords_TSNE_13.csv")
png("plot_red13.png")
plot_red
dev.off()

#### CICLO 2.3 ####
# In this cycle we run the groups simulation with five types of groups, 
# reduce the dimension with UMAP and make clusters with Louvain. 
# We keep the hyper parameters of UMAP, the vector
# of the cluster generated and the original clusters. At the end we will have the
# Jacard index to evaluate how the dimension reduction method works.

# Set the times to run the cycle
k <- 20
#Set the range of the components in the reduced dimension
a1 <- 2
b1 <- 3
n_comp <- seq(a1,b1)

jaccard1 <- NULL
hyper1 <- NULL
clusters_ori1 <- matrix(0, nrow = 1000, ncol = k)
clusters_gen1 <- matrix(0, nrow = 1000, ncol = k)
for (i in 1:k) {
  SIM <- splatSimulateGroups(batchCells = 1000,
                             group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2),                      de.prob = c(0.05, 0.05, 0.08, 0.05, 0.06),                      de.facLoc = 0.1,
                             de.facScale = 0.4,
                             dropout.type = "experiment",
                             seed = i*29)
  SIM <- logNormCounts(SIM)
  jaccard1[i] <- -50
  g <- colData(SIM)$Group
  grupos <- NULL
  for (j in 1:length(g)) {
    if(g[j]=="Group1"){
      grupos[j] <- 1
    } else { if(g[j]=="Group2"){
      grupos[j] <- 2
    } else { if(g[j]=="Group3"){
      grupos[j] <- 3
    } else { if(g[j]=="Group4"){
      grupos[j] <- 4
    } else { if(g[j]=="Group5"){
      grupos[j] <- 5
    }    
    }
      
    }
      
    }
      
    }
  }
  for (a in n_comp) {
        SIM_tsne <- runTSNE(SIM, 
                            ncomponents=a,
                            altexp = NULL, 
                            name = "TSNE")
        clust.louvain <- clusterCells(SIM_tsne, use.dimred = "TSNE", BLUSPARAM=NNGraphParam(cluster.fun="louvain", k=25))
        jaccard_temp <- mclustcomp(as.numeric(as.vector(clust.louvain)), grupos,
                                   types = 'jaccard')
        if(jaccard_temp$scores > jaccard1[i]){
          jaccard1[i] <- jaccard_temp$scores
          hyper1[i] <- a
          clusters_ori1[,i] <- grupos
          clusters_gen1[,i] <- clust.louvain
    }
  }
}

coords_red1 <- SIM_tsne$data[,1:2]
colLabels(SIM_tsne) <- clust.louvain
plot_red1 <- plotReducedDim(SIM_tsne, "TSNE", colour_by="label")

#jaccard1
#hyper1
#head(clusters_gen1)
#head(coords_red1)
write.table(jaccard1, file = "Jaccard_tsne_23.csv")
write.table(hyper1, file =  "Hyper_TSNE_23.csv")
write.table(clusters_ori1, file = "Cluster_originales_tsne_23.csv")
write.table(clusters_gen1, file = "Cluster_generados_tsne_23.csv")
write.table(coords_red1, file =  "Coords_TSNE_23.csv")
png("plot_red2.png")
plot_red1
dev.off()


#### CICLO 3.3 ####
# In this cycle we run the paths simulation with three paths, 
# reduce the dimension with UMAP and make clusters with Louvain. 
# We keep the hyper parameters of UMAP, the vector
# of the cluster generated and the original clusters. At the end we will have the
# Jacard index to evaluate how the dimension reduction method works.

# Set the times to run the cycle
k <- 20
#Set the range of the components in the reduced dimension
a1 <- 2
b1 <- 3
n_comp <- seq(a1,b1)

jaccard2 <- NULL
hyper2 <- NULL
clusters_ori2 <- matrix(0, nrow = 500, ncol = k)
clusters_gen2 <- matrix(0, nrow = 500, ncol = k)
for (i in 1:k) {
  SIM <- splatSimulatePaths(batchCells = 500,
                            group.prob = c(0.33, 0.33, 0.34),
                            path.from = c(0, 1, 1),
                            de.facLoc = 0.3,
                            de.facScale = 0.6,
                            dropout.type = "none",
                            seed = i*29)
  SIM <- logNormCounts(SIM)
  jaccard2[i] <- -50
  g <- colData(SIM)$Group
  grupos <- NULL
  for (j in 1:length(g)) {
    if(g[j]=="Path1"){
      grupos[j] <- 1
    } else { if(g[j]=="Path2"){
      grupos[j] <- 2
    } else { if(g[j]=="Path3"){
      grupos[j] <- 3
    } 
    }
    }
  }
  for (a in n_comp) {
        SIM_tsne <- runTSNE(SIM, 
                            ncomponents=a,
                            altexp = NULL, 
                            name = "TSNE")
        clust.louvain <- clusterCells(SIM_tsne, use.dimred = "TSNE", BLUSPARAM=NNGraphParam(cluster.fun="louvain", k=85))
        jaccard_temp <- mclustcomp(as.numeric(as.vector(clust.louvain)), grupos,
                                   types = 'jaccard')
        if(jaccard_temp$scores > jaccard2[i]){
          jaccard2[i] <- jaccard_temp$scores
          hyper2[i] <- a
          clusters_ori2[,i] <- grupos
          clusters_gen2[,i] <- clust.louvain
    }
  }
}

coords_red2 <- SIM_tsne$data[,1:2]
colLabels(SIM_tsne) <- clust.louvain
plot_red2 <- plotReducedDim(SIM_tsne, "TSNE", colour_by="label")
#jaccard2
#hyper2
#head(clusters_gen2)
#head(coords_red2)
write.table(jaccard2, file = "Jaccard_tsne_33.csv")
write.table(hyper2, file =  "Hyper_TSNE_33.csv")
write.table(clusters_ori2, file = "Cluster_originales_tsne_33.csv")
write.table(clusters_gen2, file = "Cluster_generados_tsne_33.csv")
write.table(coords_red2, file =  "Coords_TSNE_33.csv")
png("plot_red3.png")
plot_red2
dev.off()
