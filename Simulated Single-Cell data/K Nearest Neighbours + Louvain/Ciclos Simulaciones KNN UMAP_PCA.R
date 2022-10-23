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

#### CICLO 1.3 ####
# In this cycle we run the batches simulation with two types of batches, 
# reduce the dimension with UMAP and make clusters with Louvain. 
# We keep the hyper parameters of UMAP, the vector
# of the cluster generated and the original clusters. At the end we will have the
# Jacard index to evaluate how the dimension reduction method works.

# Set the times to run the cycle
k <- 20
#Set the range of the components in the reduced dimension
a1 <- 13
b1 <- 14
n_comp <- seq(a1,b1)
# Set the range of neighbors
a2 <- 22
b2 <- 26
n_neig <- seq(a2,b2)
#Set the minimum distance
a3 <- 0.05
b3 <- 0.3
#Desired length of the sequence
c3 <- 5
m_dis <- seq(a3,b3,length.out=c3) 

jaccard <- NULL
hyper <- matrix(0, nrow = 3, ncol = k)
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
    for (b in n_neig) {
      for (c in m_dis) {
        SIM_umap <- runUMAP(SIM, 
                            ncomponents=a,
                            n_neighbors=b,
                            min_dist=c,
                            metric="euclidean",
                            altexp = NULL, 
                            name = "UMAP")
        clust.louvain <- clusterCells(SIM_umap, use.dimred = "UMAP", BLUSPARAM=NNGraphParam(cluster.fun="louvain", k=80))
        jaccard_temp <- mclustcomp(as.numeric(as.vector(clust.louvain)), grupos,
                                   types = 'jaccard')
        if(jaccard_temp$scores > jaccard[i]){
          jaccard[i] <- jaccard_temp$scores
          hyper[,i] <- c(a,b,c)
          clusters_ori[,i] <- grupos
          clusters_gen[,i] <- clust.louvain
        }
      }
    }
  }
}

coords_red <- SIM_umap$data[,1:2]
colLabels(SIM_umap) <- clust.louvain
plot_red <- plotReducedDim(SIM_umap, "UMAP", colour_by="label")
#jaccard
#hyper
#head(clusters_gen)
#head(coords_red)
write.table(jaccard, file = "Jaccard_ind_13.csv")
write.table(hyper, file =  "Hyper_UMAP_13.csv")
write.table(clusters_ori, file = "Cluster_originales_13.csv")
write.table(clusters_gen, file = "Cluster_generados_13.csv")
write.table(coords_red, file =  "Coords_UMAP_13.csv")
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
a1 <- 14
b1 <- 16
n_comp <- seq(a1,b1)
# Set the range of neighbors
a2 <- 23
b2 <- 28
n_neig <- seq(a2,b2)
#Set the minimum distance
a3 <- 0.05
b3 <- 0.4
#Desired length of the sequence
c3 <- 5
m_dis <- seq(a3,b3,length.out=c3) 

jaccard1 <- NULL
hyper1 <- matrix(0, nrow = 3, ncol = k)
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
    for (b in n_neig) {
      for (c in m_dis) {
        SIM_umap <- runUMAP(SIM, 
                            ncomponents=a,
                            n_neighbors=b,
                            min_dist=c,
                            metric="euclidean",
                            altexp = NULL, 
                            name = "UMAP")
        clust.louvain <- clusterCells(SIM_umap, use.dimred = "UMAP", BLUSPARAM=NNGraphParam(cluster.fun="louvain", k=25))
        jaccard_temp <- mclustcomp(as.numeric(as.vector(clust.louvain)), grupos,
                                   types = 'jaccard')
        if(jaccard_temp$scores > jaccard1[i]){
          jaccard1[i] <- jaccard_temp$scores
          hyper1[,i] <- c(a,b,c)
          clusters_ori1[,i] <- grupos
          clusters_gen1[,i] <- clust.louvain
        }
      }
    }
  }
}

coords_red1 <- SIM_umap$data[,1:2]
colLabels(SIM_umap) <- clust.louvain
plot_red1 <- plotReducedDim(SIM_umap, "UMAP", colour_by="label")

#jaccard1
#hyper1
#head(clusters_gen1)
#head(coords_red1)
write.table(jaccard1, file = "Jaccard_ind_23.csv")
write.table(hyper1, file =  "Hyper_UMAP_23.csv")
write.table(clusters_ori1, file = "Cluster_originales_23.csv")
write.table(clusters_gen1, file = "Cluster_generados_23.csv")
write.table(coords_red1, file =  "Coords_UMAP_23.csv")
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
a1 <- 14
b1 <- 18
n_comp <- seq(a1,b1)
# Set the range of neighbors
a2 <- 24
b2 <- 28
n_neig <- seq(a2,b2)
#Set the minimum distance
a3 <- 0.05
b3 <- 0.3
#Desired length of the sequence
c3 <- 5
m_dis <- seq(a3,b3,length.out=c3) 

jaccard2 <- NULL
hyper2 <- matrix(0, nrow = 3, ncol = k)
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
    for (b in n_neig) {
      for (c in m_dis) {
        SIM_umap <- runUMAP(SIM, 
                            ncomponents=a,
                            n_neighbors=b,
                            min_dist=c,
                            metric="euclidean",
                            altexp = NULL, 
                            name = "UMAP")
        clust.louvain <- clusterCells(SIM_umap, use.dimred = "UMAP", BLUSPARAM=NNGraphParam(cluster.fun="louvain", k=85))
        jaccard_temp <- mclustcomp(as.numeric(as.vector(clust.louvain)), grupos,
                                   types = 'jaccard')
        if(jaccard_temp$scores > jaccard2[i]){
          jaccard2[i] <- jaccard_temp$scores
          hyper2[,i] <- c(a,b,c)
          clusters_ori2[,i] <- grupos
          clusters_gen2[,i] <- clust.louvain
        }
      }
    }
  }
}

coords_red2 <- SIM_umap$data[,1:2]
colLabels(SIM_umap) <- clust.louvain
plot_red2 <- plotReducedDim(SIM_umap, "UMAP", colour_by="label")
#jaccard2
#hyper2
#head(clusters_gen2)
#head(coords_red2)
write.table(jaccard2, file = "Jaccard_ind_33.csv")
write.table(hyper2, file =  "Hyper_UMAP_33.csv")
write.table(clusters_ori2, file = "Cluster_originales_33.csv")
write.table(clusters_gen2, file = "Cluster_generados_33.csv")
write.table(coords_red2, file =  "Coords_UMAP_33.csv")
png("plot_red3.png")
plot_red2
dev.off()
