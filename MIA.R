library(R.matlab)
library(pheatmap)
mat <- readMat("MIA_example.mat")

ST <- mat[["ST.region.genes"]]
ST <- as.data.frame(ST)
colnames(ST) <- unlist(mat[["tissue.region.labels"]])             
rownames(ST) <- unlist(mat[["gene.names"]])


SC <- mat[["cell.type.genes"]]
SC <- as.data.frame(SC)
colnames(SC) <- unlist(mat[["cell.type.labels"]])
rownames(SC) <- unlist(mat[["gene.names"]])


M <- matrix(ncol = ncol(SC),nrow = ncol(ST))
M_enr <- matrix(ncol = ncol(SC),nrow = ncol(ST))
M_dep <- matrix(ncol = ncol(SC),nrow = ncol(ST))
geneList <- rownames(ST)
b <- length(geneList)


for (region in 1:ncol(ST)) {
  #region=1
  G <- rownames(ST[which(ST[,region]==1),])
  
  for (celltype in 1:ncol(SC)) {
     #celltype=2
     C <-  rownames(SC[which(SC[,celltype]==1),])
     a <- length(intersect(G,C))
    c <- length(G)
    d <- length(C)
    
    M_enr[region,celltype] <-  -1*log10(phyper(c-a-1,c,b-c,b-d))
    M[region,celltype] <- M_enr[region,celltype]
      M_dep[region,celltype] <- -1*log10(1-(phyper(c-a-1,c,b-c,b-d)))
      
      if (M[region,celltype] < M_dep[region,celltype]) {
       
        
         M[region,celltype] <-  -1*M_dep[region,celltype]
      }
    
      }  
  
}

M <- t(M)
rownames(M) <- colnames(SC)
colnames(M) <- colnames(ST)
M[M > 10] = 10
M[M < -10] = -10

pheatmap::pheatmap(M,scale = "none",cluster_cols = F,cluster_rows = F)
