#******************************************************************************************************************************************************************
# UPGMA – Cluster's methods applied to all replicates and evaluation using UPGMA
#__________________________________________________________________________________________________________________________________________________________________
# R version 3.6.1
# author: Jana Schwarzerová

#==================================================================================================================================================================
#Path:
setwd('C:/Users/....')

#Libraries: 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("phangorn")      #https://cran.r-project.org/web/packages/phangorn/phangorn.pdf
BiocManager::install("DESeq2")   #https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
BiocManager::install('PCAtools') #https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

library("DESeq2")
library("phangorn")
library("PCAtools")
#==================================================================================================================================================================

#Loading data:  – Replicates A,B,C,D,E,F,G: 
time = c("T1","T2","T3","T4","T5","T6")
samples <- DataFrame(rep(c("T1","T2","T3","T4","T5","T6"),7), time)

#Normalisation using DESeq2 & rlog
x <- round(read.csv2(file = "CountTableABCDEFG.csv"))
dds<-DESeqDataSetFromMatrix(countData=x,colData=samples, design =~ time)
dds = DESeq(dds)
#rld <- assay(rlog(dds,blind = FALSE))
rld <- rlog(dds,blind = FALSE)

#Computing UPGMA and visualization using UPGMA method
Pca <- DESeq2::plotPCA(rld, intgroup = 'time')
PC <- cbind(Pca[["data"]][["PC1"]],Pca[["data"]][["PC2"]])
rownames(PC)<- c("A1","A2","A3","A4","A5","A6",
                 "B1","B2","B3","B4","B5","B6",
                 "C1","C2","C3","C4","C5","C6",
                 "D1","D2","D3","D4","D5","D6",
                 "E1","E2","E3","E4","E5","E6",
                 "F1","F2","F3","F4","F5","F6",
                 "G1","G2","G3","G4","G5","G6")
Dist_Matrix <- dist(PC,method = "euclidean",diag = TRUE) # https://stat.ethz.ch/R-manual/R-patched/library/stats/html/dist.html
UPGMA <- upgma(Dist_Matrix,method ="average")
plot(UPGMA,show.tip.label = TRUE,use.edge.length = TRUE,node.depth = 2)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)
title(xlab ="Euclidean distance between PC1 and PC2")
