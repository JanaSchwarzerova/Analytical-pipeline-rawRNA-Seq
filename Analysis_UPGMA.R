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
rld <- assay(rlog(dds,blind = FALSE))

#Computing UPGMA and visualization using UPGMA method
Dist_Matrix <- dist(t(rld),method = "manhattan",diag = TRUE) # https://stat.ethz.ch/R-manual/R-patched/library/stats/html/dist.html
UPGMA <- upgma(Dist_Matrix,method ="average")
plot(UPGMA)
