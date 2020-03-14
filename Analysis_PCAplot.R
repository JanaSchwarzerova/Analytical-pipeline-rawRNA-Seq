#******************************************************************************************************************************************************************
# Script PCA plots
#__________________________________________________________________________________________________________________________________________________________________
# R version 3.6.1
# author: Jana Schwarzerová

#==================================================================================================================================================================
#Path:
setwd('C:/Users/...')

#Libraries: 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("DESeq2")   #https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
BiocManager::install('PCAtools') #https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

library(DESeq2)
library(PCAtools)
#==================================================================================================================================================================
## Replicates (A,B,C,D,E,F,G): 

# Points-time samples: 
time=as.factor(rep(c("T1","T2","T3","T4","T5","T6")))
samples <- DataFrame(rep(c("T1","T2","T3","T4","T5","T6"),7), time)
x <- round(read.csv2(file = "CountTableBCDEFG_new.csv"))
dds<-DESeqDataSetFromMatrix(countData=x,colData=samples, design =~ time)
dds = DESeq(dds)
rld <- rlog(dds, blind=FALSE)
Pca <- DESeq2::plotPCA(rld, intgroup = 'time')
names = c("A1","A2","A3","A4","A5","A6",
          "B1","B2","B3","B4","B5","B6","C1","C2","C3","C4","C5","C6",
          "D1","D2","D3","D4","D5","D6","E1","E2","E3","E4","E5","E6",
          "F1","F2","F3","F4","F5","F6","G1","G2","G3","G4","G5","G6")
PCA <- Pca + geom_text(aes(label = names), position = position_nudge(y = 1))
PCA

#Divided to Butanol shock transcriptome & Standart cultivation transcriptome
Pca_table <- DESeq2::plotPCA(rld, intgroup = 'time', return=TRUE)
PCA_table <- Pca_table #Created new variable
PCA_table[,3] <- matrix("",nrow = nrow(Pca_table),ncol = 1) 
i <- 1   #initialization i
j <- 1   #initialization j
for (i in 1:nrow(Pca_table)) {
  for (j in 1:length(time)){
    #Standart cultivation vol1
    if(Pca_table[i,5]==paste("A",".",toString(j),".bam",sep="") || Pca_table[i,5]==paste("B",".",toString(j),".bam",sep="") 
     || Pca_table[i,5]==paste("C",".",toString(j),".bam",sep="")) {
    PCA_table[i,3]= "Standard cultivation transcriptome";
    }
    #Standart cultivation vol2
    else if(Pca_table[i,5]==paste("D",".",toString(j),".bam",sep="") || Pca_table[i,5]==paste("E",".",toString(j),".bam",sep="")){ #;
    PCA_table[i,3]= "Standard cultivation transcriptome";
    }
    #Butanol shock 
    else if(Pca_table[i,5]==paste("F",".",toString(j),".bam",sep="") || Pca_table[i,5]==paste("G",".",toString(j),".bam",sep="")) { #;
    PCA_table[i,3]= "Butanol shock transcriptome";
    }
    else { #;
    }
  }
}

#Visualization: 
df_out <- as.data.frame(PCA_table)
Final_PCA<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
Final_PCA<-Final_PCA + geom_point()+ geom_text(aes(label = names)) + xlab(PCA[["labels"]][["x"]]) + ylab(PCA[["labels"]][["y"]])
Final_PCA

#save as svg: 
ggsave(file="Final_PCA.svg",width=10, height=8)



