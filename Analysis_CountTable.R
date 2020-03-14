#******************************************************************************************************************************************************************
# Script for creating Count Tables 
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

BiocManager::install("Rsubread") #https://bioconductor.org/packages/release/bioc/manuals/Rsubread/man/Rsubread.pdf

library(Rsubread)
#==================================================================================================================================================================

##CountTable from REPLICATE A: 
count_table_A =featureCounts(
  files=c("A_1.bam","A_2.bam", "A_3.bam","A_4.bam","A_5.bam","A_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_A[["counts"]], file = "CountTableA.csv")

##CountTable from REPLICATE B: 
count_table_B = featureCounts(
  files=c("B_1.bam","B_2.bam", "B_3.bam","B_4.bam","B_5.bam","B_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_B[["counts"]], file = "CountTableB.csv")

##CountTable from REPLICATE C: 
count_table_C = featureCounts(
  files=c("C_1.bam","C_2.bam", "C_3.bam","C_4.bam","C_5.bam","C_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_C[["counts"]], file = "CountTableC.csv")

##CountTable from REPLICATE D:
count_table_D = featureCounts(
  files=c("D_1.bam","D_2.bam", "D_3.bam","D_4.bam","D_5.bam","D_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_D[["counts"]], file = "CountTableD.csv")

##CountTable from REPLICATE E:
count_table_E = featureCounts(
  files=c("E_1.bam","E_2.bam", "E_3.bam","E_4.bam","E_5.bam","E_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_E[["counts"]], file = "CountTableE.csv")

##CountTable from REPLICATE F:
count_table_F = featureCounts(
  files=c("F_1.bam","F_2.bam", "F_3.bam","F_4.bam","F_5.bam","F_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_F[["counts"]], file = "CountTableF.csv")

##CountTable from REPLICATE G:
count_table_G = featureCounts(
  files=c("G_1.bam","G_2.bam", "G_3.bam","G_4.bam","G_5.bam","G_6.bam"),
  annot.ext="genome.gff3",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="gene",
  GTF.attrType ="locus_tag",
  useMetaFeatures=FALSE)
#Save as *.csv
write.csv(count_table_G[["counts"]], file = "CountTableG.csv")
