## Pre-processing raw data RNA-Seq 
## *******************************************************************************************************************
## F-replicate
## quality assesment (QA) raw dat ************************************************************************************
##____________________________________________________________________________________________________________________
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF

# Add module fastQC 
module add fastQC-0.11.5
# Take all files "*.gz" and is done quality check; results will be write in to "raw_data_qa"
fastqc -o raw_data_qa *.gz 
# Add multiQC
module add python36-modules-gcc
pip freeze | grep network
#networkx==2.0
# Go to raw_data_qa file 
cd raw_data_qa
# run Multiqc
multiqc . 
# add necessarilly modules 
export LC_ALL=C.UTF-8
export LANG=C.UTF-8 
# run Multiqc
multiqc . 

cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF
gunzip *.gz

## Delete rRNA ******************************************************************************************************
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/sortmerna-2.1-linux-64

#samples: F01, F02, F03, F04, F05, F06
for i in 1 2 3 4 5 6
do
./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\
./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db\
 --reads /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/raw_F0${i}.fastq\
 --aligned /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/sortmeRNA_results/aligned/F_${i} --fastx\
 --other /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/sortmeRNA_results/non_RNA/F_${i}_non_RNA --log -v -a 10 -m 4096
done 

## quality assesment nonrRNA sekvenci *******************************************************************************
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/sortmeRNA_results/non_RNA
# Take all files "*.fasta" and is done quality check; results will be write in to "non_RNA_qa"
fastqc -o /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/non_RNA_qa *.fastq
# Go to non_RNA_qa file 
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/non_RNA_qa
# run Multiqc
multiqc . 

## trimming dat *****************************************************************************************************
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/sortmeRNA_results/non_RNA
module add trimmomatic-0.36

for i in 1 2 3 4 5 6
do
java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar SE -threads 10 F_${i}_non_RNA.fastq F_${i}_non_RNA_trim.fq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

## quality assesment (QA) nonrRNAtrim *******************************************************************************

cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/sortmeRNA_results/non_RNA
# Take all files "*..fasta" and is done quality check; results will be write in to "non_RNA_trim_qa"
fastqc -o /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/non_RNA_trim_qa *.fq
# Go to non_RNA_trim_qa file 
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/non_RNA_trim_qa
# run Multiqc
multiqc . 

## Mapping sequences *************************************************************************************************

# Add module STAR
module add star-2.5.2b

## ___________________________ These was done in Pre-processing_A.sh 
# Go to genome_annotation file 
# cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation
# module add cufflinks-2.2.1
# gffread -E -O -T genome.gff3 -o genome.gtf
## Index genome ... ONLY ONCE!
# STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation --genomeFastaFiles /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation/genome.fasta --sjdbGTFfile /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation/genome.gtf --sjdbOverhang 48 --sjdbGTFfeatureExon CDS

# Itself mapping to index genome
for i in 1 2 3 4 5 6
do
STAR --runThreadN 10 --genomeDir /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation --readFilesIn /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/sortmeRNA_results/non_RNA/F_${i}_non_RNA_trim.fq --outFileNamePrefix /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/index_genome/F${i} --outFilterMultimapNmax 5 --outReadsUnmapped Fastx
done

## Quality Assesment (QA) mapping **********************************************************************************

# Go to Mapping_sequences file
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/index_genome

# run Multiqc
multiqc . 

## Sort SAM to BAM  *************************************************************************************************
## before building count table is necessary SAM sorts
module add samtools-1.4

for i in 1 2 3 4 5 6
do
samtools sort -l 9 -o F_${i}.bam F${i}Aligned.out.sam
done

## Creating count table  *********************************************************************************************
module add subread-1.5.2
subread-buildindex
subread-align --help

featureCounts -T 2 -a /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation/genome.gff3 -o /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repF/Count_table_F.txt -t gene -g locus_tag -O F_1.bam F_2.bam F_3.bam F_4.bam F_5.bam F_6.bam