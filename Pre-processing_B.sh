## Pre-processing raw data RNA-Seq 
## *******************************************************************************************************************
## B-replicate
## quality assesment (QA) raw dat ************************************************************************************
##____________________________________________________________________________________________________________________
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB

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

cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB
gunzip *.gz

## Delete rRNA ******************************************************************************************************
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/sortmerna-2.1-linux-64

#samples: B01, B02, B03, B04, B05, B06
for i in 1 2 3 4 5 6
do
./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\
./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db\
 --reads /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/raw_B0${i}.fastq\
 --aligned /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/sortmeRNA_results/aligned/B_${i} --fastx\
 --other /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/sortmeRNA_results/non_RNA/B_${i}_non_RNA --log -v -a 10 -m 4096
done 

## quality assesment nonrRNA sekvenci *******************************************************************************
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/sortmeRNA_results/non_RNA
# Take all files "*.fasta" and is done quality check; results will be write in to "non_RNA_qa"
fastqc -o /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/non_RNA_qa *.fastq
# Go to non_RNA_qa file 
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/non_RNA_qa
# run Multiqc
multiqc . 

## trimming dat *****************************************************************************************************
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/sortmeRNA_results/non_RNA
module add trimmomatic-0.36

for i in 1 2 3 4 5 6
do
java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar SE -threads 10 B_${i}_non_RNA.fastq B_${i}_non_RNA_trim.fq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

## quality assesment (QA) nonrRNAtrim *******************************************************************************

cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/sortmeRNA_results/non_RNA
# Take all files "*..fasta" and is done quality check; results will be write in to "non_RNA_trim_qa"
fastqc -o /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/non_RNA_trim_qa *.fq
# Go to non_RNA_trim_qa file 
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/non_RNA_trim_qa
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
STAR --runThreadN 10 --genomeDir /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation --readFilesIn /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/sortmeRNA_results/non_RNA/B_${i}_non_RNA_trim.fq --outFileNamePrefix /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/index_genome/B${i} --outFilterMultimapNmax 5 --outReadsUnmapped Fastx
done

## Quality Assesment (QA) mapping **********************************************************************************

# Go to Mapping_sequences file
cd /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/index_genome

# run Multiqc
multiqc . 

## Sort SAM to BAM  *************************************************************************************************
## before building count table is necessary SAM sorts
module add samtools-1.4

for i in 1 2 3 4 5 6
do
samtools sort -l 9 -o B_${i}.bam B${i}Aligned.out.sam #sam tool index
done

## Creating count table  *********************************************************************************************
module add subread-1.5.2
subread-buildindex
subread-align --help

featureCounts -T 2 -a /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/genome_annotation/genome.gff3 -o /auto/brno6/home/xschwa16/RNASeq_C_beijerinckii_NRRL_589/RNASeq_repB/Count_table_B.txt -t gene -g locus_tag -O B_1.bam B_2.bam B_3.bam B_4.bam B_5.bam B_6.bam