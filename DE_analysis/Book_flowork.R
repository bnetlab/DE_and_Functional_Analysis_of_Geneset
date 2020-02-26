
#tutrial Libro

#obtiene todos lo datos de ENA

cut - f11 samples_at_ENA . txt | xargs wget 

# obterner muestras espacificas 
for ACC_NR in ERR458493 ERR458494 ERR458495 ERR458496 ERR458497 ERR458498ERR458499; do  egrep ${ACC_NR} ERP004763_sample_mapping.tsv | cut -f 11 PRJEB5348.txt | xargs wget; done

#fastqce analises 

../FastQC/fastqc ERR458506.fastq.gz --extract -o fastqc_results/

  #to see all the results together need run multiqc in the directory with all fastqc results 
multiqc -d fastqc_results/ -o QC
  # see the results report 
firefox QC/multiqc_report.html 


#Qaligment 

# Download genome sequence of S. cerevisiae from UCSC
 wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2 bit
#convert from bit  to fasta 
 ../UCSC_tools/twoBitToFa sacCer3.2bit sacCer3.fa

 #download the gtf and bed file 
 
 
 #download the table with all field from UCSU and process 
 
# Confirm the head

 head -n1 sacCer3_allfields.gtf
 
 cut -f 2- sacCer3_allfields.gtf | sed '1d' |\
 ../UCSC_tools/genePredToGtf file stdin sacCer3_Refseq.gtf

 head -n1 sacCer3_Refseq.gtf
 
 
 #Generated the genomic Index
 REF_DIR=~/Data_RNA-seq/Referece_Genome
 runSTAR=~/Data_RNA-seq/STAR-2.7.1a/bin/Linux_x86_64_static/STAR
 
 ${runSTAR} --runMode genomeGenerate --genomeDir STARindex/ --genomeFastaFiles ${REF_DIR}/sacCer3.fa   --sjdbGTFfile ${REF_DIR}/sacCer3.gtf --sjdbOverhang 49 --runThreadN 2
 
 #Aligment 
 #defines the files fastq routes 
 # This step has to be done for each individual FASTQ file
 FILES=`ls -m Data_raw/WT/*fastq.gz| sed 's/ //g'`
 FILES=`echo $FILES|sed 's/ //g'`

 ${runSTAR} --genomeDir STARindex --readFilesIn $FILES --readFilesCommand zcat --outFileNamePrefix alignment_STAR/WT_1_ --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 2

 #BAM file indexing Most downstream applications will require a .BAM.BAI file together with every BAM  file to quickly access the BAM files without having to load them into memory

 samtools index alignment_STAR/WT_6_Aligned.sortedByCoord.out.bam 


# a traves del archivo log.final se puede revisar la calidad del alimiento the unique reads is the more importan t

#to see the aligment quality we can used 
 
 samtools flagstat WT_1_Aligned.sortedByCoord.out.bam #use the sam tool to see how many reads were mapped
#RSeQC you can see the number of reads mapped
 
 bam_stat.py -i WT_1_Aligned.sortedByCoord.out.bam 

 #and analize the quantity of red mappes (unique )
 samtools flagstat WT_1_Aligned.sortedByCoord.out.bam > flagstat_WT_1.txt
 
 
 #See the gene aligment distibution 
 
 read_distribution.py -r /home/edianfranco/Data_RNA-seq/Referece_Genome/sacCer3.bed -i alignment_STAR/SNF2_1_Aligned.sortedByCoord.out.bam
 
 #Genes Body Gene body coverage To assess possible 3’ or 5’ biases, 
 
 geneBody_coverage.py  -i alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam -r /home/edianfranco/Data_RNA-seq/Referece_Genome/sacCer3.bed  -o GeneBody_coverage_WT_1
 
#Gene-based read counting###### 
 
 
# Quality control with QoRTs
 #count the read for run the program 
 for FASTQ in Data_raw/SNF2/ERR458500*gz; do zcat $FASTQ | wc -l ; done | paste -sd+ | bc | awk '{print $1/4}' 
 #run the prgram 
 java -Xmx4g -jar hartleys-QoRTs-099881f/QoRTs.jar QC --singleEnded --seqReadCt 1885330 --generatePdfReport alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam /home/edianfranco/Data_RNA-seq/Referece_Genome/sacCer3.gtf  ./QoRTs_output_WT_1
 

 #Read Quantification
#for this we use the tool subread  to make a raw count of reads  
 subread-1.6.4-source/bin/featureCounts -a /home/edianfranco/Data_RNA-seq/Referece_Genome/sacCer3_Refseq.gtf -o features_count_results.txt alignment_STAR/*bam
 #summary.txt a summa aboyt teh process and result content the coutn 
 
 #count the exons 
 
 subread-1.6.4-source/bin/featureCounts -a /home/edianfranco/Data_RNA-seq/Referece_Genome/sacCer3_Refseq.gtf -f -t exon -O -o features_counts_exons.txt alignment_STAR/*bam
#remove repetitive exons

 sort -k2,2n -k3,3n features_counts_exons.txt | uniq > features_counts_exons_unique.txt
 ##################R######################
 
 # read data from count_feature to meake de DE analises 
 library("magrittr")
 
read.counts2<- read.table("/home/edianfranco/Data_RNA-seq/features_count_results.txt", header = TRUE)
row.names(read.counts)<-read.counts$Geneid
read.counts<- read.counts[,-c(1:6)]
orig_names <- names(read.counts )
names(read.counts ) <- c("SNF2 _1", "SNF2 _2", "SNF2 _3", "SNF2 _4", "SNF2 _5", "SNF2 _6","SNF2 _7","WT_1", "WT_2", "WT_3", "WT_4", "WT_5","WT_6","WT_7")

#names(read.counts) <- gsub(".*(SNF2|WT)(_[0 -9]+) .*", "\\1\\2",orig_names) automatic name

#Now that we have the read counts, we also need some information about the samples, which will be stored in colData
sample_info <- data.frame(condition = gsub("_[0 -9]+", "", names (read.counts)),row.names = names (read.counts))

library(DESeq2)

#generate the DESeqDataSet

DESeq_dataset<-DESeqDataSetFromMatrix(countData = read.counts, colData = sample_info, design = ~ condition)


#Dataframe verfication
DESeq_dataset %>% head
assay(DESeq_dataset) %>% head
rowRanges(DESeq_dataset) %>% head
#test what counts () returns
counts(DESeq_dataset)%>%str

# remove genes without any counts

DESeq_dataset<-DESeq_dataset[rowSums(counts(DESeq_dataset))>0,]

colSums(counts(DESeq_dataset)) # should be the same as colSums ( readcounts )

# calculate the size factor and add it to the data set
DESeq_dataset<-estimateSizeFactors(DESeq_dataset)
sizeFactors(DESeq_dataset)

#normalized
count.sf_nomarlized<-counts(DESeq_dataset, normalized= TRUE)

#log transformation

log.norm.count<-log2(count.sf_nomarlized + 1)


#plot the results
par(mfrow=c(2,1))  #to plot the following two images underneath each other
# first , boxplots of non - transformed read counts (one per sample )
boxplot(count.sf_nomarlized, notch= TRUE, main= "Unstrafomed read counts", ylab="Read counst")

# box plots of log2 - transformed read counts
boxplot(log.norm.count, notch= TRUE, main= "log2-trnsformed read counts", ylab=" log2 (Read counst)")

#Transformation of read counts including variance shrinkage
# obtain regularized log - transformed values
DESeq_dS_rlog<-rlog(DESeq_dataset,blind = TRUE)
rlog_norm_count<-assay(DESeq_dS_rlog)


# mean -sd plot for rlog - transformed data
library(vsn)
library(ggplot2)
msd_plot<-meanSdPlot(rlog_norm_count,ranks = FALSE,plot =FALSE)

msd_plot$gg+ggtitle("rlog-tranformed read counts") + ylab ("Standard desviation")


#####Differential Gene Expression Analysis#######


#1-DESeq2 workflow

#DGE analysis is performed on the raw data
str (colData(DESeq_dataset)$condition)

# set WT as the first -level - factor

colData (DESeq_dataset)$condition <- relevel(colData(DESeq_dataset)$condition , "WT")

# sequencing depth normalization between the samples
DESeq_dataset_<-DESeq(DESeq_dataset)
# gene - wise dispersion estimates across all samples
DESeq_dataset_<-estimateSizeFactors(DESeq_dataset_)
# this fits a negative binomial GLM andapplies Wald statistics to each gene

DESeq_dataset_<-nbinomWaldTest(DESeq_dataset_)

#The results() function lets you extract the base means across samples
DGE_results<-results(DESeq_dataset_,independentFiltering = TRUE,alpha = 0.5)

#the DESeqResult object can basically be handled like a data . frame
head (DGE_results)
table (DGE_results$padj<0.05)

rownames(subset(DGE_results, padj<0.05))

#Exploratory plots following DGE analysis

hist (DGE_results$pvalue,col = " grey ", border = " white ", xlab = "", ylab = "",
      main = " frequencies of p- values ")

plotMA (DGE_results, alpha = 0.05 , main = "WT vs. SNF2 mutants ",
       ylim = c( -4 ,4))














