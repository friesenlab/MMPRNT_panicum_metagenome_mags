###############################################################################################################################
# This is R script to compare functional contig orf count data in DESeq2 for the all data comparing fertilization switchgrass #
# based on site soil microbiomes using metagenomics                                                                           #                      
# CAZy DESeq contig level and barplot/MDS                                                                                     #
###############################################################################################################################  

# Paper title "Genome-resolved metagenomics reveals niche partitioning within the switchgrass rhizosphere microbiome"
# Published date - TBD
# Journal - Phytobiomes
# Figure 5
# Written by Dr. Richard Allen White III
# Date Nov 1st, 2019

#load libraries
library("ggplot2"); packageVersion("ggplot2") #version 3.2.1
library("grid"); packageVersion("grid") #version 3.6.2
library("gridExtra"); packageVersion("gridExtra") #version 2.3
library("DESeq2");packageVersion("DESeq2") #version 1.26.0
library("reshape2");packageVersion("reshape2") #version 1.4.3

###############
###Load Data### 
###############

countData <- read.csv("CAZy_contig_gene_id_counts.csv",header=TRUE)
rownames(countData) = countData[,1]
countData <- countData[,-1]

countData <- read.csv("CAZy_contig_funtaxa.csv",header=TRUE)
rownames(countData) = countData[,1]
countData <- countData[,-1]

#Load MetaData from Samples
colData <- read.delim("mmpnrt_mag_metadata.txt",header=T)
rownames(colData) = colData[,1]

################
###DESeq2 run###
################

#Paired by sample block 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~Block + Type)
dds <- DESeq(dds, fitType="mean")
res <- results(dds, pAdjustMethod="none")
summary(res)
res = res[order(res$padj,na.last=NA),]

#Write table sigtab at 0.1
alpha = 0.1
sigtab = res[(res$padj<alpha),]
write.table(sigtab, "CAZy_DeSeq2_sigtab_contigs_paired_new.txt", sep = "\t", row.names=T) #pval <0.10

#Write table sigtab at 0.05
alpha = 0.05
sigtab = res[(res$padj<alpha),]
write.table(sigtab, "CAZy_DeSeq2_sigtab_contigs_paired_pval0.05_new.txt", sep = "\t", row.names=T) #pval <0.05

#Write table sigtab at 0.01
alpha = 0.01
sigtab = res[(res$padj<alpha),]
write.table(sigtab, "CAZy_DeSeq2_sigtab_contigs_paired_pval0.01_new.txt", sep = "\t", row.names=T) #pval <0.01

#Adjust Variance Stabilizing Transformation of the data frame
vsd <- varianceStabilizingTransformation(dds)

################
####plot A######
################

A <- plotPCA(vsd, intgroup=c("Type", "Block"), returnData=TRUE)
percentVar <- round(100 * attr(A, "percentVar"))
A <- ggplot(A, aes(PC1, PC2, color=Type, shape=Block)) 
A <- A + geom_point(size=9) 
A <- A + labs(x=paste0("PC1: ",percentVar[1],"% variance"), y=paste0("PC2: ",percentVar[2],"% variance"), title="A")
A <- A + theme_classic(base_family="Helvetica")
A <- A + scale_colour_brewer(palette = "Set2")
A <- A + scale_fill_brewer(palette = "Set2")
A <- A + theme(plot.title = element_text(hjust = 0, size=50), 
               axis.text.x = element_text(size = 25),
               axis.title.x = element_text(size = 25),
               axis.text.y = element_text(size = 25),
               legend.title = element_text(size = 16),
               legend.text = element_text(size = 16),
               axis.title.y = element_text(size = 25))
print(A)

###############
###load data###
###############

cazy_sig <- read.delim("Annotated_CAZy_DeSeq2_sigtab_contigs_paired_pval0.05_new.txt")

################
####plot B######
################

B <- ggplot(cazy_sig, aes(x = reorder(funtaxa, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0))
B <- B + geom_bar(stat = "identity")
B <- B + coord_flip()
B <- B + theme_classic(base_family="Helvetica")
B <- B + labs(x="CAZy number w/taxonomy", y="DESeq2 - Fold Change (log2, pval <0.05)", title="B")
B <- B + theme(plot.title = element_text(hjust = 0, size=50), 
               axis.text.x = element_text(size = 25),
               axis.title.x = element_text(size = 25),
               axis.text.y = element_text(size = 25),
               legend.position = "none",
               axis.title.y = element_text(size = 25))
B <- B + scale_fill_manual(values = c("#CC6666","#660099"))
print(B)

##################
###Arrange plot###
##################

grid.arrange(A,B, ncol=2)




