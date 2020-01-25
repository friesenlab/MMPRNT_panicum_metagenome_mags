###############################################################################################################################
# This is R script to compare Alpha/Beta Diversity phyloseq plots for the all data comparing fertilization switchgrass        #
# based on site soil microbiomes using metagenomics]                                                                          #                      
# Barplots for bin statistics                                                                                                 #
###############################################################################################################################  

# Paper title "Genome-resolved metagenomics reveals niche partitioning within the switchgrass rhizosphere microbiome"
# Published date - TBD
# Journal - Phytobiomes
# Supplemental Figure 1
# Written by Dr. Richard Allen White III
# Date November 1st, 2019

#Load Libraries
library("ggplot2"); packageVersion("ggplot2") #version 3.2.1
library("reshape2"); packageVersion("reshape2") #version 1.4.3

#Load data
stats <- read.delim("bin_comparsion.txt")

#melt to format for ggplot2
mm <- melt(stats)

#ggplot2 barplot
c <- ggplot(mm, aes(x=reorder(stats,-value), value, fill=variable)) 
c <- c + geom_bar(stat="identity", position=position_dodge())
c <- c + theme_bw(base_size=15, base_family="helvetica") #removes background/sets up text size
c <- c + theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
c <- c + scale_fill_brewer(palette = "Set2")
c <- c + labs(x="", y="Count", title="")
c <- c + guides(fill=guide_legend(title=NULL))
c <- c + theme(legend.justification=c(1,0), legend.position=c(1,.75))
c <- c + theme(text = element_text(size=38))
c