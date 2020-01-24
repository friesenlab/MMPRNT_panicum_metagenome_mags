###############################################################################################################################
# This is R script to compare MAG abundances the all data comparing fertilization switchgrass                                 #
# based on site soil microbiomes using metagenomics                                                                           #                      
#                                                                                                                             #
###############################################################################################################################  

# Paper title "Genome-resolved metagenomics reveals niche partitioning within the switchgrass rhizosphere microbiome"
# Published date - TBD
# Journal Phytobiomes
# Figure 6
# Written by Dr. Richard Allen White III
# Date Nov 1st, 2019

#load libraries
library("ggplot2"); packageVersion("ggplot2") #version 3.2.1
library("grid"); packageVersion("grid") #version 3.6.2
library("gridExtra"); packageVersion("gridExtra") #version 2.3
library("DESeq2");packageVersion("DESeq2") #version 1.26.0
library("reshape2");packageVersion("reshape2") #1.4.3

###############
###Load data###
###############

phy <- read.delim("bin_phyla_count.txt")
adun <- read.delim("bin_abundances_combined.txt")

#melt 
phy.m <- melt(phy)
adun.m <- melt(adun)

adun.m$Sqrt.abundance <- sqrt(adun.m$value)
adun.m$log2 <- log2(adun.m$value)
adun.m$log10 <- log10(adun.m$value)
adun.m$scale <- scale(adun.m$value)


d$Team2 <- factor(d$Team1, c("Cowboys", "Giants", "Eagles", "Redskins"))
# re-order the levels in the order of appearance in the data.frame
adun.m$no <- factor(adun.m$no, c("mag_I1", "mag_I2", "mag_I3", "mag_I4", "mag_I5",
                                 "mag_I6", "mag_I7", "mag_I8", "mag_I9", "mag_I10",
                                 "mag_I11","mag_I12","mag_I13","mag_I14","mag_P1",
                                 "mag_P2", "mag_P3", "mag_P4", "mag_P5", "mag_P6",
                                 "mag_P7", "mag_P8", "mag_P9", "mag_P10", "mag_P11",
                                 "mag_P12","mag_P13","mag_P14","mag_P15"))

#ggplot2 barplot values side by side
A <- ggplot(phy.m,aes(x=reorder(phyla,-value), y=value)) 
A <- A + geom_bar(stat="identity",position="dodge")
A <- A + theme_classic(base_family="Helvetica")
A <- A + labs(x="GTDB Phyla", y="MAG Count", title="A")
A <- A + theme(plot.title = element_text(hjust = 0, size=50), 
                 axis.text.x = element_text(size = 25, angle=90),
                 axis.title.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 legend.title = element_text(size = 16),
                 legend.text = element_text(size = 16),
                 axis.title.y = element_text(size = 25))
print(A)


B <- ggplot(adun.m, mapping = aes(x = no, y = variable, fill = log2)) 
B <- B + geom_tile() 
B <- B + theme_classic(base_family="helvetica") 
B <- B + labs(x="MAGs", y="Sample plots", title="B")
B <- B + theme(plot.title = element_text(hjust = 0, size=50), 
               axis.text.x = element_text(size = 25, angle=90),
               axis.title.x = element_text(size = 25),
               axis.text.y = element_text(size = 25),
               legend.title = element_text(size = 16),
               legend.text = element_text(size = 16),
               axis.title.y = element_text(size = 25))
B <- B + scale_fill_gradient(name = "Log2", low = "white", high = "darkred", na.value="yellow")
print(B)

##################
###Arrange plot###
##################

grid.arrange(A, B, ncol=2)


