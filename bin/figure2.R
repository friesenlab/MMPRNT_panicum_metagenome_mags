######################################################################################################################
# This is R script to make a bar phyloseq plots for the all data comparing nitrogen fertilization switchgrass        #
# based on site soil microbiomes using metagenomics and mOTU v2                                                      #                      
#                                                                                                                    #
######################################################################################################################  

# Paper title "Genome-resolved metagenomics reveals niche partitioning within the switchgrass rhizosphere microbiome"
# Published date - TBD
# Journal - Phytobiomes
# Figure 2
# Written by Dr. Richard Allen White III
# Date Nov 1st, 2019

#Load Libraries
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("reshape2"); packageVersion("reshape2")
library("grid"); packageVersion("grid")
library("gridExtra"); packageVersion("gridExtra")
library("ape"); packageVersion("ape")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("vegan"); packageVersion("vegan")
library("ggpubr");packageVersion("ggpubr")
library("DESeq2");packageVersion("DESeq2")

#set ggplot2 to no background
theme_set(theme_bw())

#load data
otus <- read.delim("mmpnrt_mags_mOTUsv2_OTUs.txt") 
tax <- read.delim("mmpnrt_mag_mOTUsv2_taxonomyids.txt") 
meta <- import_qiime_sample_data("mmpnrt_mag_metadata.txt") 

#####################
#Format for phyloseq#
#####################

#fix row.names in OTU and taxa table
rownames(otus) <- paste0("OTU", 1:nrow(otus))
otus <- otus[, -c(1)]
rownames(tax) <- paste0("OTU", 1:nrow(tax))
tax <- tax[, -c(1)]

#convert data.frames to matrix
otu_f <- as.matrix(otus)
tax_f <- as.matrix(tax)

#combine for phyloseq
OTU_m <- otu_table(otu_f, taxa_are_rows = T)
TAX_m <- tax_table(tax_f)

#print for notebook
head(OTU_m)
head(TAX_m) 

#merge data for phyloseq
mOTU <- phyloseq(OTU_m, TAX_m, meta)
print(mOTU)

#build random tree for phyloseq object for unfiltered mOTU
random_tree = rtree(ntaxa(mOTU), rooted=TRUE, tip.label=taxa_names(mOTU))
plot(random_tree)

#merge mOTUs with random tree
full_mOTUs <- merge_phyloseq(mOTU, random_tree)
print(full_mOTUs)

# mOTUs by Phyla 
glom <- tax_glom(full_mOTUs, taxrank = "phylum")
glom <- transform_sample_counts(glom, function(x) x / sum(x))
# Select top Phyla
glom.top <- prune_taxa(names(sort(taxa_sums(glom), TRUE))[0:10], glom)
# Melt for graphing and strip prefix
glom.top.melt <- psmelt(glom.top)
glom.top.melt <- arrange(glom.top.melt, Type, phylum)
glom.top.melt$Phylum <- gsub('p__', '', gsub('Cyanobacteria/Melainabacteria group', 'Cyanobacteria', glom.top.melt$phylum))
#glom.top.melt$Phylum <- gsub('p__', '', glom.top.melt$phylum)
#glom.top.melt$Phylum1 <- gsub('Cyanobacteria/', '', glom.top.melt$phylum)


glom.top.melt$Type = factor(glom.top.melt$Type, levels=c("Pre", "Post"))
#glom.top.melt$Type_f <- as.character(glom.top.melt$Type_f)

#glom.top.melt$Type_f = ifelse(is.na(glom.top.melt$Type_f), 'Post', glom.top.melt$Type_f)

# Graph
glom.top.gg <- ggplot(glom.top.melt, aes(x = Sample, y = Abundance, fill = Phylum))
glom.top.gg <- glom.top.gg + geom_bar(stat = "identity") +
  facet_grid(~Type, scales = "free_x") +
  labs(fill = "Phylum", y = "Abundance", title = "A") +
  theme(plot.margin=unit(c(5,5,-25,5), units = "pt")) +
  scale_fill_brewer(palette = "Spectral") +
  theme(strip.background = element_blank()) 
glom.top.gg <- glom.top.gg + theme(text = element_text(size=26), plot.title = element_text(hjust = 0, size=50))
print(glom.top.gg)


# mOTUs by class
glom1 <- tax_glom(full_mOTUs, taxrank = "class")
glom1 <- transform_sample_counts(glom1, function(x) x / sum(x))
# Select top class
glom.top1 <- prune_taxa(names(sort(taxa_sums(glom1), TRUE))[0:10], glom1)
# Melt for graphing and strip prefix
glom.top.melt1 <- psmelt(glom.top1)
glom.top.melt1 <- arrange(glom.top.melt1, Type, class)
glom.top.melt1$Class <- gsub('c__', '', glom.top.melt1$class)

glom.top.melt1$Type = factor(glom.top.melt1$Type, levels=c("Pre", "Post"))
# Graph
glom.top.gg1 <- ggplot(glom.top.melt1, aes(x = Sample, y = Abundance, fill = Class))
glom.top.gg1 <- glom.top.gg1 + geom_bar(stat = "identity") +
  facet_grid(~Type, scales = "free_x") +
  labs(fill = "Class", y = "Abundance", title = "B") +
  theme(plot.margin=unit(c(5,5,-25,5), units = "pt")) +
  scale_fill_brewer(palette = "Spectral") +
  theme(strip.background = element_blank()) 
glom.top.gg1 <- glom.top.gg1 + theme(text = element_text(size=26), plot.title = element_text(hjust = 0, size=50))
print(glom.top.gg1)

###################
###Write outputs###
###################

write.table(glom.top.melt1, "class_mOTU_abundances.txt", sep ='\t')
write.table(glom.top.melt, "phyla_mOTU_abundances.txt", sep ='\t')

##################
###Arrange plot###
##################

grid.arrange(glom.top.gg, glom.top.gg1, ncol=1)