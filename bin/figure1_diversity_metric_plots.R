###############################################################################################################################
# This is R script to compare Alpha/Beta Diversity phyloseq plots for the all data comparing fertilization switchgrass        #
# based on site soil microbiomes using metagenomics and mOTU v2                                                               #                      
# Barplots                                                                                                                    #
###############################################################################################################################  

# Paper title "Genome-resolved metagenomics reveals niche partitioning within the switchgrass rhizosphere microbiome"
# Published date - TBD
# Journal - Phytobiomes
# Figure 1
# Written by Dr. Richard Allen White III
# Date November 1st, 2019

#Load Libraries
library("phyloseq"); packageVersion("phyloseq") #version 1.30.0
library("ggplot2"); packageVersion("ggplot2") #version 3.2.1
library("reshape2"); packageVersion("reshape2") #version 1.4.3
library("grid"); packageVersion("grid") #version 3.6.2
library("gridExtra"); packageVersion("gridExtra") #version 2.3
library("ape"); packageVersion("ape") #version 5.3
library("plyr"); packageVersion("plyr") #version 1.8.5
library("dplyr"); packageVersion("dplyr") #version 0.8.3
library("vegan"); packageVersion("vegan") #version 2.5.6
library("ggpubr");packageVersion("ggpubr") #version 0.2.4
library("DESeq2");packageVersion("DESeq2") #version 1.26.0

#set ggplot2 to no background
theme_set(theme_bw())

#load data
otus <- read.delim("mmpnrt_mags_mOTUsv2_OTUs.txt") 
tax <- read.delim("mmpnrt_mag_mOTUsv2_taxonomyids.txt") 
meta <- import_qiime_sample_data("mmpnrt_mag_metadata.txt") 

###################################
###Function for extended metrics###
###################################

# This function extends phyloseq::estimate_richness() function by 
# implimenting two evenness metrics.
# See this PR https://github.com/joey711/phyloseq/pull/575

estimate_richness_mod <- function(physeq, split=TRUE, measures=NULL){
  
  if( !any(otu_table(physeq)==1) ){
    # Check for singletons, and then warning if they are missing.
    # These metrics only really meaningful if singletons are included.
    warning(
      "The data you have provided does not have\n",
      "any singletons. This is highly suspicious. Results of richness\n",
      "estimates (for example) are probably unreliable, or wrong, if you have already\n",
      "trimmed low-abundance taxa from the data.\n",
      "\n",
      "We recommended that you find the un-trimmed data and retry."
    )
  }
  
  # If we are not splitting sample-wise, sum the species. Else, enforce orientation.
  if( !split ){
    OTU <- taxa_sums(physeq)		
  } else if( split ){
    OTU <- as(otu_table(physeq), "matrix")
    if( taxa_are_rows(physeq) ){ OTU <- t(OTU) }
  }
  
  # Define renaming vector:
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Pielou", "Simpson", "InvSimpson", "SimpsonE", "Fisher")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "pielou", "simpson", "invsimpson", "simpsone", "fisher")
  # If measures was not explicitly provided (is NULL), set to all supported methods
  if( is.null(measures) ){
    measures = as.character(renamevec)
  }
  # Rename measures if they are in the old-style
  if( any(measures %in% names(renamevec)) ){
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
  }
  
  # Stop with error if no measures are supported
  if( !any(measures %in% renamevec) ){
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  
  # Initialize to NULL
  outlist = vector("list")
  # Some standard diversity indices
  estimRmeas = c("Chao1", "Observed", "ACE")
  if( any(estimRmeas %in% measures) ){ 
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
  }
  if( "Shannon" %in% measures ){
    outlist <- c(outlist, list(shannon = diversity(OTU, index="shannon")))
  }
  if( "Pielou" %in% measures){
    #print("Starting Pielou")
    outlist <- c(outlist, list(pielou = diversity(OTU, index = "shannon")/log(estimateR(OTU)["S.obs",])))
  }
  if( "Simpson" %in% measures ){
    outlist <- c(outlist, list(simpson = diversity(OTU, index="simpson")))
  }
  if( "InvSimpson" %in% measures ){
    outlist <- c(outlist, list(invsimpson = diversity(OTU, index="invsimpson")))
  }
  if( "SimpsonE" %in% measures ){
    #print("Starting SimpsonE")
    outlist <- c(outlist, list(simpsone = diversity(OTU, index="invsimpson")/estimateR(OTU)["S.obs",]))
  }
  if( "Fisher" %in% measures ){
    fisher = tryCatch(fisher.alpha(OTU, se=TRUE), 
                      warning=function(w){
                        warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
                        suppressWarnings(fisher.alpha(OTU, se=TRUE)[, c("alpha", "se")])
                      }
    )
    if(!is.null(dim(fisher))){
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    } else {
      outlist <- c(outlist, Fisher=list(fisher))
    }
  }
  out = do.call("cbind", outlist)
  # Rename columns per renamevec
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  # Final prune to just those columns related to "measures". Use grep.
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), ignore.case=TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop=FALSE]
  # Make sure that you return a data.frame for reliable performance.
  out <- as.data.frame(out)
  return(out)
}

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

######################
#Plot alpha diversity#
######################

#using unrarefied data
metrics <- c("Observed", "ACE", "Shannon","SimpsonE")
rich.mOTUs <- estimate_richness_mod(full_mOTUs, measures = metrics)

# merge richness with metadata
DF.mOTUs <- merge(rich.mOTUs, sample_data(full_mOTUs), by = 0)
DF.mOTUs %>% head
mdf.mOTUs = melt(DF.mOTUs, measure.vars = metrics)
mdf.mOTUs %>% head

#plot richness box plot treatment - closed boxes
A <- ggplot(mdf.mOTUs, aes(Type, value)) + 
  geom_boxplot(aes(fill = Type), outlier.size = 0) +
  geom_jitter(width = 0.3) +
  facet_wrap(facets = ~variable, nrow = 1, scales = "free_y") +
  labs(x = "Category", y = "Alpha Diversity Metric", title = "A") +
  scale_colour_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 1),
        strip.background = element_blank(),
        axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0, size=50),
        text = element_text(size=26))
A <- A + scale_fill_brewer(palette = "Set2")
print(A)

######################
#Plot beta diversity#
######################

#beta diversity - Weighted Unifrac
ordu = ordinate(full_mOTUs, "PCoA", "unifrac", weighted=T)
B <- plot_ordination(full_mOTUs, ordu, color="Type", shape="Block") + geom_point(size=8)
#plot first to obtain axis values then add
B <- B + labs(x="PC1 [28.4%]", y="PC2 [16.3%]", title="B")
B <- B + theme(text = element_text(size=26), plot.title = element_text(hjust = 0, size=50))
B <- B + scale_colour_brewer(palette = "Set2")
B <- B + scale_fill_brewer(palette = "Set2")
print(B)

#beta diversity - Unweighted Unifrac
ordu1 = ordinate(full_mOTUs, "PCoA", "unifrac", weighted=F)
C <- plot_ordination(full_mOTUs, ordu1, color="Type", shape="Block") + geom_point(size=8)
#plot first to obtain axis values then add
C <- C + labs(x="PC1 [24.5%]", y="PC2 [16.7%]", title="C")
C <- C + theme(text = element_text(size=26), plot.title = element_text(hjust = 0, size=50))
C <- C + scale_colour_brewer(palette = "Set2")
C <- C + scale_fill_brewer(palette = "Set2")
print(C)

#beta diversity - Unweighted Unifrac
ordu2 = ordinate(full_mOTUs, "PCoA", "bray")
D <- plot_ordination(full_mOTUs, ordu2, color="Type", shape="Block") + geom_point(size=8)
#plot first to obtain axis values then add
D <- D + labs(x="PC1 [32.5%]", y="PC2 [14.7%]", title="D")
D <- D + theme(text = element_text(size=26), plot.title = element_text(hjust = 0, size=50))
D <- D + scale_colour_brewer(palette = "Set2")
D <- D + scale_fill_brewer(palette = "Set2")
print(D)

##################
###Arrange plot###
##################

grid.arrange(A, arrangeGrob(B, C, D, ncol=3), heights=c(1/4, 1/4), ncol=1)

####################
###adonis testing###
####################

#by type unifrac unweighted
df = as(sample_data(full_mOTUs), "data.frame")
d = phyloseq::distance(full_mOTUs, "unifrac", weighted = F)
adonis(d ~ Block + Type, df, permutations = 999)

#by type unifrac weighted
df1 = as(sample_data(full_mOTUs), "data.frame")
d1 = phyloseq::distance(full_mOTUs, "unifrac", weighted = T)
adonis(d1 ~ Block + Type, df1, permutations = 999)

#by type brays
df2 = as(sample_data(full_mOTUs), "data.frame")
d2 = phyloseq::distance(full_mOTUs, "bray", paired=T, weighted=F)
adonis(d2 ~ Block + Type, df2, permutations = 999)

#by type brays
df3 = as(sample_data(full_mOTUs), "data.frame")
d3 = phyloseq::distance(full_mOTUs, "bray", paired=T, weighted=T)
adonis(d3 ~ Block + Type, df3, permutations = 999)




