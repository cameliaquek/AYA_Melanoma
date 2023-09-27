#######################################
# YIM score
#######################################

# Load library
library(singscore)
library(GSEABase)
library(dplyr)
library(ggplot2)


# Load data
aya <- read.table("~/aya28_singscore_df.txt", header = 1, row.names = 1)
gUP <- read.table("upreg_YIM_genes.txt", header = 1)[,1]
gDN <- read.table("downreg_YIM_genes.txt", header = 1)[,1]


# Gene ranking
rankData <- rankGenes(aya, tiesMethod = "min")


# Given the ranked data and gene signature, simpleScore returns the scores and dispersions for each sample
scoredf <- simpleScore(rankData, upSet = gUP, downSet = gDN)


# Save plot
png("AYAsingscore_dispersion.png")
plotDispersion(scoredf,annot = ayaAnnot$Type,isInteractive = FALSE)
dev.off()
