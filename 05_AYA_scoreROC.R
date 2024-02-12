devtools::install_github('DavisLaboratory/singscore')

setwd("/Users/xinyubai/Documents/01_Research projects/AYA_immunotherapy/02analysis/01NR_G1G2_analysis/Gene-signatures_singscore/tcga_singscore")
library(singscore)
library(GSEABase)
library(dplyr)
library(ggplot2)


aya <- read.table("aya28_singscore_df.txt", header = 1, row.names = 1)
#aya$AYA4 <- NULL
#aya$AYA5 <- NULL
#aya$AYA6 <- NULL
#aya$AYA28 <- NULL
head(aya,6)

aya<- read.table("aya47_lognorm.txt", header = 1, row.names = 1)
tcga <- aya [,1:19]

gUP <- read.table("upreg_YIM_genes.txt", header = 1)[,1]
gDN <- read.table("downreg_YIM_genes.txt", header = 1)[,1]


# The recommended method for dealing with ties in ranking is 'min', you can
# change by specifying 'tiesMethod' parameter for rankGenes function.
rankData <- rankGenes(aya, tiesMethod = "min")
rankData <- rankGenes(tcga, tiesMethod = "min")

# Given the ranked data and gene signature, simpleScore returns the scores and 
# dispersions for each sample
scoredf <- simpleScore(rankData, upSet = gUP, downSet = gDN)
write.table(scoredf, 
            file="aya28_YIM_singscore2.txt", 
            sep="\t", col.names = NA)

write.table(scoredf,
            file = "tcga_YIM_singscore",
            sep="\t", col.names = NA)



# Note that, when only one gene set is available in a gene signature, one can 
# only input values for the upSet argument. In addition, a knownDirection 
# argument can be set to FALSE if the direction of the gene set is unknown.
ipres <- read.table("ipres.txt",header = T)
ipres <- unique(ipres)
ipres <- ipres$gene
scoredf_ipres <- simpleScore(rankData, upSet = ipres, knownDirection = FALSE)


write.table(scoredf_ipres, 
            file="aya28_IPRESsingscore.txt", 
            sep="\t", col.names = NA)
###
write.table(scoredf_ipres, 
            file="tcga_IPRESsingscore.txt", 
            sep="\t", col.names = NA)

immuUP <- read.table("immuUP.txt", header = 1)[,1]
immuDN <- read.table("immuDN.txt", header = 1)[,1]

scoredf_immu <- simpleScore(rankData, upSet = immuUP, downSet = immuDN)
##
write.table(scoredf_immu, 
            file="aya28_IMMUsingscore.txt", 
            sep="\t", col.names = NA)
##
write.table(scoredf_immu, 
            file="tcga_IMMUsingscore.txt", 
            sep="\t", col.names = NA)

#  Get the annotations of samples by their sample names
ayaAnnot <- read.table("aya_Annot.txt", header = 1)
# Sample annotations
ayaAnnot$Type

png("AYAsingscore_dispersion.png")
plotDispersion(scoredf,annot = ayaAnnot$Type,isInteractive = FALSE)
dev.off

scoreLS <- plotScoreLandscape(scoredf, scoredf_ipres, 
                   scorenames = c('IPRES', 'G1G2_DE'),hexMin = 8)



#######################################
# ROC analysis
#######################################

# Load library
library(pROC)

# Load and prepare data
dat <- read.table("~/aya28_all_scores.txt", na.strings = "NA")
roc1 <- roc(dat$Response, dat$YIM_Score,
            levels=c("CR_PR", "PD")) 
roc2 <- roc(dat$Response, dat$IPRES_Score,
            levels=c("PD", "CR_PR"))
roc3 <- roc(dat$Response, dat$IMMU_Score,
            levels=c("CR_PR", "PD"))
roc1; roc2; roc3 ##get AUC values



#########################################################

# Save plots
png("AYA28gene-signatures_ROC_suppl.png", unit = "cm", width = 15, height = 15, res = 500)
rocobj1 <- plot.roc(dat$Response, dat$YIM_Score,
                    #main="Gene signatures for prediction of AYA poor immunogenicity",
                    percent=TRUE,
                    col="#1c61b6",
                    print.auc=TRUE,
                    print.auc.x=20)

rocobj2 <- plot.roc(dat$Response, dat$IPRES_Score, 
                    percent=TRUE, 
                    col="#fc8d62",
                    print.auc=TRUE,
                    print.auc.x=20, print.auc.y=40,
                    add = T)

rocobj3 <- plot.roc(dat$Response, dat$IMMU_Score, 
                    percent=TRUE, 
                    col="#008600",
                    print.auc=TRUE,  
                    print.auc.x=20, print.auc.y=30,
                    add = T)

legend("bottomright", 
       legend=c("YIM score", "IPRES score", "IMMU_score"), 
       col=c("#fc8d62", "#1c61b6", "#008600"), 
       lwd=3)

dev.off()
