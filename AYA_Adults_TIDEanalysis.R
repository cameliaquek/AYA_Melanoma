setwd("/Users/xinyubai/Documents/01_Research projects/AYA_immunotherapy/02analysis/TIDE_analysis")

#get normalised counts (log2)
Data <- read.table("tcga_aya_adults_log2norm.txt")

dim(Data)
#[1] 18743   120

head(Data, 3)

## all adult counts

## all aya counts
aya <- Data[,1:47]
head(aya)
write.table(aya, 
            file="aya47_lognorm.txt", 
            sep="\t", col.names = NA)

##tcga aya counts
tcga <- Data[,1:19]
head(tcga)
write.table(tcga, 
            file="tcga19_lognorm.txt", 
            sep="\t", col.names = NA)


##mia aya counts
mia <- Data[,20:47]
head(mia)
dim(mia)
write.table(mia, 
            file="mia28_lognorm.txt", 
            sep="\t", col.names = NA)



## write function for normalising with row mean
# center with 'colMeans()'
center_rowmeans <- function(x) {
  xcenter = rowMeans(x)
  x - rep(xcenter, rep.int(ncol(x), nrow(x)))
}


# apply it
res <- center_rowmeans(Data)
write.table(res, 
            file="aya.adult_lognorm_rowcentered-log2FC.txt", 
            sep="\t", col.names = NA)

#
aya_res <- center_rowmeans(aya)

write.table(aya_res, 
            file="aya47_lognorm_rowcentered-log2FC.txt", 
            sep="\t", col.names = NA)

#
tcga_res <- center_rowmeans(tcga)

write.table(tcga_res, 
            file="tcga19_lognorm_rowcentered-log2FC.txt", 
            sep="\t", col.names = NA)

#
mia_res <- center_rowmeans(mia)

write.table(mia_res, 
            file="mia28_lognorm_rowcentered-log2FC.txt", 
            sep="\t", col.names = NA)


################################################################################
### correlation Treg vs T.exclusion score
library("ggpubr")

## all aya - pos corr Treg and exclusion
dat <- read.table("corr_Treg-T.exclusion.txt", header = TRUE) 
head(dat)
dim(dat)

png("aya47_TregvsT.Exclusion.png", units="in", width = 3, height = 3, res = 300)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat, x = "Treg", y = "T.exclusion_score", size = 2, color="#f1b8af",
          add = "reg.line", conf.int = TRUE, cor.coef.coord = c(0,2.0),
          add.params = list(color = "#f1b8af", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman", 
          cor.coef.size = 5,
          xlab = "Treg_CIBERSORTproportion", ylab = "T.exclusion_score")
dev.off()




##################################################################################
## mia aya
#setwd("/Users/xinyubai/Desktop/TIDE_analysis")
dat2 <- read.table("corr_IHCTreg-T.exclusion.txt", header = TRUE) 
head(dat2)
dim(dat2)

############## plot treg vs dysfunction score ################################
png("total.treg_dysfunction.png", unit = "cm", width = 12, height = 12, res = 500)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat2, x = "log2total.treg", y = "T.dysfunction_score", size = 2, color="black",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 6,
          xlab = "log2[mIF total Treg cells/mm2]", ylab = "T cell dysfunction score")
dev.off()

png("peri.treg_dysfunction.png", unit = "cm", width = 12, height = 12, res = 500)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat2, x = "log2peri.treg", y = "T.dysfunction_score", size = 2, color="black",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 6,
          xlab = "log2[mIF peritumoural Treg cells/mm2]", ylab = "T cell dysfunction score")
dev.off()

png("tumour.treg_dysfunction.png", unit = "cm", width = 12, height = 12, res = 500)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat2, x = "log2tumour.treg", y = "T.dysfunction_score", size = 2, color="black",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 6,
          xlab = "log2[mIF intratumoural Treg cells/mm2]", ylab = "T cell dysfunction score")
dev.off()


#########################################################################################################
##################higher Treg cell density is associated with lower exclusion scores#####################
png("mia_total.TregvsT.Exclusion.png", units="in", width = 3, height = 3, res = 300)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat2, x = "log2total.treg", y = "T.exclusion_score", size = 2, color="#e98d7f",
          add = "reg.line", conf.int = TRUE, cor.coef.coord = c(0,2.0),
          add.params = list(color = "#e98d7f", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman", 
          cor.coef.size = 5,
          xlab = "total.Treg_IHCdensity", ylab = "T.exclusion_score")
dev.off()

png("mia_tumour.TregvsT.Exclusion.png", units="in", width = 3, height = 3, res = 300)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat2, x = "Treg.tumour.density", y = "T.exclusion_score", size = 2, color="#e98d7f",
          add = "reg.line", conf.int = TRUE, cor.coef.coord = c(0,2.0),
          add.params = list(color = "#e98d7f", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman", 
          cor.coef.size = 5,
          xlab = "tumour.Treg_IHCdensity", ylab = "T.exclusion_score")
dev.off()

png("mia_peri.TregvsT.Exclusion.png", units="in", width = 3, height = 3, res = 300)
setwd("/Users/xinyubai/Desktop/TIDE_analysis/all_corr_plots")
ggscatter(dat2, x = "Treg.peritumour.density", y = "T.exclusion_score", size = 2, color="#e98d7f",
          add = "reg.line", conf.int = TRUE, cor.coef.coord = c(0,2.0),
          add.params = list(color = "#e98d7f", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman", 
          cor.coef.size = 5,
          xlab = "peri.Treg_IHCdensity", ylab = "T.exclusion_score")
dev.off()
