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
