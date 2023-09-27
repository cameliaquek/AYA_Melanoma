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
# 先绘制1条ROC曲线
plot(roc1, 
     print.auc=TRUE, # 图像上输出AUC的值
     print.auc.x=0.4, print.auc.y=0.5, # 设置AUC值坐标为（x，y）
     #auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     #auc.polygon.col="#fff7f7",  # 设置ROC曲线下填充色
     #grid=c(0.5, 0.2), # 设置两轴网格线的间隔为0.5，0.2
     #grid.col=c("black", "black"),  # 设置两轴间隔线条的颜色
     #print.thres=TRUE, # 图像上输出最佳截断值
     main="Gene signatures for prediction of AYA poor immunogenicity",  # 添加图形标题
     col="#fc8d62",    # 设置ROC曲线颜色
     #legacy.axes=TRUE  # 使x轴从0到1，表示为1-特异度
     percent=TRUE)   

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
