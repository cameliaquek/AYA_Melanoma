####################################################
##
##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
##   Version 1.0 08.07.2016
##   Needs packages ggplot2,grid,gridExtra
##
####################################################

library(ggplot2)
library(grid)
library(gridExtra)

## calculate Immunophenoscore
ipsmap<- function (x) {
	if (x<=0) {
		ips<-0
	} else {
		if (x>=3) {
		 ips<-10
		} else {
			ips<-round(x*10/3, digits=0)
		}
	}
	return(ips)
}

## Assign colors 
my_palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(n = 1000)
mapcolors<-function (x) {
		za<-NULL
		if (x>=20) {
			za=1000
		} else {
			if (x<=-20) {
				za=1
			} else {
				za=round(166.5*x+500.5,digits=0)
			}
		}
		return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
		za2<-NULL
		if (x>=10) {
			za2=1000
		} else {
			if (x<=-10) {
				za2=1
			} else {
				za2=round(249.75*x+500.5,digits=0)
			}
		}
		return(my_palette2[za2])
}


## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
gene_expression<-read.table("EXPR.txt",row.names=1,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
sample_names<-names(gene_expression)

## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
# For different 
IPSG<-read.table("IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))

IPS<-NULL
MHCI<-NULL
MHCII<-NULL
ICP<-NULL
SCP<-NULL
AZ<-NULL

# Gene names in expression file
GVEC<-row.names(gene_expression)
# Genes names in IPS genes file
VEC<-as.vector(IPSG$GENE)
# Match IPS genes with genes in expression file
ind<-which(is.na(match(VEC,GVEC)))
# List genes missing or differently named
MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if (length(MISSING_GENES)>0) {
	cat("differently named or missing genes: ",MISSING_GENES,"\n")
}
for (x in 1:length(ind)) {
  print(IPSG[ind,])
}

for (i in 1:length(sample_names)) {	
	GE<-gene_expression[[i]]
	mGE<-mean(GE)
	sGE<-sd(GE)
	Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
	W1<-IPSG$WEIGHT
	WEIGHT<-NULL
	MIG<-NULL
	k<-1
	for (gen in unique_ips_genes) {
		MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
		WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
		k<-k+1
	}
	WG<-MIG*WEIGHT
	MHCI[i]<-mean(WG[1:9])
	MHCII[i]<-mean(WG[10:14])
	ICP[i]<-mean(WG[15:24])
	SCP[i]<-mean(WG[25:30])
	#AZ[i]<-sum(MHCI[i],MHCII[i],ICP[i],SCP[i])
	#IPS[i]<-ipsmap(AZ[i])
	
	## Plot Immunophenogram
	data_a <- data.frame (start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30), end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40), y1=c(rep(2.6,30),rep(0.4,4)),y2=c(rep(5.6,30),rep(2.2,4)),z=c(MIG[c(15:30,10:14,1:9)],ICP[i],SCP[i],MHCII[i],MHCI[i]),vcol=c(unlist(lapply(MIG[c(15:30,10:14,1:9)],mapcolors)), unlist(lapply(c(ICP[i],SCP[i],MHCII[i],MHCI[i]),mapbw))), label = c(unique_ips_genes[c(15:30,10:14,1:9)],"ICP","SCP","MHCII","MHCI"))
	data_a$label <- factor(data_a$label, levels=unique(data_a$label))
	plot_a1<-ggplot() + geom_rect(data=data_a, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +  coord_polar() + scale_y_continuous(limits = c(0, 6)) + scale_fill_manual(values =as.vector(data_a$vcol),guide=FALSE) + theme_bw() + theme(panel.margin = unit(0, 'mm'), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "white"), axis.text=element_blank(), axis.ticks= element_blank()) + geom_text(aes(x=5, y=1.3, label="ICP"), size=4) + geom_text(aes(x=15, y=1.3, label="SCP"), size=4) + geom_text(aes(x=25, y=1.3, label="MHC-II"), size=4) + geom_text(aes(x=35, y=1.3, label="MHC-I"), size=4)
	plot_a2<-plot_a1+geom_text(aes(x=10.25, y=4.1, label="ICOS +"), angle=85.5, size=4)+geom_text(aes(x=11.75, y=4.1, label="GITR +"),angle=76.25, size=4)+geom_text(aes(x=12.25, y=4.1, label="OX40 +"), angle=63.75,size=4)+geom_text(aes(x=14.75, y=4.1, label="CD137 +"), angle=51.25,size=4)+geom_text(aes(x=17.5, y=4.1, label="CD27 +"), angle=47.5,size=4)+geom_text(aes(x=12.5, y=4.1, label="CD28 +"), angle=32.5,size=4)		
	plot_a3<-plot_a2+geom_text(aes(x=20.5, y=4.1, label="PD-1 -"), angle=85.5, size=4)+geom_text(aes(x=21.5, y=4.1, label="CTLA4 -"), angle=76.5, size=4)+geom_text(aes(x=22.5, y=4.1, label="LAG3 -"), angle=67.5, size=4)+geom_text(aes(x=23.5, y=4.1, label="TIGIT -"), angle=58.5, size=4)+geom_text(aes(x=24.5, y=4.1, label="TIM3 -"), angle=49.5, size=4)+geom_text(aes(x=25.5, y=4.1, label="PD-L1 -"), angle=40.5, size=4)+geom_text(aes(x=26.5, y=4.1, label="PD-L2 -"), angle=31.5, size=4)+geom_text(aes(x=27.5, y=4.1, label="CD27 +"), angle=22.5, size=4)+geom_text(aes(x=28.5, y=4.1, label="ICOS +"), angle=13.5, size=4)+geom_text(aes(x=29.5, y=4.1, label="IDO1 -"), angle=4.5, size=4)
	plot_a4<-plot_a3+geom_text(aes(x=30.5, y=4.1, label="B2M +"), angle=-4.5, size=4)+geom_text(aes(x=31.5, y=4.1, label="TAP1 +"), angle=-13.5, size=4)+geom_text(aes(x=32.5, y=4.1, label="TAP2 +"), angle=-22.5, size=4)+geom_text(aes(x=33.5, y=4.1, label="HLA-A +"), angle=-31.5, size=4)+geom_text(aes(x=34.5, y=4.1, label="HLA-B +"), angle=-40.5, size=4)+geom_text(aes(x=35.5, y=4.1, label="HLA-C +"), angle=-49.5, size=4)+geom_text(aes(x=36.5, y=4.1, label="HLA-E +"), angle=-58.5, size=4)+geom_text(aes(x=37.5, y=4.1, label="HLA-F +"), angle=-67.5, size=4)
	plot_a5<-plot_a4+geom_text(aes(x=0, y=6, label=paste("Immunophenoscore: ",IPS[i],sep="")), angle=0,size=6,vjust=-0.5)+ theme(axis.title=element_blank())
	plot_a <-plot_a5 + theme(plot.margin=unit(c(0,0,0,0),"mm")) + geom_text(vjust=1.15,hjust=0,aes(x=25.5, y=6,label="\n\n\n\n   MHC-I: Major Histocompatibility Complex Class I                                 ICP: Inhibitory Checkpoints\n   MHCII: Major Histocompatibility Complex Class II              SCP: Stimulatory Checkpoints\n\n", hjust = 0), size=4)

	## Legend sample-wise (averaged) z-scores
	data_b <- data.frame (start = rep(0,23), end = rep(0.7,23), y1=seq(0,22,by=1), y2=seq(1,23,by=1),z=seq(-3,3,by=6/22),vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors))), label = LETTERS[1:23])
	data_b_ticks <- data.frame(x = rep(1.2, 7), value = seq(-3,3, by=1), y = seq(0,6, by=1)*(22/6) +0.5)
	legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"), panel.margin = unit(0,"null"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "white"), axis.text=element_blank(), axis.ticks= element_blank(), axis.title.x=element_blank())
	plot_b<-ggplot(hjust=0) + geom_rect(data=data_b, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) + scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + scale_fill_manual(values =as.vector(data_b$vcol),guide=FALSE) + geom_text(data=data_b_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +theme_bw() + legendtheme + ylab("Sample-wise (averaged) z-score") 
	
	## Legend weighted z-scores
	data_c <- data.frame (start = rep(0,23), end = rep(0.7,23), y1=seq(0,22,by=1), y2=seq(1,23,by=1),z=seq(-2,2,by=4/22),vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw))), label = LETTERS[1:23])
	data_c_ticks <- data.frame(x = rep(1.2, 5), value = seq(-2,2, by=1), y = seq(0,4, by=1)*(22/4) +0.5)
	plot_c<-ggplot() + geom_rect(data=data_c, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) + scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + scale_fill_manual(values =as.vector(data_c$vcol),guide=FALSE) + geom_text(data=data_c_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +theme_bw() + legendtheme + ylab("Weighted z-score") 
	
	## Save plot to file (1 pdf file for each sample)
	file_name<-paste("IPS_",sample_names[i],".pdf",sep="")
	pdf(file_name, width=10, height=8)
	grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
	dev.off()
}
DF<-data.frame(SAMPLE=sample_names,MHCI=MHCI,ICP=ICP,SCP=SCP,MHCII=MHCII,AZ=AZ,IPS=IPS)
write.table(DF,file="IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")
barplot(DF[[7]],DF[[1]],col=c(rep("yellow",7),rep("green",4)))