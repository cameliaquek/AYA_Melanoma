#######################################
# Variant identification and analysis
#######################################

# Load library
library(maftools)
library(dplyr)
library(stringr)


# Load and prepare data
var.annovar.maf <- annovarToMaf(annovar = "annoVar_aya24.txt", refBuild = "hg19", sep ="\t", table = "refGene", tsbCol = 'SAMPLE', MAFobj = FALSE)
var.annovar.maf$Tumor_Sample_Barcode <- var.annovar.maf$sampleID
var.annovar.maf$Tumor_Sample_Barcode <- gsub("\\_.*","",var.annovar.maf$Tumor_Sample_Barcode)
var.annovar.maf$Tumor_Sample_Barcode <- str_remove(var.annovar.maf$Tumor_Sample_Barcode, pattern = "^.*?(_|-)")
write.table(var.annovar.maf, file="aya24_maf.txt", sep="\t", quote=F, row.names = F)

# Load clinical information (optional)
aya.clin <- read.delim('aya24_clinical.txt', na.strings =  "NA" )

# Read the MAF files, with clinical data 
vc_nonsyn = c("Frame_Shift_Del", 
              "Frame_Shift_Ins", 
              "Splice_Site", 
              "Translation_Start_Site",
              "Nonsense_Mutation", 
              "Nonstop_Mutation", 
              "In_Frame_Del",
              "In_Frame_Ins", 
              "Missense_Mutation")

aya <- read.maf(maf = "aya24.maf", clinicalData = aya.clin,
                vc_nonSyn = c(vc_nonsyn, "5'Flank", "5'UTR","3'Flank","Intron"))


# Prepare legends for oncoprint plot
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  "5'Flank",
  "Missense_Mutation",
  "Nonsense_Mutation", 
  "5'UTR",
  "In_Frame_Del",
  "Frame_Shift_Del", 
  "Frame_Shift_Ins", 
  "Splice_Site"
)

aya_genes <- read.table("aya_variants.txt")

gcol = c("R_CR" = "#E38D7F",
         "R_PR" = "#414487",
         "NR_Group1" = "#7BA63C",
         "NR_Group2" = "#F2A516",
         "Unknown" = "light grey")
gcol = list(Response_group = gcol)


gorder = c("R_CR","R_PR","NR_Group1","NR_Group2","Unknown")

gcol = RColorBrewer::brewer.pal(n = 5,name = 'Spectral')
names(gcol) = c("R_CR","R_PR","NR_Group1","NR_Group2","Unknown")


# Save Oncoplot
png("Oncoprint_aya24_final2.png",  units = "cm", width = 20, height = 17, res = 300)
oncoplot(maf = aya, 
         genes = aya_genes,
         colors = vc_cols, 
         clinicalFeatures = c('Response_group' 
                              #'M_stage', 
                              #'Sex',
                              #'LDH',
                              #'StageIV_met',
                              #'X1st_line_immux',
                              #'Response'
                              ), 
         sortByAnnotation = TRUE,
         annotationOrder = gorder, 
         annotationColor = gcol,
         draw_titv = TRUE, 
         showTumorSampleBarcodes = TRUE,
         removeNonMutated = FALSE)
dev.off()

png("Oncoprint_aya24.png",  units = "cm", width = 20, height = 17, res = 300)
oncoplot(maf = aya,  colors = vc_cols,
         showTumorSampleBarcodes=TRUE, top=100, removeNonMutated = FALSE)
dev.off()


# Get summary of all your data
getSampleSummary(aya)
getGeneSummary(aya)
getFields(aya)


# Plot a summary of the file
png("mafsummary_aya.png", width = 600, height = 500, units="px", res=NA)
plotmafSummary(maf = aya, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


# Plot oncoplot for top x mutated genes
png("Oncoprint_aya24.png",  units = "cm", width = 15, height = 15, res = 300)
oncoplot(maf = aya, showTumorSampleBarcodes=TRUE, top=100, removeNonMutated = FALSE)
dev.off()
