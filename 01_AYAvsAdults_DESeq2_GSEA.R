###################################################            
# Differential gene expression and pathway analysis
###################################################

# Load library
library(DESeq2)
library(affy)
library(DEGreport)
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggrepel)
library(dplyr)
library(enrichplot)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

# Load data
dat <- read.table("~/aya_adults_counts.txt")
mdata <- read.table("aya_conditions_new.txt", row.names=1)
ayaNR <- row.names(mdata)
dat <- dat %>% select(all_of(ayaNR))
dat <- as.data.frame(dat)
mdata <- as.data.frame(mdata)
all(row.names(mdata) == colnames(dat)) #	[1] TRUE


# Set up deseq object and filter low gene counts
dds <- DESeqDataSetFromMatrix(countData=dat, colData=mdata, design=~IPgroup)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 20) >=5
dds_filter <- dds[idx,]
dds_filter <- assay(dds_filter)


# Create deseq2 data frame using the updated dds with genes being filtered
dds_filter <- DESeqDataSetFromMatrix(countData=dds_filter, 
                                     colData=mdata, 
                                     design=~IPgroup)


# Differential gene expression analysis Group 2 (G2) vs Group 1(G1)
dds2 <-DESeq(dds_filter)
res <- results(dds2, contrast=c("IPgroup","Group2", "Group1"))
na.omit(res) -> res #remove NA from the data


# Get sig gene lists
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

sig_genes <- sig_res %>% 
  pull(gene)

length(sig_genes)


# Plot Volcano plot for significant genes
vol <- read.delim(sig_genes, header=TRUE, na.strings=c("","NA"), row.names = 1)
q <- ggplot(vol, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=-log10(pvalue)), size = 2) + scale_color_viridis_c(option="B", direction=-1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size=45),
    axis.title.y = element_text(size=45),
    legend.text=element_text(size=30),
    legend.title=element_text(size=30))

png("volplot_aya_300DPI.png", pointsize=15, width=25, height=40, res=300, units="cm")
q + geom_text_repel(data=filter(vol, padj < 6.23987020526904E-09), aes(label=gene))
dev.off()

        
# Set up input file for pathway analysis
de <- read.table(sig_genes, row.names=NULL,col.names = 1)
degene <- de$X1
dat <- subset(dat, rownames(dat) %in% degene)
dat <- tibble::rownames_to_column(dat, "Gene")
dat <- vol
geneList <- dat[,3]
names(geneList) <- as.character(dat[,1])
geneList <- sort(geneList, decreasing = TRUE)


# Convert gene symbols to EntrezID
keytypes(org.Hs.eg.db)
ids <- bitr(as.character(dat[,1]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
dat2 = dat[dat$row.names %in% dedup_ids$SYMBOL,] #$Gene is because of the first column name
dat2$Y = dedup_ids$ENTREZID
geneListID <- dat2[,3]
names(geneListID) <- dat2$Y
geneListID <-na.omit(geneListID)
geneListID <- sort(geneListID, decreasing = TRUE)


# Select C2 Kegg pathways retrieved from Broad Institute MSigDB
m_df <- msigdbr(species = "Homo sapiens", category = "C2")
try <- m_df %>% filter(gs_subcat == "CP:KEGG")
m_t2g <- try %>% dplyr::select(gs_name, entrez_gene)


# Pathway results
de <- geneListID[abs(geneList) > 2]
de <- sort(de, decreasing = TRUE)
em2 <- GSEA(geneListID, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = m_t2g)

de <- geneListID[abs(geneList) > 1.5]
de <- sort(de, decreasing = TRUE)
em3 <- GSEA(de, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = m_t2g)
dotplot(em3, showCategory = 50, title = "Enriched Pathways using MSigDB" , split=".sign", font.size=9) + facet_grid(.~.sign)
edox2 <- setReadable(em3, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
p3 <- heatplot(edox2, foldChange=geneList)

cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

kk2 <- gseKEGG(geneList     = geneListID,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

kk2<-setReadable(kk2, OrgDb="org.Hs.eg.db", keyType = "ENTREZID")

# Save plots
png("KEGG_pathways1.png", width = 900, height = 900)
edox <- setReadable(em2, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
cnetplot(edox, categorySize="p.adjust", foldChange=geneListID)
dev.off()

png("KEGG_heatplot.png", width = 900, height = 500)
p2
dev.off()

png("KEGG_msigdb_top10.png")
dotplot(em2, showCategory = 50, title = "Enriched Pathways using MSigDB" , split=".sign", font.size=9) + facet_grid(.~.sign)
dev.off()

write.table(kk2, 
            file="KEGG_enrichedpathways.txt", 
            sep="\t",
            row.names = F)