#######################################                
#       pathway analysis
#######################################
#load library
library(enrichplot)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)

setwd("/Users/xinyubai/Documents/01_Research projects/AYA_immunotherapy/02analysis/01NR_G1G2_analysis/DE_pathway_analysis")

dat <- read.table("filtered_DEgenes_aya_G2vsG1.txt", header = T)
dat <- read.table("filtered_DE_AYA_CRvsPD.txt", header = T)
dat <- subset(dat, dat$padj<0.05) #for G1 vs G2
dat <- subset(dat, dat$padj<0.1) #for CR vs PD due to small sample size
dat <- tibble::rownames_to_column(dat, "Gene")
names(dat)[names(dat) == 'gene'] <- 'Gene'

head(dat)
#Gene log2FoldChange     lfcSE      stat      pvalue        padj
#1   ACBD4      -2.283025 0.4157193 -5.491746 3.98000e-08 0.000021200
#2   ACOT1      -3.878062 1.0772814 -3.599860 3.18389e-04 0.025263130
#3    ACP2      -1.253017 0.3522745 -3.556934 3.75208e-04 0.027818237


#log2FC values start from column 2
geneList <- dat[,3]
#gene symbols start from column1
names(geneList) <- as.character(dat[,1])
#Arrange the list in decreasing order. Important for the next steps.
geneList <- sort(geneList, decreasing = TRUE)


#Check that the gene list is correct: first row is gene name and second row is expression values
head(geneList)
#      SI    OLFM3     CRB1     EYA1    DMBT1     PIGR 
#9.581688 7.069940 6.412627 6.101842 5.924687 5.619746 

keytypes(org.Hs.eg.db)
#Convert gene symbols to EntrezID
ids <- bitr(as.character(dat[,1]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
dat2 = dat[dat$Gene %in% dedup_ids$SYMBOL,] #$Gene is because of the first column name
dat2$Y = dedup_ids$ENTREZID
geneListID <- dat2[,3]
names(geneListID) <- dat2$Y
geneListID <-na.omit(geneListID)
geneListID <- sort(geneListID, decreasing = TRUE)


#Check that the gene list is correct: first row is gene name and second row is expression values
head(geneListID)
#Example:
#    6476   118427    23418     2138     1755     5284 
#9.581688 7.069940 6.412627 6.101842 5.924687 5.619746  


#If this is needed, we can define fold change greater than 1.5 as DEGs.
#Run below code.
#gene <- names(geneList)[abs(geneList) > 1.5]
#head(gene)


#mut sig database analysis. This analysis is similar to the GSEA software (https://www.gsea-msigdb.org/gsea/index.jsp)
#if ID is used: entrez_gene OR gene_symbol OR human_gene_symbol
#only select C2 Kegg pathways
m_df <- msigdbr(species = "Homo sapiens", category = "C2")
try1 <- m_df %>% filter(gs_subcat == "CP:KEGG")
m_t2g <- try1 %>% dplyr::select(gs_name, entrez_gene)
head(m_t2g)
# A tibble: 6 x 2
#  gs_name               entrez_gene
#  <chr>                       <int>
#1 KEGG_ABC_TRANSPORTERS          19
#2 KEGG_ABC_TRANSPORTERS       10349
#3 KEGG_ABC_TRANSPORTERS       26154


#Plot dotplot
#Set min and max gene size to be the same as GSEA analysis
em2 <- GSEA(geneListID, minGSSize = 3, maxGSSize = 500 ,pvalueCutoff = 0.05, TERM2GENE = m_t2g)
em2 <- GSEA(geneListID, minGSSize = 10, maxGSSize = 500 ,pvalueCutoff = 0.05, TERM2GENE = m_t2g) #for G1vsG2
em2 <- GSEA(geneListID, minGSSize = 10, maxGSSize = 400 ,pAdjustMethod = "none", pvalueCutoff = 0.05, TERM2GENE = m_t2g) #for CR vs PD
results <- setReadable(em2, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.table(results, 
            file="aya16_G2vsG1_KEGG_gsea.txt", 
            sep="\t", col.names = NA)

png("KEGG_msigdb_top-new.png",units="in", width = 6, height = 4, res = 300)
dotplot(em2, showCategory = 50, title = "GSEA Group1 vs Group2" , split=".sign", font.size=9) + facet_grid(.~.sign)
dev.off()
#dotplot(em2, showCategory = 50, title = "Enriched Pathways using MSigDB" , split=".sign") + facet_grid(.~.sign)


png("ayaG2vsG1_KEGG_networktry3.png", units="cm", width = 25, height = 20, res =500)
edox <- setReadable(em2, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
cnetplot(edox, node_label="category", 
         cex_label_category = 1.2,color_category='firebrick',
         categorySize="p.adjust", foldChange=geneListID) +
  scale_color_gradient2(name='log2FC', low="#02284f", mid = "lightblue", high = "white")
dev.off()

png("ayaG2vsG1_KEGG_networktry3.png", units="cm", width = 25, height = 20, res =500)
edox <- setReadable(em2, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
cnetplot(edox, node_label="gene", 
         cex_label_gene = 1.2,color_category='firebrick',
         categorySize="p.adjust", foldChange=geneListID) +
  scale_color_gradient2(name='log2FC', low="#02284f", mid = "lightblue", high = "white")
dev.off()

png("ayaG1vsG2_KEGG_circle.png", units="cm", width = 30, height = 30, res =500, pointsize = 3)
cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
dev.off
while (!is.null(dev.list()))  dev.off()

#heatmaps
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
#cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

png("ayaG2vsG1_KEGG_heatplot.png", units="cm", width = 50, height = 10, res = 500)
heatplot(edox, showCategory=5)
dev.off
while (!is.null(dev.list()))  dev.off()

kegg <- enrichKEGG(gene = gene,
               universe = names(geneListID),
               organism = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               minGSSize = 10, 
               maxGSSize = 500)

enrich.kegg.res <- setReadable(kegg, OrgDb = org.Hs.eg.db)


#### GO over-representation analysis
gene <- names(geneListID)[abs(geneListID) > 2]
length(gene)

# Entrez gene ID
head(gene)


# GO enrichment

kk <- enrichGO(gene = gene,
               universe = names(geneListID),
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05,
               ont = "all",
               readable = T)

head(kk) 


library(GOplot)
#########################################################################################
###change function circle_dat because go enrichment results genes are separated by '/'###
#########################################################################################
circle_dat <- function(terms, genes){
  
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), '/')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), '/')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    zsc <- c(zsc, sum(value) / sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}
#########################################################################################
#########################################################################################

# terms
kk2 <- as.data.frame(kk)
colnames(kk2) <- c("category", "ID", "term", "GeneRatio", "BgRatio", "pvalue", "adj_pval", "qvalue", "genes", "count")

# genes
dat3 <- as.data.frame(dat)
colnames(dat3) <- c("ID", "baseMean", "logFC", "lfcSE", "stat", "pvalue", "padj")

# construct GOplot data frame
circ <- circle_dat(kk2, dat3)
head(circ)

## draw bubble plot to visualise the terms
GOBubble(circ, title = 'GO terms', display = 'multiple', bg.col = T, labels = 3)  

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
png("GOenrich_bubbleplot.png", units="in", width = 12, height = 8, res = 300)
GOBubble(reduced_circ, labels = 2.8)
dev.off()


# Generate a circular visualization of the results of gene- annotation enrichment analysis
png("GOenrich_circleplot.png", units="in", width = 12, height = 8, res = 300)
GOCircle(circ)
dev.off()

# Select genes of interes
vol <- read.delim("DE_G2vsG1_volcano.txt", header=TRUE, na.strings=c("","NA"))
head(vol)

top50 <- filter(vol, top == "top50")
dim(top50)

geneIDs <- top50$genes

topgenes <- dat3[order(-abs(dat3$logFC)),]

top50 <- topgenes[1:50,]
head(top50)
top <- top50[,c(1,3)]

#### chord plot
# Generate the matrix with selected processes
chord <- chord_dat(circ, top, reduced_circ$term)

#Eliminate processes genes that contain < 5
which(rowMeans(chord) == 0)
which(colMeans(chord) == 0)

#Eliminate genes related to no processes
if(length(which(rowMeans(chord) == 0) != 0)){
  paste("one or more genes are not related any of the processes")
  chord <- chord[-which(rowMeans(chord) == 0),]
}else {
  paste("all genes are related to at least one process")
}

#Eliminate processes genes that contain no genes
if(length(which(colMeans(chord) == 0) != 0)){
  paste("one or more processes are not related any of the selected genes")
  chord <- chord[,-which(colMeans(chord) == 0)]
}else {
  paste("all processes have at least one related gene")
}

#Eliminate processes with < 5 genes
if(length(which(colSums(chord) < 5) != 0)){
  chord <- chord[,-which(colSums(chord) < 5)]
}else {
  paste("all processes have at least one related gene")
}

png("GOenrich_chordplot.png", units="in", width = 15, height = 15, res = 300)

GOChord(chord, 
        space = 0.02, 
        gene.order = 'logFC', 
        gene.space = 0.25, 
        gene.size = 2,
        nlfc = 1)
dev.off()









###C5 GO pathways

c5 <- msigdbr(species = "Homo sapiens", category = "C5")
test <- c5 %>% filter(gs_subcat == "GO:BP")
m_t3g <- test %>% dplyr::select(gs_name, entrez_gene)
head(m_t3g)

c5 <- msigdbr(species = "Homo sapiens", category = "C5")
test <- c5 %>% filter(gs_subcat == "GO:CC")
m_t3g <- test %>% dplyr::select(gs_name, entrez_gene)
head(m_t3g)

c5 <- msigdbr(species = "Homo sapiens", category = "C5")
test <- c5 %>% filter(gs_subcat == "GO:MF")
m_t3g <- test %>% dplyr::select(gs_name, entrez_gene)
head(m_t3g)

#C5:BP/CC/MF (for MF - minGSSize adjusted to 10)
em3 <- GSEA(geneListID, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = m_t3g)

#dotplot
dotplot(em3, showCategory = 20, title = "Enriched Pathways using MSigDB" , split=".sign", font.size = 8) + facet_grid(.~.sign)

#networkplot
edox <- setReadable(em3, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, categorySize="pvalue", foldChange=geneListID, showCategory = 6)

#write results
results1 <- setReadable(em3, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.table(results1, 
            file="aya16_GOBP_gsea.txt", 
            sep="\t", col.names = NA)
results1 <- setReadable(em3, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.table(results1, 
            file="aya16_GOCC_gsea.txt", 
            sep="\t", col.names = NA)
results1 <- setReadable(em3, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.table(results1, 
            file="aya16_GOMF_gsea.txt", 
            sep="\t", col.names = NA)

##### for CRvsPD YIM signature 
res <- setReadable(em2, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.table(res, 
            file="YIM_CRvsPD_KEGGpathways.txt", 
            sep="\t", col.names = NA)

res2 <- setReadable(kk, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.table(res2, 
            file="YIM_CRvsPD_GOoverrepresentation.txt", 
            sep="\t", col.names = NA)
