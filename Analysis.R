# Load libraries
library(DESeq2)
library(ReactomePA)
library(EnrichmentBrowser)
library(gplots)
library(ggvenn) 
library(RColorBrewer)
library(genefilter)
library(clusterProfiler) 
library(DOSE)
library(org.Mm.eg.db) 
library(tidyverse)
library(scales)
library(ggpubr)

##### Data preparation #####

# Load data
counts <- list()
metadata <- list()
metadata$TGFB <- read.table("metadata_TGFB.txt", header = TRUE)
metadata$TOR <- read.table("metadata_TOR.txt", header = TRUE)
counts$TGFB <- read.table("counts_TGFB", header = TRUE)
counts$TOR <- read.table("counts_TOR", header = TRUE)

# Prepare table
rownames(counts$TGFB) <- counts$TGFB$Geneid
counts$TGFB <- counts$TGFB[,-1:-6]
rownames(counts$TOR) <- counts$TOR$Geneid
counts$TOR <- counts$TOR[,-1:-6]

# Concatenate data
counts$all <- cbind(counts$TOR, counts$TGFB)
metadata$all <- rbind(metadata$TOR, metadata$TGFB)

# Filtering Data
keep <- rowSums(counts$all >= 10) >= 3 # removes genes that have less than 10 counts in total in at least 3 samples
counts$filt <- counts$all[keep,]

# Create metatable
coldata <- list()
coldata$all <- cbind(metadata$all["ID"], metadata$all["Condition"])

# Reorganize columns
counts$filt = counts$filt[, coldata$all$ID]

# Create full tables for DESeq2
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = counts$filt, colData = coldata$all, design = ~Condition)

##### DESeq2 pipeline & Z-score transformation #####

# Run DESeq2 pipeline
dds <- DESeq(ddsFullCountTable)
rld <- rlog(dds)

# Z score transformation
rld_matrix <- assay(rld)
counts$Z <- as.data.frame(t(scale(t(rld_matrix))))
counts$Z$entrezID = mapIds(org.Mm.eg.db, keys = row.names(counts$Z), keytype = "ENSEMBL", column = "ENTREZID")
counts$Z <- counts$Z[!is.na(counts$Z$entrezID),]

##### Principal Component Analysis (PCA) #####

# TOR
coloring <- c("#FF2600","yellow","#FF7800","#FEBC81","#FF0064","#FF80B2","black","grey","#007369","#00B6A6","#8952FF","#C3A6FF")
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = coloring)
}
plotPCA(rld[,c("TOR_1", "TOR_2", "TOR_3","TOR_4","TOR_5","TOR_6","TOR_7", "TOR_8", "TOR_9","TOR_10","TOR_11","TOR_12","TOR_13", "TOR_14", "TOR_15","TOR_16","TOR_17","TOR_18","TOR_19", "TOR_20", "TOR_21","TOR_22","TOR_23","TOR_24","TOR_25","TOR_26","TOR_27","TOR_28","TOR_29","TOR_30","TOR_31", "TOR_32", "TOR_33","TOR_34","TOR_35","TOR_36")], intgroup = c("Condition"), returnData=T) %>%
  separate(col=Condition, into= c("Background", "Day", "VC"), remove = FALSE, sep = "_") %>%
  ggplot(., aes(x=PC1, y=PC2, color = Condition ))+
  geom_point(size=3)
ggsave("TOR_PCA.pdf", units="mm", dpi = 300, width = 170)

# TGFB
coloring <- c("#F83FF8", "#FF2600","#2BFA2C","#0431FF","#02FDFF","black")
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = coloring)
}
plotPCA(rld[,c("TGFB_1", "TGFB_2", "TGFB_3","TGFB_7", "TGFB_8", "TGFB_9","TGFB_10","TGFB_11","TGFB_12","TGFB_13", "TGFB_14", "TGFB_15","TGFB_16","TGFB_17","TGFB_18","TGFB_19", "TGFB_20", "TGFB_21")], intgroup = c("Condition"), returnData=T) %>%
  ggplot(aes(x=PC1, y=PC2, color=Condition))+
  geom_point(size=3)
ggsave("TGFB_PCA.pdf", units="mm", dpi = 300, width = 170)

##### TOR #####

##### Results D0-D2 mut ######

# Pair-wise comparison
res <- list()
res$mutd02 <- results(dds, contrast = c("Condition","mutant_day-2_noVC","mutant_control_noVC"))
write.csv(as.data.frame(res$mutd02), file="TOR_mutD0-2_res.csv")

# MA Plot
pdf("TOR_MA-plot_mutD0-2.pdf", width = 4, height = 4) 
plotMA(res$mutd02, alpha = 0.05, ylim = c(-14, 14), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
down <- list()
up <- list()
res_sig <- res$mutd02[which(res$mutd02$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$mutd02 <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_mutD0-2_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$mutd02 <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_mutD0-2_up.csv")

# Add entrez gene IDs
res$mutd02$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$mutd02), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$mutd02) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_mutD02 = as.data.frame(gsea)
write.csv(gsea_table_mutD02, file="TOR_GSEA_mutD0-2.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_mutD02 = as.data.frame(Reactome)
write.csv(reactome_table_mutD02, file="TOR_Reactome_mutD0-2.csv")

# Select pathways
gene_sets <- getGenesets(org = "mmu", db = "go", return.type="list")
pathways <- list(
  # DNA damage terms
  "GO:0006281_DNA_repair",
  "GO:0019985_translesion_synthesis",
  "GO:0006298_mismatch_repair",
  "GO:0043504_mitochondrial_DNA_repair",
  "GO:1901255_nucleotide-excision_repair_involved_in_interstrand_cross-link_repair",
  "GO:0000724_double-strand_break_repair_via_homologous_recombination",
  "GO:0006303_double-strand_break_repair_via_nonhomologous_end_joining",
  "GO:0097681_double-strand_break_repair_via_alternative_nonhomologous_end_joining",
  "GO:0006289_nucleotide-excision_repair",
  "GO:0006284_base-excision_repair",
  # other terms
  "GO:0007179_transforming_growth_factor_beta_receptor_signaling_pathway",
  "GO:0030511_positive_regulation_of_transforming_growth_factor_beta_receptor_signaling_pathway",
  "GO:0001837_epithelial_to_mesenchymal_transition",
  "GO:0006325_chromatin_organization",
  "GO:0061621_canonical_glycolysis",
  "GO:0000278_mitotic_cell_cycle"
)


# Pathways
path <- list()
path$mutd02 <- data.frame(pathway = character(),
                          entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$mutd02 <- rbind(path$mutd02, data.frame(pathway = v1, entrezID = v2))
}
path$mutd02 <- merge(path$mutd02, as.data.frame(res$mutd02), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$mutd02, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_mutD0-2.txt")

##### Results D0-D4 mut ######

# Pair-wise comparison
res$mutd04 <- results(dds, contrast = c("Condition","mutant_day-4_noVC","mutant_control_noVC"))
write.csv(as.data.frame(res$mutd04), file="TOR_mutD0-4_res.csv")

# MA Plot
pdf("TOR_MA-plot_mutD0-4.pdf", width = 4, height = 4) 
plotMA(res$mutd04, alpha = 0.05, ylim = c(-14, 14), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$mutd04[which(res$mutd04$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$mutd04 <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_mutD0-4_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$mutd04 <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_mutD0-4_up.csv")

# Add entrez gene IDs
res$mutd04$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$mutd04), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$mutd04) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_mutD04 = as.data.frame(gsea)
write.csv(gsea_table_mutD04, file="TOR_GSEA_mutD0-4.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_mutD04 = as.data.frame(Reactome)
write.csv(reactome_table_mutD04, file="TOR_Reactome_mutD0-4.csv")

# Pathways
path$mutd04 <- data.frame(pathway = character(),
                          entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$mutd04 <- rbind(path$mutd04, data.frame(pathway = v1, entrezID = v2))
}
path$mutd04 <- merge(path$mutd04, as.data.frame(res$mutd04), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$mutd04, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_mutD0-4.txt")

##### Results D0-D2 mut VC ######

# Pair-wise comparison
res$mutd02VC <- results(dds, contrast = c("Condition","mutant_day-2_VC","mutant_control_VC"))
write.csv(as.data.frame(res$mutd02VC), file="TOR_mutD0-2VC_res.csv")

# MA Plot
pdf("TOR_MA-plot_mutD0-2VC.pdf", width = 4, height = 4) 
plotMA(res$mutd02VC, alpha = 0.05, ylim = c(-14, 14), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$mutd02VC[which(res$mutd02VC$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$mutd02VC <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_mutD0-2VC_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$mutd02VC <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_mutD0-2VC_up.csv")

# Add entrez gene IDs
res$mutd02VC$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$mutd02VC), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$mutd02VC) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_mutD02VC = as.data.frame(gsea)
write.csv(gsea_table_mutD02VC, file="TOR_GSEA_mutD0-2VC.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_mutD02VC = as.data.frame(Reactome)
write.csv(reactome_table_mutD02VC, file="TOR_Reactome_mutD0-2VC.csv")

# Pathways
path$mutd02VC <- data.frame(pathway = character(),
                            entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$mutd02VC <- rbind(path$mutd02VC, data.frame(pathway = v1, entrezID = v2))
}
path$mutd02VC <- merge(path$mutd02VC, as.data.frame(res$mutd02VC), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$mutd02VC, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_mutD0-2VC.txt")

##### Results D0-D4 mut VC ######

# Pair-wise comparison
res$mutd04VC <- results(dds, contrast = c("Condition","mutant_day-4_VC","mutant_control_VC"))
write.csv(as.data.frame(res$mutd04VC), file="TOR_mutD0-4VC_res.csv")

# MA Plot
pdf("TOR_MA-plot_mutD0-4VC.pdf", width = 4, height = 4) 
plotMA(res$mutd04VC, alpha = 0.05, ylim = c(-16, 16), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$mutd04VC[which(res$mutd04VC$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$mutd04VC <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_mutD0-4VC_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$mutd04VC <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_mutD0-4VC_up.csv")

# Add entrez gene IDs
res$mutd04VC$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$mutd04VC), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$mutd04VC) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_mutD04VC = as.data.frame(gsea)
write.csv(gsea_table_mutD04VC, file="TOR_GSEA_mutD0-4VC.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_mutD04VC = as.data.frame(Reactome)
write.csv(reactome_table_mutD04VC, file="TOR_Reactome_mutD0-4VC.csv")

# Pathways
path$mutd04VC <- data.frame(pathway = character(),
                            entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$mutd04VC <- rbind(path$mutd04VC, data.frame(pathway = v1, entrezID = v2))
}
path$mutd04VC <- merge(path$mutd04VC, as.data.frame(res$mutd04VC), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$mutd04VC, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_mutD0-4VC.txt")

##### Results D0-D2 wt ######

# Pair-wise comparison
res$wtd02 <- results(dds, contrast = c("Condition","wt_day-2_noVC","wt_control_noVC"))
write.csv(as.data.frame(res$wtd02), file="TOR_wtD0-2_res.csv")

# MA Plot
pdf("TOR_MA-plot_wtD0-2.pdf", width = 4, height = 4) 
plotMA(res$wtd02, alpha = 0.05, ylim = c(-8, 8), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$wtd02[which(res$wtd02$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$wtd02 <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_wtD0-2_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$wtd02 <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_wtD0-2_up.csv")

# Add entrez gene IDs
res$wtd02$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$wtd02), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$wtd02) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_wtD02 = as.data.frame(gsea)
write.csv(gsea_table_wtD02, file="TOR_GSEA_wtD0-2.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_wtD02 = as.data.frame(Reactome)
write.csv(reactome_table_wtD02, file="TOR_Reactome_wtD0-2.csv")

# Pathways
path$wtd02 <- data.frame(pathway = character(),
                         entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$wtd02 <- rbind(path$wtd02, data.frame(pathway = v1, entrezID = v2))
}
path$wtd02 <- merge(path$wtd02, as.data.frame(res$wtd02), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$wtd02, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_wtD0-2.txt")

##### Results D0-D4 wt ######

# Pair-wise comparison
res$wtd04 <- results(dds, contrast = c("Condition","wt_day-4_noVC","wt_control_noVC"))
write.csv(as.data.frame(res$wtd04), file="TOR_wtD0-4_res.csv")

# MA Plot
pdf("TOR_MA-plot_wtD0-4.pdf", width = 4, height = 4) 
plotMA(res$wtd04, alpha = 0.05, ylim = c(-10, 10), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$wtd04[which(res$wtd04$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$wtd04 <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_wtD0-4_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$wtd04 <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_wtD0-4_up.csv")

# Add entrez gene IDs
res$wtd04$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$wtd04), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$wtd04) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_wtD04 = as.data.frame(gsea)
write.csv(gsea_table_wtD04, file="TOR_GSEA_wtD0-4.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_wtD04 = as.data.frame(Reactome)
write.csv(reactome_table_wtD04, file="TOR_Reactome_wtD0-4.csv")

# Pathways
path$wtd04 <- data.frame(pathway = character(),
                         entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$wtd04 <- rbind(path$wtd04, data.frame(pathway = v1, entrezID = v2))
}
path$wtd04 <- merge(path$wtd04, as.data.frame(res$wtd04), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$wtd04, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_wtD0-4.txt")

##### Results D0-D2 wt VC ######

# Pair-wise comparison
res$wtd02VC <- results(dds, contrast = c("Condition","wt_day-2_VC","wt_control_VC"))
write.csv(as.data.frame(res$wtd02VC), file="TOR_wtD0-2VC_res.csv")

# MA Plot
pdf("TOR_MA-plot_wtD0-2VC.pdf", width = 4, height = 4) 
plotMA(res$wtd02VC, alpha = 0.05, ylim = c(-10, 10), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$wtd02VC[which(res$wtd02VC$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$wtd02VC <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_wtD0-2VC_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$wtd02VC <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_wtD0-2VC_up.csv")

# Add entrez gene IDs
res$wtd02VC$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$wtd02VC), keytype = "ENSEMBL", column = "ENTREZID")

list <- as.data.frame(res$wtd02VC) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_wtD02VC = as.data.frame(gsea)
write.csv(gsea_table_wtD02VC, file="TOR_GSEA_wtD0-2VC.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_wtD02VC = as.data.frame(Reactome)
write.csv(reactome_table_wtD02VC, file="TOR_Reactome_wtD0-2VC.csv")

# Pathways
path$wtd02VC <- data.frame(pathway = character(),
                           entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$wtd02VC <- rbind(path$wtd02VC, data.frame(pathway = v1, entrezID = v2))
}
path$wtd02VC <- merge(path$wtd02VC, as.data.frame(res$wtd02VC), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$wtd02VC, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_wtD0-2VC.txt")

##### Results D0-D4 wt VC ######

# Pair-wise comparison
res$wtd04VC <- results(dds, contrast = c("Condition","wt_day-4_VC","wt_control_VC"))
write.csv(as.data.frame(res$wtd04VC), file="TOR_wtD0-4VC_res.csv")

# MA Plot
pdf("TOR_MA-plot_wtD0-4VC.pdf", width = 4, height = 4) 
plotMA(res$wtd04VC, alpha = 0.05, ylim = c(-12, 12), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$wtd04VC[which(res$wtd04VC$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$wtd04VC <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TOR_wtD0-4VC_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$wtd04VC <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TOR_wtD0-4VC_up.csv")

# Add entrez gene IDs
res$wtd04VC$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$wtd04VC), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$wtd04VC) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_wtD04VC = as.data.frame(gsea)
write.csv(gsea_table_wtD04VC, file="TOR_GSEA_wtD0-4VC.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_wtD04VC = as.data.frame(Reactome)
write.csv(reactome_table_wtD04VC, file="TOR_Reactome_wtD0-4VC.csv")

# Pathways
path$wtd04VC <- data.frame(pathway = character(),
                           entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$wtd04VC <- rbind(path$wtd04VC, data.frame(pathway = v1, entrezID = v2))
}
path$wtd04VC <- merge(path$wtd04VC, as.data.frame(res$wtd04VC), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$wtd04VC, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TOR_stats_pathways_wtD0-4VC.txt")

##### TOR group comparisons #####

# Pathways trajectories
path_Z <- merge(path$wtd04VC, counts$Z, by="entrezID")
a = path_Z %>%
  group_by(pathway) %>%
  summarise_at(vars(TOR_1:TOR_36), median) %>%
  mutate(pathway = str_replace_all(pathway, ":", "")) %>%
  mutate(pathway = str_replace_all(pathway, "-", "_")) %>%
  mutate(pathway = str_replace_all(pathway, ",", "")) %>%
  column_to_rownames(var="pathway")
a = as.data.frame(t(a))
a = merge(a, metadata$all, by.x=0, by.y="ID") %>%
  mutate(VC = grepl("_VC", Condition)) %>%
  mutate(across('Info', str_replace, 'control', '0')) %>%
  mutate(Info = as.numeric(Info)) %>%
  separate(col = Condition, into = c("Background","Length","Cocktail"), sep = "_") %>%
  unite("select", Background, Cocktail, sep = "_")
coloring <- c("#FEBC81","#FF7800","#00B6A6","#007369")
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = coloring)
}
for (var in colnames(a)[2:21]) {
  plt = ggplot(a, aes(x=Info, y=.data[[var]], color=select)) +
    geom_point(size=6)+
    geom_smooth(se = F, size=3)+
    scale_x_continuous(breaks=pretty_breaks(n=3))
  ggsave(paste0("TOR_trajectory_", var, ".pdf"), units="mm", dpi = 300)
}

# Heatmap Top DEGs
heat1 = read.csv("TOR_wtD0-2_res.csv") %>%
  mutate(Condition = "wt-D2")
heat2 = read.csv("TOR_wtD0-4_res.csv") %>%
  mutate(Condition = "wt-D4")
heat3 = read.csv("TOR_wtD0-2VC_res.csv") %>%
  mutate(Condition = "wt-D2+VC")
heat4 = read.csv("TOR_wtD0-4VC_res.csv") %>%
  mutate(Condition = "wt-D4+VC")
heat5 = read.csv("TOR_mutD0-2_res.csv") %>%
  mutate(Condition = "mut-D2")
heat6 = read.csv("TOR_mutD0-4_res.csv") %>%
  mutate(Condition = "mut-D4")
heat7 = read.csv("TOR_mutD0-2VC_res.csv") %>%
  mutate(Condition = "mut-D2+VC")
heat8 = read.csv("TOR_mutD0-4VC_res.csv") %>%
  mutate(Condition = "mut-D4+VC")
heat = rbind(heat1, heat2, heat3, heat4,
             heat5, heat6, heat7, heat8) %>%
  mutate(Condition = factor(Condition, levels = c(
    "wt-D2", "wt-D4", "wt-D2+VC", "wt-D4+VC",
    "mut-D2", "mut-D4", "mut-D2+VC", "mut-D4+VC")))
sig_genes_any = unique(heat[heat$padj < 1e-18, "X"])
length(sig_genes_any)
genes <- sig_genes_any
gene_order = heat %>% 
  dplyr::filter(Condition == "mut-D2") %>%
  dplyr::filter(X %in% genes) %>%
  arrange(log2FoldChange) %>%
  pull(X)
heat %>%
  mutate(strain = ifelse(grepl("mut", Condition), "ERCC1", "WT")) %>%
  dplyr::filter(X %in% genes) %>%
  mutate(X = factor(X, levels = gene_order)) %>%
  ggplot(., aes(x=Condition, y=X, fill=log2FoldChange)) +
  geom_tile()+
  scale_fill_gradient2()+
  theme(axis.text.x = element_text(angle=45, hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  facet_grid(~strain, scales = "free")
ggsave("TOR_heatmap.pdf", units="mm", dpi = 300, height = 100)

# Check pathways for statistical significance
path_stats = list()
path_stats$mutd02 = read.table("TOR_stats_pathways_mutD0-2.txt")
path_stats$mutd04 = read.table("TOR_stats_pathways_mutD0-4.txt")
path_stats$mutd02VC = read.table("TOR_stats_pathways_mutD0-2VC.txt")
path_stats$mutd04VC = read.table("TOR_stats_pathways_mutD0-4VC.txt")
path_stats$wtd02 = read.table("TOR_stats_pathways_wtD0-2.txt")
path_stats$wtd04 = read.table("TOR_stats_pathways_wtD0-4.txt")
path_stats$wtd02VC = read.table("TOR_stats_pathways_wtD0-2VC.txt")
path_stats$wtd04VC = read.table("TOR_stats_pathways_wtD0-4VC.txt")
pathways_sig = c()
for (comp in path_stats) {
  pathways_sig = c(pathways_sig, rownames(comp[comp$padj < 0.05, ]))
}
pathways_sig = unique(pathways_sig)

path2 = path
for (p in names(path2)) {
  path2[[p]]$Condition = p
  path2[[p]] = merge(path2[[p]], path_stats[[p]], by.x="pathway", by.y=0)
}
path_cat = do.call(rbind, path2)
go_terms = c(
  "GO:0006281_DNA_repair",
  "GO:0006298_mismatch_repair",
  "GO:1901255_nucleotide-excision_repair_involved_in_interstrand_cross-link_repair",
  "GO:0000724_double-strand_break_repair_via_homologous_recombination",
  "GO:0006303_double-strand_break_repair_via_nonhomologous_end_joining",
  "GO:0097681_double-strand_break_repair_via_alternative_nonhomologous_end_joining",
  "GO:0006289_nucleotide-excision_repair",
  "GO:0006284_base-excision_repair",
  "GO:0007179_transforming_growth_factor_beta_receptor_signaling_pathway",
  "GO:0001837_epithelial_to_mesenchymal_transition",
  "GO:0006325_chromatin_organization"
)
for (gt in go_terms) {
  path_cat %>%
    dplyr::filter(pathway == gt) %>%
    ggplot(., aes(x=Condition, log2FoldChange, fill=pval<0.05)) +
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("lightgray", "red"), name="Statistically significant")+
    lims(y=quantile(path_cat$log2FoldChange, c(0.09, 0.91)))
  short_name = substr(gt, 4, 40)
  ggsave(paste0("TOR_", short_name, ".pdf"), units="mm", dpi = 300, height = 100)
}

# GSEA
gsea_table_mutD02$comparison = "Ercc1 D2"
gsea_table_mutD02VC$comparison = "Ercc1 D2+VC"
gsea_table_mutD04$comparison = "Ercc1 D4"
gsea_table_mutD04VC$comparison = "Ercc1 D4+VC"
gsea_table_wtD02$comparison = "WT D2"
gsea_table_wtD02VC$comparison = "WT D2+VC"
gsea_table_wtD04$comparison = "WT D4"
gsea_table_wtD04VC$comparison = "WT D4+VC"
gsea_all = rbind(
  gsea_table_mutD02,
  gsea_table_mutD02VC,
  gsea_table_mutD04,
  gsea_table_mutD04VC,
  gsea_table_wtD02,
  gsea_table_wtD02VC,
  gsea_table_wtD04,
  gsea_table_wtD04VC
)
rownames(gsea_all) = 1:nrow(gsea_all)
gsea_all %>% 
  group_by(comparison) %>%
  slice_max(abs(NES), n = 8) %>%
  ggplot(., aes(x=comparison, y=Description, fill=NES, size=-log(p.adjust)))+
  geom_point(pch=21)+
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave("TOR_Dotplot_GSEA.pdf", units="mm", dpi = 300, width = 200, height=150)

# Reactome
reactome_table_mutD02$comparison = "Ercc1 D2"
reactome_table_mutD02VC$comparison = "Ercc1 D2+VC"
reactome_table_mutD04$comparison = "Ercc1 D4"
reactome_table_mutD04VC$comparison = "Ercc1 D4+VC"
reactome_table_wtD02$comparison = "WT D2"
reactome_table_wtD02VC$comparison = "WT D2+VC"
reactome_table_wtD04$comparison = "WT D4"
reactome_table_wtD04VC$comparison = "WT D4+VC"
reactome_all = rbind(
  reactome_table_mutD02,
  reactome_table_mutD02VC,
  reactome_table_mutD04,
  reactome_table_mutD04VC,
  reactome_table_wtD02,
  reactome_table_wtD02VC,
  reactome_table_wtD04,
  reactome_table_wtD04VC
)
rownames(reactome_all) = 1:nrow(reactome_all)
reactome_all %>% 
  group_by(comparison) %>%
  slice_max(abs(NES), n = 8) %>%
  ggplot(., aes(x=comparison, y=Description, fill=NES, size=-log(p.adjust)))+
  geom_point(pch=21)+
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave("TOR_Dotplot_Reactome.pdf", units="mm", dpi = 300, width = 240, height=150)


##### TGFB ######

##### Results Repsox ######

# Pair-wise comparison
res$repsox <- results(dds, contrast = c("Condition","mutant_repsox_0","mutant_control_0"))
write.csv(as.data.frame(res$repsox), file="TGFB_Repsox_res.csv")

# MA Plot
pdf("TGFB_MA-plot_Repsox.pdf", width = 4, height = 4) 
plotMA(res$repsox, alpha = 0.05, ylim = c(-8, 8), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$repsox[which(res$repsox$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$repsox <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TGFB_Repsox_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$repsox <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TGFB_Repsox_up.csv")

# Add entrez gene IDs
res$repsox$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$repsox), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$repsox) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_repsox = as.data.frame(gsea)
write.csv(gsea_table_repsox, file="TGFB_GSEA_Repsox.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_repsox = as.data.frame(Reactome)
write.csv(reactome_table_repsox, file="TGFB_Reactome_Repsox.csv")

# Pathways
path <- list()
path$repsox <- data.frame(pathway = character(),
                   entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$repsox <- rbind(path$repsox, data.frame(pathway = v1, entrezID = v2))
}
path$repsox <- merge(path$repsox, as.data.frame(res$repsox), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$repsox, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TGFB_stats_pathways_Repsox.txt")

##### Results Vactoserib ######

# Pair-wise comparison
res$vactoserib <- results(dds, contrast = c("Condition","mutant_vactoserib_0","mutant_control_0"))
write.csv(as.data.frame(res$vactoserib), file="TGFB_Vactoserib_res.csv")

# MA Plot
pdf("TGFB_MA-plot_Vactoserib.pdf", width = 4, height = 4) 
plotMA(res$vactoserib, alpha = 0.05, ylim = c(-8, 8), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$vactoserib[which(res$vactoserib$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$vactoserib <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TGFB_Vactoserib_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$vactoserib <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TGFB_Vactoserib_up.csv")

# Add entrez gene IDs
res$vactoserib$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$vactoserib), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$vactoserib) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_vactoserib = as.data.frame(gsea)
write.csv(gsea_table_vactoserib, file="TGFB_GSEA_Vactoserib.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_vactoserib = as.data.frame(Reactome)
write.csv(reactome_table_vactoserib, file="TGFB_Reactome_Vactoserib.csv")

# Pathways
path$vactoserib <- data.frame(pathway = character(),
                   entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$vactoserib <- rbind(path$vactoserib, data.frame(pathway = v1, entrezID = v2))
}
path$vactoserib <- merge(path$vactoserib, as.data.frame(res$vactoserib), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$vactoserib, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TGFB_stats_pathways_Vactoserib.txt")

##### Results A83-1 ######

# Pair-wise comparison
res$a83 <- results(dds, contrast = c("Condition","mutant_A83-01_0","mutant_control_0"))
write.csv(as.data.frame(res$a83), file="TGFB_A83-01_res.csv")

# MA Plot
pdf("TGFB_MA-plot_A83-01.pdf", width = 4, height = 4) 
plotMA(res$a83, alpha = 0.05, ylim = c(-10, 10), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$a83[which(res$a83$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$a83 <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TGFB_A83-01_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$a83 <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TGFB_A83-01_up.csv")

# Add entrez gene IDs
res$a83$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$a83), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$a83) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_a83 = as.data.frame(gsea)
write.csv(gsea_table_a83, file="TGFB_GSEA_A83-01.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_a83 = as.data.frame(Reactome)
write.csv(reactome_table_a83, file="TGFB_Reactome_A83-01.csv")

# Pathways
path$a83 <- data.frame(pathway = character(),
                   entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$a83 <- rbind(path$a83, data.frame(pathway = v1, entrezID = v2))
}
path$a83 <- merge(path$a83, as.data.frame(res$a83), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$a83, pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TGFB_stats_pathways_A83-01.txt")

##### Results DMH-1 ######

# Pair-wise comparison
res$dmh <- results(dds, contrast = c("Condition","mutant_DMH-1_0","mutant_control_0"))
write.csv(as.data.frame(res$dmh), file="TGFB_DMH-1_res.csv")

# MA Plot
pdf("TGFB_MA-plot_DMH-1.pdf", width = 4, height = 4) 
plotMA(res$dmh, alpha = 0.05, ylim = c(-8, 8), xlab = "Mean of Normalized Counts", ylab = "Fold Change" )
garbage = dev.off()

# List with significant genes
res_sig <- res$dmh[which(res$dmh$padj < 0.05 ), ] 
res_down <- res_sig[res_sig$log2FoldChange < 0,] 
down$dmh <- rownames(res_down)
write.csv(as.data.frame(res_down), file="TGFB_DMH-1_down.csv")
res_up <- res_sig[res_sig$log2FoldChange > 0,] 
up$dmh <- rownames(res_up)
write.csv(as.data.frame(res_up), file="TGFB_DMH-1_up.csv")

# Add entrez gene IDs
res$dmh$entrezID = mapIds(org.Mm.eg.db, keys = row.names(res$dmh), keytype = "ENSEMBL", column = "ENTREZID")
list <- as.data.frame(res$dmh) %>%
  drop_na(entrezID, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, entrezID)

# GSEA GO
gc()
gsea = gseGO(geneList=list, 
             ont ="BP", 
             keyType = "ENTREZID",
             verbose = TRUE, 
             OrgDb = get("org.Mm.eg.db"))
gsea_table_dmh = as.data.frame(gsea)
write.csv(gsea_table_dmh, file="TGFB_GSEA_DMH-1.csv")

# GSEA REACTOME 
Reactome <- gsePathway(list,
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
reactome_table_dmh = as.data.frame(Reactome)
write.csv(reactome_table_dmh, file="TGFB_Reactome_DMH-1.csv")

# Pathways
path$dmh <- data.frame(pathway = character(),
                   entrezID = character())
for (p in pathways) {
  v1 <- rep(p, length(gene_sets[p]))
  v2 <- gene_sets[[p]]
  path$dmh  <- rbind(path$dmh , data.frame(pathway = v1, entrezID = v2))
}
path$dmh  <- merge(path$dmh , as.data.frame(res$dmh), by="entrezID")

# Statistics
path_stats = data.frame()
for (p in pathways) {
  subs = subset(path$dmh , pathway == p)
  test = wilcox.test(subs$log2FoldChange)
  path_stats[p, "n"] = nrow(subs)
  path_stats[p, "median"] = median(subs$log2FoldChange, na.rm = T)
  path_stats[p, "pval"] = round(test$p.value, digits = 6)
}
path_stats$padj = p.adjust(path_stats$pval, method="BH")
write.table(path_stats, file="TGFB_stats_pathways_DMH-1.txt")

##### TGFB group comparisons #####

# GSEA
gsea_table_a83$comparison = "A38-01"
gsea_table_dmh$comparison = "DMH-1"
gsea_table_repsox$comparison = "Repsox"
gsea_table_vactoserib$comparison = "Vactoserib"
gsea_all = rbind(
  gsea_table_a83,
  gsea_table_dmh,
  gsea_table_repsox,
  gsea_table_vactoserib
)
rownames(gsea_all) = 1:nrow(gsea_all)
gsea_all %>% 
  group_by(comparison) %>%
  slice_max(abs(NES), n = 8) %>%
  ggplot(., aes(x=comparison, y=Description, fill=NES, size=-log(p.adjust)))+
  geom_point(pch=21)+
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave("TGFB_Dotplot_GSEA.pdf", units="mm", dpi = 300, width = 220)

# Reactome
reactome_table_a83$comparison = "A38-01"
reactome_table_dmh$comparison = "DMH-1"
reactome_table_repsox$comparison = "Repsox"
reactome_table_vactoserib$comparison = "Vactoserib"
reactome_all = rbind(
  reactome_table_a83,
  reactome_table_dmh,
  reactome_table_repsox,
  reactome_table_vactoserib
)
rownames(reactome_all) = 1:nrow(reactome_all)
reactome_all %>% 
  group_by(comparison) %>%
  slice_max(abs(NES), n = 8) %>%
  ggplot(., aes(x=comparison, y=Description, fill=NES, size=-log(p.adjust)))+
  geom_point(pch=21)+
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave("TGFB_Dotplot_Reactome.pdf", units="mm", dpi = 300, width = 220)

# Venn diagram downregulated genes
venn <-list('Repsox'= down$repsox,'Vactoserib'= down$vactoserib,'A83-01'= down$a83,'DMH-1'= down$dmh)
pdf("TGFB-Venn-down.pdf", width = 10, height = 4) 
ggvenn(venn, fill_color = c("blue", "lightblue", "pink", "red"),set_name_size = 4, text_size = 2.5)
garbage = dev.off()
TGFB_down_all <- Reduce(intersect,list(down$repsox,down$vactoserib,down$a83,down$dmh))
write.csv(TGFB_down_all, file="TGFB_down_all.csv")

# Venn diagram upregulated genes
venn <-list('Repsox'= up$repsox,'Vactoserib'= up$vactoserib,'A83-01'= up$a83,'DMH-1'= up$dmh)
pdf("TGFB-Venn-up.pdf", width = 10, height = 4) 
ggvenn(venn, fill_color = c("blue", "lightblue", "pink", "red"),set_name_size = 4, text_size = 2.5)
garbage = dev.off()
TGFB_up_all <- Reduce(intersect,list(up$repsox,up$vactoserib,up$a83,up$dmh))
write.csv(TGFB_up_all, file="TGFB_up_all.csv")

# Create background list and add entrez IDs
counts$filt$entrezID = mapIds(org.Mm.eg.db, keys = row.names(counts$filt), keytype = "ENSEMBL", column = "ENTREZID")
list_background = as.data.frame(counts$filt) %>%
  drop_na(entrezID)
list_background <- list_background$entrezID

# Add entrez IDs
TGFB_up_all$entrezID = mapIds(org.Mm.eg.db, keys = TGFB_up_all, keytype = "ENSEMBL", column = "ENTREZID")
TGFB_up_all_entrez = as.data.frame(TGFB_up_all) %>%
  drop_na(entrezID)
TGFB_up_all_entrez <- TGFB_up_all_entrez$entrezID

# GO analysis
ego <- enrichGO(gene          = TGFB_up_all_entrez,
                universe      = list_background,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)
GO_table_TGFB_up_all = as.data.frame(ego)
write.csv(GO_table_TGFB_up_all, file="TGFB_GO_up_all.csv")

# Add entrez IDs
TGFB_down_all$entrezID = mapIds(org.Mm.eg.db, keys = TGFB_down_all, keytype = "ENSEMBL", column = "ENTREZID")
TGFB_down_all_entrez = as.data.frame(TGFB_down_all) %>%
  drop_na(entrezID)
TGFB_down_all_entrez <- TGFB_down_all_entrez$entrezID

ego <- enrichGO(gene          = TGFB_down_all_entrez,
                universe      = list_background,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)
GO_table_TGFB_down_all = as.data.frame(ego)
write.csv(GO_table_TGFB_down_all, file="TGFB_GO_down_all.csv")

# Barplot
tail1 = GO_table_TGFB_up_all %>%
  mutate(dir = "up") %>%
  slice_min(p.adjust, n=8)
tail2 = GO_table_TGFB_down_all %>%
  mutate(dir = "down") %>%
  slice_min(p.adjust, n=8)
go_tails = rbind(tail1, tail2) %>%
  mutate(logp.adjust = -log10(p.adjust))
go_tails %>%
  mutate(Description = fct_reorder(Description, p.adjust)) %>%
  ggplot(., aes(x=Description, y=logp.adjust, fill=dir))+
  geom_bar(stat="identity")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill="none")+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  scale_fill_manual(values=c("blue","red"))
ggsave("TGFB_GO_comparison.pdf", units="mm", dpi = 300, height = 130)

# Pathways
path_stats = list()
path_stats$repsox = read.table("TGFB_stats_pathways_repsox.txt")
path_stats$vactoserib = read.table("TGFB_stats_pathways_vactoserib.txt")
path_stats$a83 = read.table("TGFB_stats_pathways_A83-01.txt")
path_stats$dmh = read.table("TGFB_stats_pathways_DMH-1.txt")
pathways_sig = c()
for (comp in path_stats) {
  pathways_sig = c(pathways_sig, rownames(comp[comp$padj < 0.05, ]))
}
pathways_sig = unique(pathways_sig)

path2 = path
for (p in names(path2)) {
  path2[[p]]$Condition = p
  path2[[p]] = merge(path2[[p]], path_stats[[p]], by.x="pathway", by.y=0)
}
path_cat = do.call(rbind, path2)

go_terms = c(
  "GO:0006281_DNA_repair",
  "GO:0006298_mismatch_repair",
  "GO:1901255_nucleotide-excision_repair_involved_in_interstrand_cross-link_repair",
  "GO:0000724_double-strand_break_repair_via_homologous_recombination",
  "GO:0006303_double-strand_break_repair_via_nonhomologous_end_joining",
  "GO:0097681_double-strand_break_repair_via_alternative_nonhomologous_end_joining",
  "GO:0006289_nucleotide-excision_repair",
  "GO:0006284_base-excision_repair",
  "GO:0007179_transforming_growth_factor_beta_receptor_signaling_pathway",
  "GO:0001837_epithelial_to_mesenchymal_transition",
  "GO:0006325_chromatin_organization",
  "GO:0061621_canonical_glycolysis",
  "GO:0000278_mitotic_cell_cycle"
)
for (gt in go_terms) {
  path_cat %>%
    dplyr::filter(pathway == gt) %>%
    ggplot(., aes(x=Condition, log2FoldChange, fill=pval<0.05)) +
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("lightgray", "red"), name="Statistically significant")+
    lims(y=quantile(path_cat$log2FoldChange, c(0.09, 0.91)))
  short_name = substr(gt, 4, 40)
  ggsave(paste0("TGFB_", short_name, ".pdf"), units="mm", dpi = 300, height = 100)
}

##### Comparison TOR and TGFB #####

# Venn Diagram downregulated
venn <-list('Repsox' = down$repsox, 'Vactoserib'= down$vactoserib, 'A83-01' = down$a83, 'D4' = down$mutd04)
pdf("TGFB_TOR_Venn_down_D4.pdf", width = 10, height = 4) 
ggvenn(venn, fill_color = c("lightblue", "blue","darkblue","red"),set_name_size = 4, text_size = 2.5)
garbage = dev.off()
TGFB_TOR_down_D4 <- Reduce(intersect,list(down$repsox,down$vactoserib,down$a83,down$mutd04))
write.csv(TGFB_TOR_down_D4, file="TGFB_TOR_down_D4.csv")

# Add entrez IDs
TGFB_TOR_down_D4$entrezID = mapIds(org.Mm.eg.db, keys = TGFB_TOR_down_D4, keytype = "ENSEMBL", column = "ENTREZID")
TGFB_TOR_down_entrez = as.data.frame(TGFB_TOR_down_D4) %>%
  drop_na(entrezID)
TGFB_TOR_down_entrez <- TGFB_TOR_down_entrez$entrezID

# GO analysis
ego <- enrichGO(gene          = TGFB_TOR_down_entrez,
                universe      = list_background,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)
GO_table_TGFB_TOR_down_D4 = as.data.frame(ego)
write.csv(GO_table_TGFB_TOR_down_D4, file="TGFB_TOR_GO_down_D4.csv")

# Reactome analysis
react <- enrichPathway(gene          = TGFB_TOR_down_entrez,
                       universe      = list_background,
                       organism = "mouse",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
Reactome_table_TGFB_TOR_down_D4 = as.data.frame(react)
write.csv(Reactome_table_TGFB_TOR_down_D4, file="TGFB_TOR_Reactome_down_D4.csv")

# Venn Diagram upregulated
venn <-list('Repsox' = up$repsox,'Vactoserib'= up$vactoserib, 'A83' = up$a83, 'D4' = up$mutd04)
pdf("TGFB_TOR_Venn_up_D4.pdf", width = 10, height = 4) 
ggvenn(venn, fill_color = c("lightblue", "blue","darkblue","red"),set_name_size = 4, text_size = 2.5)
garbage = dev.off()
TGFB_TOR_up_D4 <- Reduce(intersect,list(up$repsox,up$vactoserib,up$a83,up$mutd04))
write.csv(TGFB_TOR_up_D4, file="TGFB_TOR_up_D4.csv")

# Add entrez IDs
TGFB_TOR_up_D4$entrezID = mapIds(org.Mm.eg.db, keys = TGFB_TOR_up_D4, keytype = "ENSEMBL", column = "ENTREZID")
TGFB_TOR_up_entrez = as.data.frame(TGFB_TOR_up_D4) %>%
  drop_na(entrezID)
TGFB_TOR_up_entrez <- TGFB_TOR_up_entrez$entrezID

# GO analysis
ego <- enrichGO(gene          = TGFB_TOR_up_entrez,
                universe      = list_background,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)
GO_table_TGFB_TOR_up_D4 = as.data.frame(ego)
write.csv(GO_table_TGFB_TOR_up_D4, file="TGFB_TOR_GO_up.csv")

# Reactome analysis
react <- enrichPathway(gene          = TGFB_TOR_up_entrez,
                       universe      = list_background,
                       organism = "mouse",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
Reactome_table_TGFB_TOR_up_D4 = as.data.frame(react)
write.csv(Reactome_table_TGFB_TOR_up_D4, file="TGFB_TOR_Reactome_up_D4.csv")

# GO
tail1 = GO_table_TGFB_TOR_up_D4 %>%
  mutate(dir = "up") %>%
  slice_min(p.adjust, n=6)
tail2 = GO_table_TGFB_TOR_down_D4 %>%
  mutate(dir = "down") %>%
  slice_min(p.adjust, n=6)
go_tails = rbind(tail1, tail2) %>%
  mutate(logp.adjust = -log10(p.adjust))
go_tails %>%
  mutate(Description = fct_reorder(Description, p.adjust)) %>%
  ggplot(., aes(x=Description, y=logp.adjust, fill=dir))+
  geom_bar(stat="identity")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill="none")+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  scale_fill_manual(values=c("blue","red"))
ggsave("TGFB_TOR_GO_comparison.pdf", units="mm", dpi = 300, height = 130)

# Reactome
tail1 = Reactome_table_TGFB_TOR_up_D4 %>%
  mutate(dir = "up") %>%
  slice_min(p.adjust, n=6)
tail2 = Reactome_table_TGFB_TOR_down_D4 %>%
  mutate(dir = "down") %>%
  slice_min(p.adjust, n=6)
go_tails = rbind(tail1, tail2) %>%
  mutate(logp.adjust = -log10(p.adjust))
go_tails %>%
  mutate(Description = fct_reorder(Description, p.adjust)) %>%
  ggplot(., aes(x=Description, y=logp.adjust, fill=dir))+
  geom_bar(stat="identity")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill="none")+
  theme(axis.title.y = element_blank())+
  coord_flip()+
  scale_fill_manual(values=c("blue","red"))
ggsave("TGFB_TOR_Reactome_comparison.pdf", units="mm", dpi = 300, height = 130)

# Heatmap Top DEGs
heat1 = read.csv("TGFB_Repsox_res.csv") %>%
  mutate(Condition = "Repsox")
heat2 = read.csv("TGFB_Vactoserib_res.csv") %>%
  mutate(Condition = "Vacto")
heat3 = read.csv("TGFB_A83-01_res.csv") %>%
  mutate(Condition = "A83")
heat4 = read.csv("TGFB_DMH-1_res.csv") %>%
  mutate(Condition = "DMH")
heat5 = read.csv("TOR_mutD0-2_res.csv") %>%
  mutate(Condition = "D2")
heat6 = read.csv("TOR_mutD0-4_res.csv") %>%
  mutate(Condition = "D4")
heat7 = read.csv("TOR_mutD0-2VC_res.csv") %>%
  mutate(Condition = "D2+VC")
heat8 = read.csv("TOR_mutD0-4VC_res.csv") %>%
  mutate(Condition = "D4+VC")
heat = rbind(heat1, heat2, heat3, heat4,
             heat5, heat6, heat7, heat8)
sig_genes_any = unique(heat[heat$padj < 1e-10, "X"])
length(sig_genes_any)
genes <- sig_genes_any
gene_order = heat %>% 
  dplyr::filter(Condition == "D2") %>%
  dplyr::filter(X %in% genes) %>%
  arrange(log2FoldChange) %>%
  pull(X)
heat %>%
  mutate(Condition = factor(Condition, levels = c(
    "D2", "D4", "D2+VC", "D4+VC",
    "Repsox", "Vacto", "A83", "DMH"))) %>%
  dplyr::filter(X %in% genes) %>%
  mutate(X = factor(X, levels = gene_order)) %>%
  ggplot(., aes(x=Condition, y=X, fill=log2FoldChange)) +
  geom_tile()+
  scale_fill_gradient2()+
  theme(axis.text.x = element_text(angle=45, hjust=1),axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave("TGFB_TOR_heatmap.pdf", units="mm", dpi = 300, height = 100)

# Correlation matrix
tmp = list()
comps = c("mutd02", "mutd04", "mutd02VC", "mutd04VC", "repsox", "vactoserib", "a83", "dmh")
for (comp in comps) {
  tmp[[comp]] = as.data.frame(res[[comp]])$log2FoldChange
}
tmp = do.call(cbind, tmp)
trt_cor = data.frame(row.names = comps)
trt_pval = data.frame(row.names = comps)
for (c1 in comps) {
  for (c2 in comps) {
    test = cor.test(tmp[, c1], tmp[, c2])
    trt_cor[c1, c2] = test$estimate
    trt_pval[c1, c2] = test$p.value
  }
}
trt_cor = trt_cor%>%
  mutate(Trt1 = rownames(.)) %>%
  pivot_longer(cols=-Trt1, names_to = "Trt2", values_to = "Cor") %>%
  mutate(Trt1 = factor(Trt1, levels = comps))%>%
  mutate(Trt2 = factor(Trt2, levels = rev(comps)))
ggplot(trt_cor, aes(x=Trt1, y=Trt2, fill=Cor))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  geom_text(aes(label=round(Cor, 2)))
ggsave("TGFB_TOR_Correlation-matrix.pdf", units="mm", dpi = 300, width = 220, height=600)


