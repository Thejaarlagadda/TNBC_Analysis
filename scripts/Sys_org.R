# TNBC Immune-Related Gene Analysis
# Author: Theja
# Description: Complete pipeline for DEG, Immune Gene Filtering, WGCNA, Enrichment, and PPI

# =========================
# 1. LOAD LIBRARIES
# =========================
library(GEOquery)
library(DESeq2)
library(limma)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(msigdbr)
library(STRINGdb)
library(ggplot2)
library(ggrepel)
library(dplyr)

# =========================
# 2. LOAD & PROCESS GEO DATA
# =========================
gse <- getGEO("GSE76250", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
pdata <- pData(gse[[1]])
pdata$condition <- ifelse(grepl("TNBC", pdata$`tissue type:ch1`, ignore.case = TRUE), "TNBC", "Normal")

# Align expr with pdata
expr <- expr[, rownames(pdata)]

# =========================
# 3. DIFFERENTIAL EXPRESSION (limma)
# =========================
group <- factor(pdata$condition)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, makeContrasts(TNBC-Normal, levels=design))
fit2 <- eBayes(fit2)
deg_results <- topTable(fit2, number=Inf, adjust="fdr")
deg_sig <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)

# =========================
# 4. ANNOTATE GENE SYMBOLS
# =========================
gpl <- getGEO("GPL17586", AnnotGPL = TRUE)
gpl_table <- Table(gpl)
extract_gene_symbol <- function(x) {
  parts <- strsplit(x, " /// ")[[1]]
  first_part <- strsplit(parts[1], " // ")[[1]]
  if (length(first_part) >= 2) return(first_part[2]) else return(NA)
}
gpl_table$Gene_Symbol <- sapply(gpl_table$gene_assignment, extract_gene_symbol)
id_map <- gpl_table[, c("ID", "Gene_Symbol")]
colnames(id_map) <- c("TC_ID", "Gene_Symbol")
deg_sig$TC_ID <- rownames(deg_sig)
deg_sig_mapped <- merge(deg_sig, id_map, by = "TC_ID")

# =========================
# 5. FILTER IMMUNE GENES
# =========================
immune_genes <- read.table("data/immune_genes.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
deg_immune <- deg_sig_mapped[deg_sig_mapped$Gene_Symbol %in% immune_genes, ]
write.csv(deg_immune, "results/immune_related_DEGs.csv", row.names = FALSE)

# =========================
# 6. PLOT DEGS AND ANNOTATE TOP GENES
# =========================
deg_sig$Significance <- "Not Significant"
deg_sig$Significance[deg_sig$logFC > 1 & deg_sig$adj.P.Val < 0.05] <- "Upregulated"
deg_sig$Significance[deg_sig$logFC < -1 & deg_sig$adj.P.Val < 0.05] <- "Downregulated"
deg_sig$Gene <- rownames(deg_sig)
top_genes <- rbind(
  head(deg_sig[deg_sig$Significance == "Upregulated", ][order(-deg_sig$logFC), ], 10),
  head(deg_sig[deg_sig$Significance == "Downregulated", ][order(deg_sig$logFC), ], 10)
)

# Volcano plot
ggplot(deg_sig, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot with Top Annotated Genes", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")

# =========================
# 7. GO/KEGG ENRICHMENT FOR IMMUNE GENES
# =========================
entrez_ids <- bitr(deg_immune$Gene_Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

# =========================
# 8. STRING PPI NETWORK
# =========================
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")
mapped <- string_db$map(data.frame(gene=deg_immune$Gene_Symbol), "gene", removeUnmappedRows = TRUE)
string_db$plot_network(mapped$STRING_id)

png("STRING_network.png", width = 1000, height = 800)
string_db$plot_network(mapped$STRING_id)
dev.off()

# =========================
# 9. WGCNA MODULE DETECTION
# =========================
datExpr <- as.data.frame(t(expr))
gene_variances <- apply(datExpr, 2, var)
datExpr_top <- datExpr[, names(sort(gene_variances, decreasing = TRUE))[1:5000]]
trait_data <- pdata[, c("ki67(%):ch1", "age (yrs):ch1", "tumor grade:ch1")]
trait_data <- as.data.frame(lapply(trait_data, function(x) as.numeric(as.character(x))))
trait_data <- trait_data[rownames(datExpr_top), ]
softPower <- 6
adjacency <- adjacency(datExpr_top, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
moduleColors <- labels2colors(dynamicMods)
MEList <- moduleEigengenes(datExpr_top, colors = moduleColors)
MEs <- MEList$eigengenes
moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr_top))

# =========================
# 10. EXTRACT HUB IMMUNE GENES
# =========================
turquoise_genes <- names(datExpr_top)[which(moduleColors == "turquoise")]
turquoise_gene_symbols <- deg_sig_mapped$Gene_Symbol[deg_sig_mapped$TC_ID %in% turquoise_genes]
hub_immune_genes <- intersect(toupper(turquoise_gene_symbols), immune_genes)
write.table(hub_immune_genes, "results/hub_immune_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
