
BiocManager::install("sva")
library(TCGAbiolinks)
library(GEOquery)
library(DESeq2)
library(limma)
library(WGCNA)
library(clusterProfiler)
library(survival)
library(survminer)
library(org.Hs.eg.db)
library(sva)
library(ggplot2)
library(dplyr)
library(ggplot2)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  sample.type = "Primary Tumor"
)

setwd("/Users/thejaarlagadda/Library/CloudStorage/OneDrive-IndianaUniversity/Semester_3/System_Biology/Project/new")

####---- Step1- Download GEO data ----####
gse <- getGEO("GSE76250", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])  # Expression matrix
pdata <- pData(gse[[1]]) # Clinical/phenotype data

expr <- expr[, rownames(pdata)]  # Ensure columns align with sample info

# Extract sample condition (TNBC vs Normal)
pdata$condition <- ifelse(grepl("TNBC", pdata$`tissue type:ch1`, ignore.case = TRUE), "TNBC", "Normal")

# Sanity check
table(pdata$condition)  # Should be 165 TNBC and 33 Normal

library(limma)

# Design matrix
group <- factor(pdata$condition)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit model
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(TNBC-Normal, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
deg_results <- topTable(fit2, number=Inf, adjust="fdr")
deg_sig <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)

head(deg_sig)

# Download platform annotation for GPL17586
gpl <- getGEO("GPL17586", AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# Preview available columns
head(colnames(gpl_table))

# Find the mapping columns
# Typically: "ID" → TC ID, "gene_assignment" or "Gene Symbol" → gene names
# Extract gene symbols from the gene_assignment column
extract_gene_symbol <- function(x) {
  # Split on ' /// ' to separate multiple entries
  parts <- strsplit(x, " /// ")[[1]]
  # Extract the second field from the first entry
  first_part <- strsplit(parts[1], " // ")[[1]]
  if (length(first_part) >= 2) return(first_part[2]) else return(NA)
}

# Apply this function to the gene_assignment column
gpl_table$Gene_Symbol <- sapply(gpl_table$gene_assignment, extract_gene_symbol)

# Create mapping table
id_map <- gpl_table[, c("ID", "Gene_Symbol")]
colnames(id_map) <- c("TC_ID", "Gene_Symbol")

# Add TC_ID to your deg_sig object
deg_sig$TC_ID <- rownames(deg_sig)

# Merge DEG results with gene symbols
deg_sig_mapped <- merge(deg_sig, id_map, by = "TC_ID")




nrow(deg_immune)
head(deg_immune[, c("Gene_Symbol", "logFC", "adj.P.Val")])
table(is.na(deg_sig_mapped$Gene_Symbol))
head(deg_sig_mapped$Gene_Symbol, 20)
head(immune_genes, 20)


# Loading immune geneslist
immune_genes <- read.table("/Users/thejaarlagadda/Library/CloudStorage/OneDrive-IndianaUniversity/Semester_3/System_Biology/Project/immune_genes.txt", 
                           header = FALSE, stringsAsFactors = FALSE)[,1]
# Filter for immune-related genes
deg_immune <- deg_sig_mapped[deg_sig_mapped$Gene_Symbol %in% immune_genes, ]
nrow(deg_immune)
head(deg_immune)


# Adjust thresholds if needed
logFC_cutoff <- 1
pval_cutoff <- 0.05

# Upregulated genes
deg_up <- deg_sig[deg_sig$logFC > logFC_cutoff & deg_sig$adj.P.Val < pval_cutoff, ]
up_genes <- rownames(deg_up)

# Downregulated genes
deg_down <- deg_sig[deg_sig$logFC < -logFC_cutoff & deg_sig$adj.P.Val < pval_cutoff, ]
down_genes <- rownames(deg_down)

# basic barplot of DEG counts
png("DEG_barplot.png", width = 800, height = 600)
barplot(c(nrow(deg_up), nrow(deg_down)),
        names.arg = c("Upregulated", "Downregulated"),
        col = c("red", "blue"),
        main = "Differentially Expressed Genes",
        ylab = "Number of Genes")
dev.off()

library(ggplot2)
library(ggrepel)  # For smart text labels

# Define thresholds
logfc_cutoff <- 1
padj_cutoff <- 0.05
top_n <- 10  # Number of top genes to annotate

# Add significance label
deg_sig$Significance <- "Not Significant"
deg_sig$Significance[deg_sig$logFC > logfc_cutoff & deg_sig$adj.P.Val < padj_cutoff] <- "Upregulated"
deg_sig$Significance[deg_sig$logFC < -logfc_cutoff & deg_sig$adj.P.Val < padj_cutoff] <- "Downregulated"

# Add gene names as a column
deg_sig$Gene <- rownames(deg_sig)

# Select top N genes for annotation
top_up <- deg_sig[deg_sig$Significance == "Upregulated", ]
top_down <- deg_sig[deg_sig$Significance == "Downregulated", ]

top_genes <- rbind(
  head(top_up[order(-top_up$logFC), ], top_n),
  head(top_down[order(top_down$logFC), ], top_n)
)
library(ggplot2)
# Volcano plot with annotations
ggplot(deg_sig, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  geom_text_repel(data = top_genes,
                  aes(label = Gene),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color = 'black') +
  theme_minimal() +
  labs(title = "Volcano Plot with Top Annotated Genes",
       x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")

dev.off()
####---- To extract immune genelist ----####
# Install msigdbr if not already installed
install.packages("msigdbr")

library(msigdbr)

# Get immune-related gene sets from humans
immune_msig <- msigdbr(species = "Homo sapiens", category = "C7")

# Optional: filter specific immune processes like T-cells, cytokines, etc.
immune_genes <- unique(immune_msig$gene_symbol)
write.table(immune_genes, file = "immune_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


####----Step 3 - 
library(ggplot2)

# Volcano Plot with Immune Gene Highlight
deg_sig_mapped$Immune <- ifelse(deg_sig_mapped$Gene_Symbol %in% immune_genes, "Immune", "Other")

ggplot(deg_sig_mapped, aes(x = logFC, y = -log10(adj.P.Val), color = Immune)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Immune" = "red", "Other" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Immune vs Non-Immune DEGs",
       x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")

# Export Top Immune DEGs
write.csv(deg_immune, "immune_related_DEGs.csv", row.names = FALSE)


####----step 4 Pathway Enrichment (GO/KEGG) for Immune DEGs

# Required packages
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(deg_immune$Gene_Symbol, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# GO: Biological Process
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# KEGG Pathway
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                    organism = "hsa",
                    pvalueCutoff = 0.05)

# Visualize
dotplot(ego, showCategory = 10, title = "Top GO Biological Processes")
dotplot(ekegg, showCategory = 10, title = "Top KEGG Pathways")

####----Protein-Protein Interaction Network via STRING----####

library(STRINGdb)

string_db <- STRINGdb$new(version="11.5", species=9606,
                          score_threshold=400, input_directory="")

# Map symbols to STRING IDs
mapped <- string_db$map(data.frame(gene=deg_immune$Gene_Symbol),
                        "gene", removeUnmappedRows = TRUE)

# Plot STRING network
string_db$plot_network(mapped$STRING_id)

png("STRING_network.png", width = 1000, height = 800)
string_db$plot_network(mapped$STRING_id)
dev.off()

####---- Gene Co-Expression Network (WGCNA) ----####
# Load WGCNA
BiocManager::install("WGCNA")
library(WGCNA)

# Allow WGCNA to use all CPU cores
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Transpose expression data (WGCNA expects genes as columns)
datExpr <- as.data.frame(t(expr))  # Dimensions: samples x genes

# Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Match traits
trait_data <- pdata[, c("ki67(%):ch1", "age (yrs):ch1", "tumor grade:ch1")]
rownames(trait_data) <- rownames(pdata)
trait_data <- trait_data[rownames(datExpr), ]
trait_data <- as.data.frame(lapply(trait_data, function(x) as.numeric(as.character(x))))

# Pick soft-threshold power
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot soft thresholding results
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",
     type="n", main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")


#Filter Genes Before WGCNA Only use the most variable genes, e.g., top 5000–10000
# Calculate variance for each gene
gene_variances <- apply(datExpr, 2, var)

# Select top 5000 most variable genes
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:5000]

# Subset datExpr
datExpr_top <- datExpr[, top_genes]


# Choose power (e.g., power=6) based on above plot
softPower <- 6
adjacency <- adjacency(datExpr, power = softPower)

# Replace previous adjacency line
adjacency <- adjacency(datExpr_top, power = softPower)


# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Cluster genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Identify modules using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)

# Convert numeric labels into colors
moduleColors <- labels2colors(dynamicMods)

# Plot dendrogram
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03)

# Calculate module eigengenes
MEList <- moduleEigengenes(datExpr_top, colors = moduleColors)
MEs <- MEList$eigengenes

# Correlate modules with traits
moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr_top))

# Heatmap of module-trait relationships
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_data),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleTraitCor, 2),
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-Trait Relationships")

# Extract genes from turquoise module
turquoise_genes <- names(datExpr_top)[which(moduleColors == "turquoise")]

# Find overlap with immune-related DEGs
hub_immune_genes <- intersect(turquoise_genes, deg_immune$Gene_Symbol)

# Save hub gene list
write.table(hub_immune_genes, "hub_immune_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)



# Use the mapping table you already made
# deg_sig_mapped has: TC_ID (e.g., TC01000563.hg.1) and Gene_Symbol

# Get gene symbols from turquoise TC IDs
turquoise_gene_symbols <- deg_sig_mapped$Gene_Symbol[deg_sig_mapped$TC_ID %in% turquoise_genes]

# Now match to immune DEGs
hub_immune_genes <- intersect(toupper(turquoise_gene_symbols), immune_genes)

# Print them
hub_immune_genes

colnames(pdata)
head(pdata)
