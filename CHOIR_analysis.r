suppressPackageStartupMessages({
  library(CHOIR)
  library(countsplit)
  library(ggnewscale)
  library(ragg)
  library(scRNAseq)
  library(Seurat)
  library(tictoc)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(pheatmap)
  library(SingleCellExperiment)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(dplyr)
  library(plyr)
})

#Pull in the Chen et al., 2017 hypothalamus dataset as a SingleCellExperiment object
sce <- ChenBrainData(ensembl = TRUE, location = TRUE)

sce

# Cluster ID assigned by the original study
unique(sce$SVM_clusterID)

gene_annot <- rowData(sce)

# Check what's there
head(gene_annot)

gene_info <- data.frame(
  gene = rownames(sce),
  gene_symbol = gene_annot$originalName,
  ensembl_id = rownames(sce),
  stringsAsFactors = FALSE
)

head(gene_info)

# Extract the raw counts matrix (unnormalized and unfiltered)
# Rows: genes
# Columns: cells
counts_matrix <- as.matrix(assay(sce, "counts"))

# Assign column names explicitly from sce colnames if missing
if (is.null(colnames(counts_matrix))) {
  colnames(counts_matrix) <- colnames(sce)
}

# Create Seurat object: filter genes and cells
object <- CreateSeuratObject(
  counts = counts_matrix,
  min.features = 100,   # min genes per cell
  min.cells = 5         # min cells per gene
)

length(colnames(object)) #I think count is less because of the filtering step above
length(colnames(sce))

# Make sure the size of sce is same as the seurat object
sce_filtered <- sce[, colnames(object)]
sce <- sce_filtered

# Verify
all(colnames(sce) == colnames(object))  # Should be TRUE

# Add gene annotations from rowData(sce) to the Seurat object
# Store gene annotations separately inside Seurat object's misc slot
object@misc$gene_annotations <- as.data.frame(rowData(sce))

# Normalize the data (log-normalization by default)
object <- NormalizeData(object)

# Seurat object now has gene annotations stored in: object@misc$gene_annotations
object

#Parallelize analysis
options(future.globals.maxSize = 2.0 * 1e9)
n_cores = 8

# COMMMENTING BECAUSE ALREADY CALCULATED

# #Use CHOIR func with default parameters
# object <- CHOIR(object, n_cores = n_cores)

# #Save CHOIR object
# saveRDS(object, file = "results/choir_object.rds")

#Read CHOIR object
object <- readRDS("results/choir_object.rds")

# Add original cluster labels to Seurat meta.data
object$Original_Cluster <- sce$SVM_clusterID

head(object@meta.data)

names(object@misc$CHOIR$var_features)

object <- runCHOIRumap(object, reduction = "P0_reduction")
names(object@misc$CHOIR$reduction)

options(repr.plot.width = 8, repr.plot.height = 6)

# Display inline plot
plotCHOIR(object, key = "CHOIR", reduction = "P0_reduction_UMAP")

#Save to file
pdf(file = "results/CHOIR_umap_plot.pdf", width = 8, height = 6)
plotCHOIR(object, key = "CHOIR", reduction = "P0_reduction_UMAP")
dev.off()

# Get top markers per CHOIR cluster using Seurat
# Set the cluster identities as the active identity of each cell
Idents(object) <- "CHOIR_clusters_0.05"

# Double check 
table(Idents(object))

# COMMMENTING BECAUSE ALREADY CALCULATED

# all_markers <- FindAllMarkers(
#   object,             # your Seurat object
#   only.pos = TRUE,    # only markers more expressed in cluster
#   min.pct = 0.25,     # gene expressed in at least 25% of cells in cluster
#   logfc.threshold = 0.25  # minimum log fold change of 0.25
# )

# # Now annotate the marker genes before saving
# all_markers <- all_markers %>%
#   left_join(gene_info, by = "gene")

# head(all_markers)

# # Save the annotated markers instead of raw
# saveRDS(all_markers, file = "results/all_markers.rds")

#Read all_markers object
all_markers <- readRDS("results/all_markers.rds")

all_markers_annotated <- all_markers

# Extract the gene symbols and map from your original sce rowData
gene_map <- data.frame(
  feature = rownames(sce),            # Ensembl IDs like ENSMUSG...
  gene_symbol = rowData(sce)$originalName,
  stringsAsFactors = FALSE
)

# Mapping from FeatureXXXX to Ensembl IDs.
# Extract numeric part of the Feature IDs in all_markers
all_markers_annotated$feature_num <- as.integer(gsub("Feature", "", all_markers_annotated$gene))

# Map gene symbols from sce rowData by row index
all_markers_annotated$gene_symbol <- gene_map$gene_symbol[all_markers_annotated$feature_num]

head(all_markers_annotated)

# Save the annotated markers instead of raw
saveRDS(all_markers_annotated, file = "results/all_markers_annotated.rds")
write.csv(all_markers_annotated, "results/all_markers_annotated.csv", row.names = FALSE)

#Read all_markers object
all_markers_annotated <- readRDS("results/all_markers_annotated.rds")

# Original Clusters from the study
table(object$Original_Cluster)

# CHOIR clusters
table(object$CHOIR_clusters_0.05)

length(unique(object$Original_Cluster))
length(unique(object$CHOIR_clusters_0.05))

# Quick Overview - how the new clusters map onto the old ones
cross_tab <- table(object$CHOIR_clusters_0.05, object$Original_Cluster)
print(cross_tab)

# Convert the contingency table (cross-tabulation) to a data frame
cross_tab_df <- as.data.frame(cross_tab)

# Rename the columns of the data frame for clarity
# CHOIR_cluster: the cluster label assigned by CHOIR
# Original_cluster: the original cluster label (e.g., from Seurat or manual annotation)
# count: number of cells that belong to both the CHOIR and Original cluster
names(cross_tab_df) <- c("CHOIR_cluster", "Original_cluster", "count")

# For each CHOIR cluster, find the Original cluster that has the most overlapping cells
best_matches <- cross_tab_df %>%
  group_by(CHOIR_cluster) %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup()

best_matches

# Extract the UMAP coords matrix from misc slot
umap_coords <- object@misc$CHOIR$reduction$P0_reduction_UMAP

# Add it as a DimReduc object to your Seurat object
object[["CHOIR_UMAP"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "CHOIRUMAP_",
  assay = DefaultAssay(object)
)

options(repr.plot.width = 8, repr.plot.height = 6)

# UMAP colored by original cluster labels
DimPlot(object, reduction = "CHOIR_UMAP", group.by = "Original_Cluster", label = TRUE, repel = TRUE) + 
  ggtitle("Original Clusters")

# UMAP colored by CHOIR clusters
DimPlot(object, reduction = "CHOIR_UMAP", group.by = "CHOIR_clusters_0.05", label = TRUE, repel = TRUE) + 
  ggtitle("CHOIR Clusters")

# Get markers for CHOIR clusters 4 and 7
# tells which genes are most upregulated in each of those two CHOIR subclusters that contain OPCs
opc_markers <- all_markers_annotated %>% 
  dplyr::filter(cluster %in% c(4, 7)) %>%    # cluster represents CHOIR cluster
  arrange(cluster, desc(avg_log2FC))

opc_markers

# Extract top genes per cluster
top_4 <- opc_markers %>% dplyr::filter(cluster == 4 & p_val_adj < 0.05) #%>% pull(gene_symbol) %>% head(50)
top_7 <- opc_markers %>% dplyr::filter(cluster == 7 & p_val_adj < 0.05) #%>% pull(gene_symbol) %>% head(50)

top_4
top_7

# Extract gene lists for the above
genes_4 <- top_4$gene_symbol
genes_7 <- top_7$gene_symbol

head(genes_4)
head(genes_7)

# Do GO enrichment analysis for the 2 gene lists
library(clusterProfiler)
library(org.Mm.eg.db)

ego_4 <- enrichGO(gene = genes_4,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",  # Biological Process
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

ego_7 <- enrichGO(gene = genes_7,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)


dotplot(ego_4, showCategory = 12, title = "Cluster 4 GO Enrichment")
dotplot(ego_7, showCategory = 12, title = "Cluster 7 GO Enrichment")

# View top GO enrichment processes/pathways
head(ego_4$Description)
head(ego_7$Description)

# Create a mapping vector
gene_module_map <- setNames(all_markers_annotated$gene, all_markers_annotated$gene_symbol)

# Convert gene symbols to Seurat gene names for each list
top4_genes <- gene_module_map[unlist(top_4)]
top4_genes <- top4_genes[!is.na(top4_genes)]

top7_genes <- gene_module_map[unlist(top_7)]
top7_genes <- top7_genes[!is.na(top7_genes)]

# AddModuleScore() summarizes a set of gene expressions (i.e., markers) into a single score per cell.
# Add module scores separately
module_object <- AddModuleScore(object, features = list(top4_genes, top7_genes), name = "ClusterScore")

# This gives you: ClusterScore1 and ClusterScore2
# ClusterScore1 = how strongly a cell expresses Cluster 4 genes.
# ClusterScore2 = how strongly it expresses Cluster 7 genes.

module_object@meta.data

options(repr.plot.width = 15, repr.plot.height = 6)

VlnPlot(module_object, features = "ClusterScore1", group.by = "Original_Cluster") + ggtitle("Cluster 4 Markers")
VlnPlot(module_object, features = "ClusterScore2", group.by = "Original_Cluster") + ggtitle("Cluster 7 Markers")

# Extract the module scores and cluster assignments for all cells
# Extract module scores
score1 <- module_object@meta.data$ClusterScore1
score2 <- module_object@meta.data$ClusterScore2

# Extract cluster labels (replace with the actual cluster column name in your object)
clusters <- module_object@meta.data$Original_Cluster 

# Check lengths to confirm everything matches
length(score1)
length(score2)
length(clusters)

# Subset cells belonging to cluster 4 based on CHOIR clustering
cluster4_cells <- subset(module_object, subset = CHOIR_clusters_0.05 == "4")

# Subset cells belonging to cluster 7 based on CHOIR clustering
cluster7_cells <- subset(module_object, subset = CHOIR_clusters_0.05 == "7")

# Extract ClusterScore1 (Cluster 4 marker expression) for cluster 4 cells
score1_cluster4 <- cluster4_cells@meta.data$ClusterScore1

# Extract ClusterScore2 (Cluster 7 marker expression) for cluster 4 cells
score2_cluster4 <- cluster4_cells@meta.data$ClusterScore2

# Wilcoxon signed-rank test (paired) to check if cluster 4 cells express their own markers (score1) significantly differently than cluster 7 markers (score2)
wilcox.test(score1_cluster4, score2_cluster4, paired = TRUE)

# Extract ClusterScore1 (Cluster 4 marker expression) for cluster 7 cells
score1_cluster7 <- cluster7_cells@meta.data$ClusterScore1

# Extract ClusterScore2 (Cluster 7 marker expression) for cluster 7 cells
score2_cluster7 <- cluster7_cells@meta.data$ClusterScore2

# Wilcoxon signed-rank test (paired) to check if cluster 7 cells express their own markers (score2) significantly differently than cluster 4 markers (score1)
wilcox.test(score2_cluster7, score1_cluster7, paired = TRUE)

# Create df for ClusterScore1
df1 <- data.frame(
  Cluster = as.factor(clusters),
  ModuleScore = score1
)

# Create df for ClusterScore2
df2 <- data.frame(
  Cluster = as.factor(clusters),
  ModuleScore = score2
)

# Median for ClusterScore1
aggregate(ModuleScore ~ Cluster, data = df1, FUN = median)

# Median for ClusterScore2
aggregate(ModuleScore ~ Cluster, data = df2, FUN = median)


library(ggplot2)

# Combine data for easier plotting with a new column indicating the score type
df1$ScoreType <- "ClusterScore1"
df2$ScoreType <- "ClusterScore2"
combined_df <- rbind(df1, df2)

# Density plot
ggplot(combined_df[combined_df$Cluster %in% c("OPC"), ], aes(x = ModuleScore, fill = Cluster)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ScoreType) +
  theme_minimal() +
  labs(title = "Density plots of Module Scores by Cluster",
       x = "Module Score",
       y = "Density")

# Boxplot
ggplot(combined_df[combined_df$Cluster %in% c("OPC"), ], aes(x = Cluster, y = ModuleScore, fill = Cluster)) +
  geom_boxplot() +
  facet_wrap(~ScoreType) +
  theme_minimal() +
  labs(title = "Boxplots of Module Scores by Cluster",
       y = "Module Score")


# Extract metadata
meta <- module_object@meta.data

# make original cluster column a factor
meta$original_cluster <- as.factor(meta$Original_Cluster)

# Filter only OPC cells
opc_cells <- rownames(meta[meta$original_cluster == "OPC", ])

# Plot
ggplot(meta[opc_cells, ], aes(x = ClusterScore1, y = ClusterScore2)) +
  geom_point(alpha = 0.7, size = 1.8, color = "#1f78b4") +
  labs(
    title = "Module Scores in OPC Cells",
    x = "Cluster 4 Gene Module Score",
    y = "Cluster 7 Gene Module Score"
  ) +
  theme_minimal(base_size = 14) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40")


