
# ============================================================
# PROJECT 1: SINGLE-CELL RNA-SEQ ANALYSIS
# Dataset: GSE132771
# Sample: NML1
# ============================================================

# STEP 1: LOAD AND UNDERSTAND THE RAW DATA

# 1.1 Load required packages
library(Seurat)
library(Matrix)
library(ggplot2)

# 1.2 Load the raw gene-expression matrix
nml1 <- Read10X(
  data.dir = "GSE132771_RAW/NML1"
)

# 1.3 Check the dimensions of the matrix
dim(nml1)

# 1.4 Inspect the first five genes and first five cells
nml1[1:5, 1:5]


# ============================================================
# STEP 2: CREATE THE SEURAT OBJECT
# ============================================================

seu_nml1 <- CreateSeuratObject(
  counts = nml1,
  project = "GSE132771_NML1",
  min.cells = 0,
  min.features = 0
)

# Display the Seurat object summary
seu_nml1

# Inspect the first six cells
head(seu_nml1@meta.data)

# Summarize cell-level measurements
summary(
  seu_nml1@meta.data[
    ,
    c("nCount_RNA", "nFeature_RNA")
  ]
)

# ============================================================
# STEP 3: QUALITY-CONTROL ASSESSMENT
# ============================================================

# 3.1 Calculate mitochondrial percentages
seu_nml1[["percent.mt"]] <- PercentageFeatureSet(
  seu_nml1,
  pattern = "^MT-"
)

# 3.2 Summarize QC measurements
summary(
  seu_nml1@meta.data[
    ,
    c("nCount_RNA", "nFeature_RNA", "percent.mt")
  ]
)

# 3.3 Create QC violin plots using raw counts
qc_plot <- VlnPlot(
  seu_nml1,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1,
  layer = "counts"
)

# Display the plots
qc_plot

# 3.4 Create the figures folder if needed
dir.create("figures", showWarnings = FALSE)

# 3.5 Save the QC plots
ggsave(
  filename = "figures/quality_control_violin_plots.png",
  plot = qc_plot,
  width = 12,
  height = 5,
  dpi = 300
)


# ============================================================
# STEP 4: FILTER CELLS
# ============================================================

# 4.1 Filter cells using the initial QC criteria
seu_nml1_filtered <- subset(
  seu_nml1,
  subset = nFeature_RNA >= 200 & percent.mt < 10
)

# 4.2 Count cells before and after filtering
cells_before <- ncol(seu_nml1)
cells_after <- ncol(seu_nml1_filtered)

# 4.3 Calculate the number and percentage of cells removed
cells_removed <- cells_before - cells_after
percent_removed <- (cells_removed / cells_before) * 100

# 4.4 Display the filtering results
cat("Cells before filtering:", cells_before, "\n")
cat("Cells after filtering:", cells_after, "\n")
cat("Cells removed:", cells_removed, "\n")
cat("Percentage removed:", round(percent_removed, 2), "%\n")

# 4.5 Save the filtered Seurat object
saveRDS(
  seu_nml1_filtered,
  file = "seu_nml1_filtered.rds"
)

# ============================================================
# STEP 5: NORMALIZATION AND HIGHLY VARIABLE GENES
# ============================================================

# 5.1 Normalize gene expression
seu_nml1_filtered <- NormalizeData(
  seu_nml1_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# 5.2 Identify the 2,000 highly variable genes
seu_nml1_filtered <- FindVariableFeatures(
  seu_nml1_filtered,
  selection.method = "vst",
  nfeatures = 2000
)

# 5.3 Count the selected variable genes
number_variable_genes <- length(
  VariableFeatures(seu_nml1_filtered)
)

cat(
  "Number of highly variable genes:",
  number_variable_genes,
  "\n"
)

# 5.4 Display the first 10 highly variable genes
head(
  VariableFeatures(seu_nml1_filtered),
  10
)

# 5.5 Save the filtered and normalized Seurat object
saveRDS(
  seu_nml1_filtered,
  file = "seu_nml1_filtered.rds"
)


# ============================================================
# STEP 6: DIMENSIONALITY REDUCTION AND CLUSTERING
# ============================================================

# 6.1 Scale the 2,000 highly variable genes
seu_nml1_filtered <- ScaleData(
  seu_nml1_filtered,
  features = VariableFeatures(seu_nml1_filtered)
)

# 6.2 Check the dimensions of the scaled expression matrix
dim(seu_nml1_filtered[["RNA"]]$scale.data)


# 6.3 Run principal component analysis (PCA)
seu_nml1_filtered <- RunPCA(
  seu_nml1_filtered,
  features = VariableFeatures(seu_nml1_filtered)
)

# 6.4 Check the dimensions of the PCA results
dim(
  Embeddings(
    seu_nml1_filtered,
    reduction = "pca"
  )
)
# 6.5 Examine the principal components
elbow_plot <- ElbowPlot(
  seu_nml1_filtered,
  ndims = 50
)

# Display the elbow plot
elbow_plot

# 6.6 Save the elbow plot
ggsave(
  filename = "figures/pca_elbow_plot.png",
  plot = elbow_plot,
  width = 8,
  height = 5,
  dpi = 300
)
# 6.7 Find neighboring cells using the first 20 PCs
seu_nml1_filtered <- FindNeighbors(
  seu_nml1_filtered,
  dims = 1:20
)

# 6.8 Check the resulting Seurat object
seu_nml1_filtered

# 6.9 Identify cell clusters
seu_nml1_filtered <- FindClusters(
  seu_nml1_filtered,
  resolution = 0.5
)

# 6.10 Count cells in each cluster
table(Idents(seu_nml1_filtered))

# 6.11 Count the number of clusters
length(
  levels(Idents(seu_nml1_filtered))
)

# 6.12 Calculate UMAP coordinates
seu_nml1_filtered <- RunUMAP(
  seu_nml1_filtered,
  dims = 1:20
)

# 6.13 Plot the cells by cluster
umap_clusters <- DimPlot(
  seu_nml1_filtered,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

# Display the UMAP plot
umap_clusters

# 6.14 Save the UMAP plot
ggsave(
  filename = "figures/umap_clusters.png",
  plot = umap_clusters,
  width = 9,
  height = 7,
  dpi = 300
)

# 6.15 Save the Seurat object after clustering and UMAP
saveRDS(
  seu_nml1_filtered,
  file = "seu_nml1_filtered.rds"
)


# ============================================================
# STEP 7: IDENTIFY MARKER GENES
# ============================================================

# 7.1 Find marker genes for each cluster
cluster_markers <- FindAllMarkers(
  seu_nml1_filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 7.2 Inspect the marker-gene results
head(cluster_markers)

# 7.3 Count the marker genes identified per cluster
table(cluster_markers$cluster)


# 7.4 Select the top 10 marker genes per cluster
library(dplyr)

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  arrange(cluster, desc(avg_log2FC))

# 7.5 Display the top markers
print(
  top10_markers[, c(
    "cluster",
    "gene",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
  )],
  n = 180
)

# 7.6 Save the marker table
write.csv(
  top10_markers,
  file = "figures/top10_markers_per_cluster.csv",
  row.names = FALSE
)
#=============================================================================#
# 7.7 Keep statistically significant marker genes
significant_markers <- cluster_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25
  )

# 7.8 Count significant markers per cluster
table(significant_markers$cluster)

# 7.9 Save the significant marker table
write.csv(
  significant_markers,
  file = "figures/significant_cluster_markers.csv",
  row.names = FALSE
)


# 7.10 Visualize selected cell-type marker genes
genes_to_check <- c(
  "MARCO", "LYZ",       # Myeloid / macrophage
  "SFTPC", "AGER",      # Alveolar epithelial
  "CLDN5", "PECAM1",    # Endothelial
  "DCN", "LUM",          # Fibroblast / stromal
  "ACTA2", "MYH11",     # Smooth muscle
  "TPSB2", "CPA3",      # Mast cells
  "MKI67", "TOP2A",     # Cycling cells
  "JCHAIN", "MZB1"      # Plasma cells
)

# Keep only genes present in the dataset
genes_to_check <- genes_to_check[
  genes_to_check %in% rownames(seu_nml1_filtered)
]

# Create a dot plot across clusters
marker_dotplot <- DotPlot(
  seu_nml1_filtered,
  features = genes_to_check,
  group.by = "seurat_clusters"
) + RotatedAxis()

# Display the plot
marker_dotplot

# Save the plot
ggsave(
  filename = "figures/selected_marker_dotplot.png",
  plot = marker_dotplot,
  width = 13,
  height = 7,
  dpi = 300
)

# 7.13 Create provisional cell-type annotations

cluster_annotations <- data.frame(
  cluster = as.character(0:17),
  cell_type = c(
    "Macrophages",
    "Alveolar type 2 epithelial cells",
    "Endothelial cells",
    "Macrophages",
    "Fibroblasts / stromal cells",
    "Inflammatory monocytes",
    "Myeloid cells, subtype unresolved",
    "Alveolar type 1 epithelial cells",
    "Macrophages, subtype unresolved",
    "Endothelial cells",
    "Fibroblasts / stromal cells",
    "Smooth muscle cells",
    "Proliferating cells",
    "Endothelial cells, identity to verify",
    "Mast cells",
    "Lymphatic endothelial cells (provisional)",
    "Plasma cells",
    "SPP1-expressing myeloid cells"
  )
)

# Cluster 13 requires further investigation.
# Both endothelial and macrophage-associated markers are expressed.
# Its cell identity remains provisional; possible explanations include
# mixed cell populations, doublets, or ambient RNA.

print(cluster_annotations)

# 7.14 Save average expression of selected marker genes by cluster

marker_expression <- AverageExpression(
  seu_nml1_filtered,
  features = genes_to_check,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)$RNA

marker_expression_table <- data.frame(
  gene = rownames(marker_expression),
  marker_expression,
  check.names = FALSE
)

write.csv(
  marker_expression_table,
  "figures/marker_expression_by_cluster.csv",
  row.names = FALSE
)

write.csv(
  cluster_annotations,
  "figures/cluster_annotations.csv",
  row.names = FALSE
)

# 7.15 Visualize QC metrics across clusters

cluster_qc_plot <- VlnPlot(
  seu_nml1_filtered,
  features = c("nCount_RNA", "nFeature_RNA"),
  group.by = "seurat_clusters",
  pt.size = 0
)

# Save the cluster QC plot
ggsave(
  filename = "figures/cluster_qc_violin.png",
  plot = cluster_qc_plot,
  width = 12,
  height = 6,
  dpi = 300
)

# 7.16 Save the final Seurat object

saveRDS(
  seu_nml1_filtered,
  file = "seu_nml1_filtered.rds"
)