# Single-cell RNA-seq pipeline ------------------------------------------------

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cellmarkeraccordion)

# Parameters -----------------------------------------------------------------

MIN_CELLS <- 3
MIN_FEATURES <- 200
MAX_FEATURES <- 6000
MAX_MT <- 20
N_VARIABLE_FEATURES <- 2000
N_PCS <- 20
CLUSTER_RESOLUTION <- 0.5

# 1. Load data ----------------------------------------------------------------

load("YKG.rds")

YKG <- CreateSeuratObject(
  counts = counts,
  project = "YKG",
  min.cells = MIN_CELLS,
  min.features = MIN_FEATURES
)

# 2. Quality control ----------------------------------------------------------

YKG[["percent.mt"]] <- PercentageFeatureSet(
  YKG,
  pattern = "^MT-"
)

VlnPlot(
  YKG,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(YKG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(YKG, feature1 = "nCount_RNA", feature2 = "percent.mt")

YKG <- subset(
  YKG,
  subset =
    nFeature_RNA > MIN_FEATURES &
    nFeature_RNA < MAX_FEATURES &
    percent.mt < MAX_MT
)

# 3. Preprocessing ------------------------------------------------------------

YKG <- NormalizeData(YKG)

YKG <- FindVariableFeatures(
  YKG,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES
)

VariableFeaturePlot(YKG)

YKG <- ScaleData(
  YKG,
  features = rownames(YKG)
)

# 4. Cell cycle scoring -------------------------------------------------------

YKG <- CellCycleScoring(
  YKG,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes
)

# 5. Dimensionality reduction -------------------------------------------------

YKG <- RunPCA(
  YKG,
  npcs = 50
)

ElbowPlot(YKG)

pcs <- 1:N_PCS

YKG <- RunUMAP(
  YKG,
  dims = pcs
)

# 6. Clustering ---------------------------------------------------------------

YKG <- FindNeighbors(
  YKG,
  dims = pcs
)

YKG <- FindClusters(
  YKG,
  resolution = CLUSTER_RESOLUTION
)

DimPlot(
  YKG,
  reduction = "umap",
  label = TRUE
)

# 7. Marker gene identification ----------------------------------------------

markers <- FindAllMarkers(
  YKG,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

DoHeatmap(
  YKG,
  features = unique(top10$gene)
)

# 8. SingleR annotation -------------------------------------------------------

ref <- HumanPrimaryCellAtlasData()

pred.cluster <- SingleR(
  test = GetAssayData(YKG, assay = "RNA", layer = "data"),
  ref = ref,
  labels = ref$label.main,
  clusters = YKG$seurat_clusters,
  de.method = "t"
)

cluster.labels <- pred.cluster$labels

YKG$SingleR_celltype <- cluster.labels[
  as.character(YKG$seurat_clusters)
]

DimPlot(
  YKG,
  reduction = "umap",
  group.by = "SingleR_celltype",
  label = TRUE,
  repel = TRUE
)

plotScoreHeatmap(pred.cluster)

# 9. CellMarkerAccordion annotation ------------------------------------------

YKG <- accordion(
  YKG,
  assay = "RNA",
  species = "Human",
  annotation_resolution = "cluster",
  max_n_marker = 30,
  include_detailed_annotation_info = TRUE,
  plot = TRUE
)

DimPlot(
  YKG,
  reduction = "umap",
  group.by = "accordion_per_cluster",
  label = TRUE,
  repel = TRUE
)

# 10. Cell type composition ---------------------------------------------------

singleR_counts <- table(YKG$SingleR_celltype)
singleR_percent <- round(prop.table(singleR_counts) * 100, 2)

singleR_counts
singleR_percent

# 11. Compare annotations -----------------------------------------------------

table(
  SingleR = YKG$SingleR_celltype,
  Accordion = YKG$accordion_per_cluster
)

# 12. Marker validation -------------------------------------------------------

FeaturePlot(
  YKG,
  features = c("CD3D", "CD3E", "MS4A1", "CD79A", "NKG7", "GNLY", "LYZ")
)

FeaturePlot(
  YKG,
  features = c("MKI67", "TOP2A", "PCNA")
)

FeaturePlot(
  YKG,
  features = c("COL2A1", "ACAN", "COL9A1")
)

# 13. Save results ------------------------------------------------------------

saveRDS(YKG, file = "YKG_processed_seurat_object.rds")
write.csv(markers, file = "YKG_cluster_markers.csv", row.names = FALSE)
write.csv(top10, file = "YKG_top10_markers.csv", row.names = FALSE)