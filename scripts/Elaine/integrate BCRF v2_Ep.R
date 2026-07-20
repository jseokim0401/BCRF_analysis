# from integrate BCRF v2.R
# from Claude/ final 

# In Seurat v5, when a merged object keeps its layers split by dataset, running NormalizeData() → FindVariableFeatures() → ScaleData() → RunPCA() on the merged object computes each of those steps per layer/per dataset automatically 
# Your code has the steps in the right order: preprocess → IntegrateLayers → JoinLayers


# integrate BCRF.R
# from TBCRCn, RAHBTn, COMETn
# Normalize, scale from BCRF set n.R 

# Combine v1 and v2 using seurat v5 rpca integration, did not do any other testing- Magda
# Combine v1 and v2, need batch correction, check w pca/ clustering - Siri
# S3S v2 has more reads


library(future)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

options(future.globals.maxSize = 1e10)  # ~10 GB
options(Seurat.object.assay.version = "v5")

# ---------------------------------------------------------------------------
# Load Seurat objects
# ---------------------------------------------------------------------------
BCRF_set1_3  <- readRDS('/Users/eguo/Downloads/BCRF RNA set 1-3, EN/BCRF RNA set 1-3 & EN.rds')
BCRF_set4_10 <- readRDS('/Users/eguo/Downloads/BCRF RNA set 4-10/BCRF RNA set 4-10.rds')



class(BCRF_set1_3[["RNA"]])
class(BCRF_set4_10[["RNA"]])

# If they're old-style Assay objects, convert first, e.g.:
#   BCRF_set1_3[["RNA"]]  <- as(BCRF_set1_3[["RNA"]],  "Assay5")
#   BCRF_set4_10[["RNA"]] <- as(BCRF_set4_10[["RNA"]], "Assay5")

# ---------------------------------------------------------------------------
# Merge: keeps counts/data layers split by dataset (v5)
# ---------------------------------------------------------------------------
merged_obj <- merge(BCRF_set1_3, y = BCRF_set4_10, add.cell.ids = c("set1_3", "set4_10"))

# Confirm these are already Assay5 objects, otherwise layers won't stay split after merge(), and IntegrateLayers will fail. 
merged_obj[["RNA"]] # should print Assay5
Layers(merged_obj) 

# Use group.by = "dataset" in DimPlot
merged_obj$dataset <- ifelse(grepl("^set1_3_", colnames(merged_obj)), "set1_3", "set4_10")

# Sanity check: confirm layers are actually split, and see their real names
merged_obj
Layers(merged_obj)



# ---------------------------------------------------------------------------
# Subset to Ep tissue only (from both datasets, post-merge)
# ---------------------------------------------------------------------------
table(merged_obj$tissue, merged_obj$dataset)

merged_obj <- subset(merged_obj, subset = tissue == "Ep")

# Subsetting a split-layer v5 object should preserve per-dataset layers, just
# filtered down - confirm this held, and check new sample counts per dataset
Layers(merged_obj)
table(merged_obj$dataset)





# ---------------------------------------------------------------------------
# Individual (per-dataset) preprocessing
# Because layers are split by dataset, these run per-dataset under the hood.
# Run PCA individually on each dataset prior to integration! [Reviewed- RPCA satija lab, last updated on Oct 31, 2023]
# ---------------------------------------------------------------------------
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

# Make sure k.weight is comfortably below the smaller count, or IntegrateLayers will error 
table(merged_obj$dataset)

# ---------------------------------------------------------------------------
# RPCA integration
# ---------------------------------------------------------------------------
merged_obj <- IntegrateLayers(
  object = merged_obj,
  method = RPCAIntegration,
  k.weight = 90,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

# Original unmodified data still resides in the 'RNA' assay.
# IntegrateLayers does NOT create a new "integrated" assay in the v5 workflow - it only adds "integrated.rpca" reduction to RNA. 

# Join layers back together for downstream analysis/visualization
merged_obj <- JoinLayers(merged_obj)

# ---------------------------------------------------------------------------
# Compare before vs. after integration
# ---------------------------------------------------------------------------
merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

p1 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "dataset") +
  ggtitle("Before Integration")
p2 <- DimPlot(merged_obj, reduction = "umap.rpca", group.by = "dataset") +
  ggtitle("After Integration (RPCA)")

p1 + p2


setwd("/Users/eguo/Downloads/integrated BCRF set 1-10 & EN/integrate BCRF v2_Ep")

ggsave("before v after integration.png", width = 12, height = 6)

# ---------------------------------------------------------------------------
# Clustering on the integrated reduction
# (No need to re-run ScaleData/RunPCA here - integrated.rpca already exists and is what clustering/UMAP should be based on.)
# ---------------------------------------------------------------------------
merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:30)

# Resolution = 2.5 is high and will produce many small clusters; typical cell-type clustering uses 0.4-1.2
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "rpca_clusters")


# ---------------------------------------------------------------------------
# Visualize
# ---------------------------------------------------------------------------






DimPlot(merged_obj, reduction = "umap.rpca", label = TRUE)
ggsave("umap.rpca.png", width = 8, height = 8)

merged_obj$tissue[is.na(merged_obj$tissue)] <- 'NA'
DimPlot(merged_obj, reduction = "umap.rpca", label = TRUE, group.by = "tissue")
ggsave("tissue.png", width = 8, height = 8)

DimPlot(merged_obj, reduction = "umap.rpca", group.by = "protocol")
ggsave("protocol.png", width = 8, height = 8)

DimPlot(merged_obj, reduction = "umap.rpca", group.by = "dxs")
ggsave("dxs.png", width = 8, height = 8)

DimPlot(merged_obj, reduction = "umap.rpca", group.by = "dx")
ggsave("dx.png", width = 12, height = 8)

n_patients <- length(unique(merged_obj$patient))
DimPlot(merged_obj, reduction = "umap.rpca", group.by = "patient") +
  ggtitle(paste("n patients = ", n_patients))
ggsave("patient.png", width = 12, height = 8)

DimPlot(merged_obj, reduction = "umap.rpca", group.by = "sequence batch number")
ggsave("sequence batch number.png", width = 12, height = 8)

# log2 of total counts (library size) per cell
merged_obj$log2counts <- log2(merged_obj$nCount_RNA + 1)
log2counts <- FeaturePlot(merged_obj, reduction = "umap.rpca", features = "log2counts", combine = TRUE)
print(log2counts)
ggsave("log2counts.png", width = 8, height = 8)

FeaturePlot(merged_obj, reduction = "umap.rpca", features = c("ESR1", "ERBB2", "MYC", "TP53", "PIK3CA"), cols = c("red", "blue"))
ggsave("genes.png", width = 8, height = 8)


saveRDS(merged_obj, file = "/Users/eguo/Downloads/integrated BCRF set 1-10 & EN/integrated_BCRF_set_1-10_&_EN_v2_Ep.rds")




#################
# FindAllMarkers run a simple cluster-vs-rest test; there's no argument to add covariates or a design formula. 
# If you want to control for patient, Seurat's wrapper can't do that — need to run DESeq2 directly with design = ~ patient + cluster 
# Since multiple samples in your data come from the same patient, use DESeq()


library(DESeq2)

Idents(merged_obj) <- "rpca_clusters"

# Pull raw counts + metadata from Seurat object 
counts_mat <- as.matrix(GetAssayData(merged_obj, assay = "RNA", layer = "counts"))
counts_mat <- round(counts_mat)  # DESeq2 requires integer counts

meta <- merged_obj@meta.data
meta$patient <- factor(as.character(merged_obj$patient))
meta$cluster <- factor(Idents(merged_obj))

# Check: every patient should span >1 cluster, or the design below will be unidentifiable (model matrix not full rank)
table(meta$patient, meta$cluster)

# Loop: one-vs-rest per cluster, correcting for patient each time
cluster_levels <- levels(meta$cluster)

cluster_markers_all <- lapply(cluster_levels, function(cl) {
  meta_bin <- meta
  meta_bin$cluster_bin <- factor(ifelse(meta_bin$cluster == cl, cl, "rest"),
                                 levels = c("rest", cl))
  
  dds_bin <- DESeqDataSetFromMatrix(countData = counts_mat,
                                    colData = meta_bin,
                                    design = ~ patient + cluster_bin)
  dds_bin <- DESeq(dds_bin)
  
  res_df <- as.data.frame(results(dds_bin, contrast = c("cluster_bin", cl, "rest")))
  res_df$gene <- rownames(res_df)
  res_df$cluster <- cl
  res_df
})

cluster_markers_all <- do.call(rbind, cluster_markers_all)
head(cluster_markers_all)

# Filter by adj pval and Log2FC
cluster_markers_filtered <- cluster_markers_all %>%
  filter(!is.na(padj), padj < 0.05)

top_genes_per_cluster <- cluster_markers_filtered %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  group_by(cluster) %>%
  slice_max(order_by = abs(log2FoldChange), n = 100) %>%
  ungroup()





