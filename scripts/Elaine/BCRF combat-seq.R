

## ============================================================================
##  1: load BCRF Ep Seurat v5 object -> run ComBat-seq w batch correction (S3S v1 vs v2)
##  2: UMAP before vs after, colored by protocol and diagnosis
##  3: save batch-corrected outputs for inferCNV
## ============================================================================

library(Seurat)       # v5
library(sva)          # ComBat_seq
library(ggplot2)
library(patchwork)    
library(uwot)         # umap
library(matrixStats)  # rowVars


## ============================================================================
## CONFIG - edit
## ============================================================================
seurat_rds_path <- "/Users/eguo/Downloads/integrated BCRF set 1-10 & EN/RE-DONE integrate BCRF v1/RE-DONE Ep/RPCA/RE-DONE integrate BCRF v1_Ep.rds"
seurat_assay    <- "RNA"   # correct, raw counts is stored in the RNA assay
diagnosis_col   <- "dxs" 
protocol_col    <- "protocol"   # S3S v1, v2

# Drop ALH and BBD_RS before ComBat-seq
# ALH and BBD_RS are completely confounded with protocol — every sample is v1, none is v2
exclude_diagnoses <- c("ALH", "BBD_RS")  


output_dir <- "/Users/eguo/Downloads/integrated BCRF set 1-10 & EN/RE-DONE integrate BCRF v1/RE-DONE Ep/combat-seq"
corrected_matrix_out <- file.path(output_dir, "RE-DONE BCRF v1_Ep_combatseq_counts.tsv") #V
corrected_seurat_out <- file.path(output_dir, "RE-DONE BCRF v1_Ep_combatseq.rds") #V
design_out           <- file.path(output_dir, "RE-DONE BCRF v1_Ep_combatseq_design.tsv") #V
umap_plot_out        <- file.path(output_dir, "Ep_umap_before_after_combatseq.png")

n_top_var_genes <- 2000   # variable genes used for the UMAP QC plots (not for ComBat-seq itself)



## ============================================================================
## 1: Load Seurat v5 object & run ComBat-seq
## ============================================================================

seu <- readRDS(seurat_rds_path)
stopifnot(inherits(seu, "Seurat"))

## --- pull raw counts ---
counts <- as.matrix(LayerData(seu, assay = seurat_assay, layer = "counts"))
counts <- round(counts)     # ComBat_seq requires integer counts
counts[counts < 0] <- 0     # Enforce non-negative counts
message("Counts matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

## --- metadata ---
meta <- seu@meta.data

# Ensure counts and meta align to the same set/order of samples
# Keep only samples present in both counts and metadata (should normally be all of them), in the same order
common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples, drop = FALSE]
meta   <- meta[common_samples, , drop = FALSE]

## --- exclude diagnoses ---
keep <- !(meta[[diagnosis_col]] %in% exclude_diagnoses)
message("Excluding ", sum(!keep), " samples with diagnosis in: ",
        paste(exclude_diagnoses, collapse = ", "))
counts <- counts[, keep, drop = FALSE]
meta   <- meta[keep, , drop = FALSE]

## drop zero-variance genes (breaks PCA/ComBat-seq)
gene_var <- matrixStats::rowVars(counts)
counts <- counts[gene_var > 0, , drop = FALSE]

meta[[protocol_col]] <- factor(meta[[protocol_col]], levels = sort(unique(meta[[protocol_col]])))


message("Final matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")
print(table(meta[[protocol_col]]))
print(table(meta[[diagnosis_col]], meta[[protocol_col]]))



## --- run ComBat-seq ---
combatseq_counts <- ComBat_seq(
  counts = counts, # raw integer matrix (genes × samples)
  batch = meta[[protocol_col]],  # batch labels (protocol v1 vs v2)
  group = meta[[diagnosis_col]]  # biological group you want to preserve (diagnosis). 
  )

# group = dxs tells ComBat-seq to protect dxs-driven differences
# in both the mean and the dispersion when it fits the negative-binomial model 
# remove protocol effects, but keep differences driven by dxs




# ComBat_seq can return non-integers due to model-based adjustment, so:
combatseq_counts <- round(combatseq_counts)   # round to integers
combatseq_counts[combatseq_counts < 0] <- 0   # clamp negatives to 0
storage.mode(combatseq_counts) <- "integer"   # store as integer matrix

message("ComBat-seq done.")




## ============================================================================
## 2: UMAP before vs after ComBat-seq, colored by protocol and diagnosis
## ============================================================================
# standard scanpy/Seurat-style steps: 
# normalize → log → highly variable genes → scale → PCA → Euclidean distance on PCs → UMAP

compute_umap <- function(count_mat, n_top = n_top_var_genes, n_pcs = 20) {
  ## 1. library-size normalize (CPM)
  cpm_mat <- sweep(count_mat, 2, colSums(count_mat), "/") * 1e6
  logmat  <- log1p(cpm_mat)
  
  ## 2. select top variable genes
  gv <- matrixStats::rowVars(logmat)
  top_genes <- order(gv, decreasing = TRUE)[seq_len(min(n_top, length(gv)))]
  mat_sub <- t(logmat[top_genes, , drop = FALSE])   # samples x genes
  
  ## 3. scale (z-score) each gene -- standard pre-PCA step
  mat_scaled <- scale(mat_sub)
  mat_scaled[is.nan(mat_scaled)] <- 0   # guard against zero-variance genes post-subset
  
  ## 4. PCA
  n_pcs_use <- min(n_pcs, nrow(mat_scaled) - 1, ncol(mat_scaled))
  pca <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)  # already scaled above
  pcs <- pca$x[, seq_len(n_pcs_use), drop = FALSE]
  
  ## 5. UMAP on PCs, Euclidean distance
  set.seed(74)
  n_neighbors <- max(2, min(15, nrow(pcs) - 1))
  um <- uwot::umap(pcs, n_neighbors = n_neighbors, min_dist = 0.3, metric = "euclidean")
  
  data.frame(
    UMAP1 = um[, 1],
    UMAP2 = um[, 2],
    protocol  = meta[[protocol_col]],
    diagnosis = meta[[diagnosis_col]]
  )
}



plot_umap <- function(df, color_col, title, legend_title, palette = NULL) {
  p <- ggplot(df, aes(UMAP1, UMAP2, color = .data[[color_col]])) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = title, color = legend_title) +
    theme(plot.title = element_text(face = "bold"))
  if (!is.null(palette)) p <- p + scale_color_brewer(palette = palette)
  p
}

message("Computing UMAP: before ComBat-seq ...")
df_before <- compute_umap(counts)

message("Computing UMAP: after ComBat-seq ...")
df_after  <- compute_umap(combatseq_counts)

## 2x2 grid:
##   row 1 = BEFORE correction  (protocol | diagnosis)
##   row 2 = AFTER  correction  (protocol | diagnosis)
p_before_protocol  <- plot_umap(df_before, "protocol",  "Before ComBat-seq -- Protocol",  "Protocol")
p_before_diagnosis <- plot_umap(df_before, "diagnosis", "Before ComBat-seq -- Diagnosis", "Diagnosis")
p_after_protocol   <- plot_umap(df_after,  "protocol",  "After ComBat-seq -- Protocol",   "Protocol")
p_after_diagnosis  <- plot_umap(df_after,  "diagnosis", "After ComBat-seq -- Diagnosis",  "Diagnosis")

combined_plot <- (p_before_protocol | p_before_diagnosis) /
  (p_after_protocol  | p_after_diagnosis)

ggsave(umap_plot_out, combined_plot, width = 10, height = 10, dpi = 300)
message("Saved before/after UMAP comparison (protocol + diagnosis) to: ", umap_plot_out)




## ============================================================================
## 3: Save batch-corrected outputs for inferCNV
## ============================================================================

## corrected counts matrix (genes x samples) -> direct inferCNV input
write.table(
  data.frame(gene = rownames(combatseq_counts), combatseq_counts, check.names = FALSE),
  file = corrected_matrix_out, sep = "\t", quote = FALSE, row.names = FALSE
)
message("Saved batch-corrected counts matrix (tsv) to: ", corrected_matrix_out)

## replaces the "counts" layer in the RNA assay with combatseq_counts
## The Seurat object now stores corrected counts as its raw counts layer

seu_corrected <- subset(seu, cells = colnames(combatseq_counts))
seu_corrected <- seu_corrected[rownames(combatseq_counts), ]
seu_corrected <- SetAssayData(seu_corrected, assay = seurat_assay, layer = "counts",
                              new.data = combatseq_counts)
saveRDS(seu_corrected, corrected_seurat_out)
message("Saved batch-corrected Seurat object (rds) to: ", corrected_seurat_out)

## matching annotation/design file (sample_id, diagnosis, protocol, ...) for inferCNV

meta$cell_id <- rownames(meta)
meta <- meta[, c("cell_id", setdiff(names(meta), "cell_id"))] 
write.table(meta, design_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved matching design/annotation file to: ", design_out)




