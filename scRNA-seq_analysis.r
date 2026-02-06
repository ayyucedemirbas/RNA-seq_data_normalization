if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("glmGamPoi", quietly = TRUE))
    BiocManager::install("glmGamPoi")

if (!requireNamespace("Seurat", quietly = TRUE))
    install.packages("Seurat")

if (!requireNamespace("sctransform", quietly = TRUE))
    install.packages("sctransform")

if (!requireNamespace("Matrix", quietly = TRUE))
    install.packages("Matrix")

if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)

data_dir <- "pbmc3k_data"
if (!dir.exists(data_dir)) dir.create(data_dir)

url <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dest_file <- file.path(data_dir, "pbmc3k_filtered_gene_bc_matrices.tar.gz")

if (!file.exists(dest_file)) {
  message("Downloading PBMC 3k dataset (~7MB)...")
  download.file(url, destfile = dest_file)
}

untar(dest_file, exdir = data_dir)

matrix_dir <- file.path(data_dir, "filtered_gene_bc_matrices/hg19")

pbmc_data <- Read10X(data.dir = matrix_dir)

seurat_object <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)

message(paste("Loaded:", nrow(seurat_object), "genes x", ncol(seurat_object), "cells"))

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(seurat_object, features = "nFeature_RNA") +
  geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
  ggtitle("Unique Genes per Cell")

VlnPlot(seurat_object, features = "nCount_RNA") +
  ggtitle("Total UMI Counts per Cell")

VlnPlot(seurat_object, features = "percent.mt") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  ggtitle("Mitochondrial %")

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

seurat_object <- SCTransform(
  seurat_object,
  method = "glmGamPoi",
  vars.to.regress = "percent.mt",
  verbose = TRUE,
  vst.flavor = "v2"
)

sct_features <- as.data.frame(seurat_object[["SCT"]]@meta.features)

print(colnames(sct_features))

raw_counts <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")

vst_out <- sctransform::vst(
  raw_counts,
  method = "glmGamPoi",
  vst.flavor = "v2",
  verbosity = 1
)

plot_data <- vst_out$gene_attr

x_col <- ifelse("gmean" %in% colnames(plot_data), "gmean", "geometric_mean")

ggplot(plot_data, aes(x = .data[[x_col]], y = residual_variance)) +
  geom_point(alpha = 0.2, size = 0.5) +
  scale_x_log10() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  labs(
    title = "Variance Stabilization (Negative Binomial Fit)",
    subtitle = "Red line (1.0) = Theoretical variance of the NB model",
    x = "Geometric Mean Expression",
    y = "Residual Variance"
  ) +
  theme_bw()

seurat_object <- RunPCA(seurat_object, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
seurat_object <- FindClusters(seurat_object, verbose = FALSE)

DimPlot(seurat_object, label = TRUE) + ggtitle("PBMC 3k, Processed with sctransform (NB)")

DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
seurat_object <- ScaleData(seurat_object, verbose = FALSE)

all_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(dplyr)
top_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

print(top_markers)

# MS4A1 = B Cells
# GNLY  = NK Cells
# CD3E  = T Cells
# CD14  = CD14+ Monocytes
# FCGR3A = FCGR3A+ Monocytes
# PPBP  = Platelets

FeaturePlot(seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCGR3A", "PPBP"))

num_clusters <- length(levels(seurat_object))
print(paste("Clusters: ", num_clusters))

cluster_ids <- c(
  "Naive CD4 T",
  "CD14+ Mono",
  "Memory CD4 T",
  "B Cells",
  "CD8 T",
  "NK Cells",
  "FCGR3A+ Mono",
  "Dendritic DC",
  "T Cell Mix",
  "Platelets",
  "Cluster 10",
  "Cluster 11"
)

names(cluster_ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, cluster_ids)

DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("PBMC 3k Annotated (NB Model)")
