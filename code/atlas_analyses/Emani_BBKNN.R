# Install if needed
install.packages("Seurat")
install.packages("SeuratDisk")

library(Seurat)
library(SeuratDisk)

# Load the RDS file
obj <- readRDS("/home/mmiihkin/Downloads/GSE261983_integrated2.rds")

obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)

# Now run Leiden at multiple resolutions (algorithm=4 is Leiden in many Seurat versions)
obj <- FindClusters(
  obj,
  algorithm = 4,
  resolution = c(0.5, 1.0, 1.5, 2.0, 2.5),
  verbose = FALSE
)

# Check what cluster columns were created
grep("snn_res", colnames(obj@meta.data), value = TRUE)

# Save as h5Seurat
SaveH5Seurat(obj, filename = "GSE261983_integrated2.h5Seurat")

# Convert to h5ad
Convert("GSE261983_integrated2.h5Seurat", dest = "h5ad")
