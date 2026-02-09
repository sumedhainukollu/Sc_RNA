library(Seurat)

sc <- readRDS("results/sc_qc_filtered.rds")

# Normalization
sc <- NormalizeData(sc, normalization.method = "LogNormalize",
                    scale.factor = 10000)

# Variable genes
sc <- FindVariableFeatures(sc, selection.method = "vst",
                            nfeatures = 2000)

# Scaling
sc <- ScaleData(sc)

# PCA
sc <- RunPCA(sc, features = VariableFeatures(sc))

# Elbow plot
ElbowPlot(sc, ndims = 30)

saveRDS(sc, "results/sc_pca.rds")
