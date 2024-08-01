library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(scGate)
library(plotly)
library(purrr)



##### User Input and Output #####
# Input File
sc_RNAseq_dir = "/single_cell_analysis_GSE205013/GSM6204110_P02/Data" 

# Project name
Prj = "GSE205013"

# Output directory
dir = getwd()

# Output File 1 cell-level meta
meta_dir = "GSM6204110_P02_meta_data_combined"

# Output File 2 cluster level annotatiion
cluster_annotation_dir = "GSM6204110_P02_cluster_annotation_percentage"



##### Load in seurat object #####
# Load the sc-RNAseq dataset
pdac.data <- Read10X(data.dir = paste0(dir, sc_RNAseq_dir))

# Initialize the Seurat object with the raw (non-normalized data).
seu_obj <- CreateSeuratObject(counts = pdac.data, project = Prj, min.cells = 3, min.features = 200)
seu_obj

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalizing the data
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- NormalizeData(seu_obj)

# Identification of highly variable features
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu_obj), 10)

# Scailing the data
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)

# Linear dimensional reduction
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
# Examine and visualize seurat object results a few different ways
print(seu_obj[["pca"]], dims = 1:5, nfeatures = 5)

# Cluster the cells 
seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(seu_obj), 5)

# Non-linear dimensional reduction 
seu_obj <- RunUMAP(seu_obj, dims = 1:10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
seu_obj.markers <- FindAllMarkers(seu_obj, only.pos = TRUE)
seu_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)



##### scGate Annoatation #####
# List of model names and corresponding scGate models
scGate_models_DB <-  get_scGateDB()
models <- list(
  CD4_TIL = scGate_models_DB$human$CD4_TIL,
  CD8_TIL = scGate_models_DB$human$CD8_TIL,
  generic = scGate_models_DB$human$generic,
  HiTME = scGate_models_DB$human$HiTME,
  PBMC = scGate_models_DB$human$PBMC,
  TME_broad = scGate_models_DB$human$TME_broad,
  TME_HiRes = scGate_models_DB$human$TME_HiRes
)

# Create a list to store new columns
new_columns_list <- list()

# Loop through each model
for (model_name in names(models)) {
  # Get the current model
  my_scGate_model <- models[[model_name]]
  
  # Copy the original Seurat object to avoid modification
  seu_obj_copy <- seu_obj
  
  # Apply scGate to the copied Seurat object
  seu_obj_processed <- scGate(seu_obj_copy, model = my_scGate_model, ncores = 4)
  
  # Assign the processed Seurat object to a variable named according to the model name
  assign(paste0("seu_obj_", model_name), seu_obj_processed)
  
  # Extract metadata
  meta_data <- seu_obj_processed@meta.data
  
  ##### Process "is.pure_" Columns #####
  # Identify columns that start with "is.pure_"
  pure_columns <- grep("^is.pure_", colnames(meta_data), value = TRUE)
  
  # Initialize the new column to store the features
  new_col_name <- paste0("scGate_", model_name, "_Celltype")
  meta_data[[new_col_name]] <- NA
  
  # Loop through each "is.pure_" column and update the new column
  for (col in pure_columns) {
    feature <- sub("^is.pure_", "", col)
    pure_barcodes <- rownames(meta_data)[meta_data[[col]] == "Pure"]
    for (barcode in pure_barcodes) {
      if (is.na(meta_data[barcode, new_col_name])) {
        meta_data[barcode, new_col_name] <- feature
      } else {
        meta_data[barcode, new_col_name] <- paste(meta_data[barcode, new_col_name], feature, sep = " ")
      }
    }
  }
  
  # Store the new column in the list
  new_columns_list[[new_col_name]] <- meta_data[, new_col_name, drop = F]
  
  ##### Process "scGate_multi" Column #####
  if ("scGate_multi" %in% colnames(meta_data)) {
    multi_col_name <- paste0("scGate_multi_", model_name)
    meta_data[[multi_col_name]] <- meta_data[["scGate_multi"]]
    
    # Store the renamed "scGate_multi" column in the list
    new_columns_list[[multi_col_name]] <- meta_data[, multi_col_name, drop = F]
  }
  
  # Save the metadata to a CSV file
  write.csv(meta_data, file = paste0(dir, "/meta_data_", model_name, ".csv"))
  View(meta_data)
}

# Combine the original metadata with the new columns
meta_data_combined <- seu_obj@meta.data

for (new_col in new_columns_list) {
  meta_data_combined <- merge(meta_data_combined, new_col, by = "row.names", all.x = TRUE)
  rownames(meta_data_combined) <- meta_data_combined$Row.names
  meta_data_combined$Row.names <- NULL
}

# Update the original Seurat object with the combined metadata
seu_obj@meta.data <- meta_data_combined
View(seu_obj@meta.data)



##### Cluster annotation #####
# Automated annotation using SingleR
reference_HPCA <- celldex::HumanPrimaryCellAtlasData()
singleR_result_HPCA <- SingleR(test = seu_obj@assays$RNA$data, ref = reference_HPCA, labels = reference_HPCA$label.main)
seu_obj$SingleR.labels.HPCA <- singleR_result_HPCA$labels
Idents(seu_obj) = "SingleR.labels.HPCA"
DimPlot(seu_obj, reduction = "umap", label = T, label.size = 3) + ggtitle("SingleR annotation HPCA")



##### find the percentage of each cell types in each clusters #####
# List of cell type columns to analyze
cell_type_columns <- c(
  "scGate_generic_Celltype", "scGate_HiTME_Celltype", "scGate_PBMC_Celltype",  
  "scGate_TME_broad_Celltype", "scGate_TME_HiRes_Celltype", "SingleR.labels.HPCA",  
  "scGate_CD4_TIL_Celltype", "scGate_CD8_TIL_Celltype"
)

# Initialize an empty list to store the results
cluster_percentage_list <- list()

# Loop through each cell type column
for (cell_type_col in cell_type_columns) {
  # Generate a df with cluster and cell type information
  cluster_celltype_df <- data.frame(
    cluster = Idents(seu_obj), 
    cell_type = seu_obj[[cell_type_col]]
  )
  cluster_celltype_df
  
  # Rename the cell_type column if it's not named correctly
  colnames(cluster_celltype_df)[2] <- "cell_type"
  
  # Convert NA values to a string so they are treated as a category
  cluster_celltype_df$cell_type <- ifelse(is.na(cluster_celltype_df$cell_type), "NA", cluster_celltype_df$cell_type)
  
  # Generate a table of cluster vs cell type counts 
  cluster_celltype_count <- table(cluster_celltype_df$cluster, cluster_celltype_df$cell_type)
  
  # Convert counts to percentages 
  cluster_celltype_percentage <- prop.table(cluster_celltype_count, margin = 1) * 100
  
  # Convert to df for easier viewing
  cluster_celltype_df <- as.data.frame(cluster_celltype_percentage)
  colnames(cluster_celltype_df) <- c('cluster', 'cell_type', 'percentage')
  
  # Cell type percentage grouped by cluster
  cluster_percentage_strings <- cluster_celltype_df %>% 
    group_by(cluster) %>% 
    arrange(desc(percentage)) %>% 
    summarise(!!paste0(cell_type_col, "percentage") := paste0(cell_type, "(", round(percentage, 2), "%)", collapse = " "))
  
  # Store the result in the list
  cluster_percentage_list[[cell_type_col]] <- cluster_percentage_strings
}

# Merge all the results into a single data frame
cluster_percentage_strings_combined <- reduce(cluster_percentage_list, full_join, by = "cluster")

# View the final combined data frame
View(cluster_percentage_strings_combined)

# Generate cluster level annotation
write.csv(cluster_percentage_strings_combined, file = paste0(dir, "/", cluster_annotation_dir, ".csv"))

# Generate the combined metadata to a csv file
write.csv(seu_obj@meta.data, file = paste0(dir, "/", meta_dir, ".csv"))


##### Generate Umap #####
# Metadata column to group
group_col = "scGate_multi_generic"
DimPlot(seu_obj, reduction = 'umap', group.by = group_col, label = T, label.size = 3) + ggtitle(group_col)
