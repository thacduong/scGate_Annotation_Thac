############ Load in packages ###############
lapply(c("dplyr","Seurat","HGNChelper","openxlsx", "patchwork", "SingleR", "scGate", "plotly", "purrr", "tidyverse"), library, character.only = T)

############ load gene set preparation function ###############

# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type excel file
#
# @params: path_to_db_file - DB file with cell types
# @cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#

gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}



############ load cell type annotation function ############

# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# sctype_score: calculate ScType scores and assign cell types
#
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells), 
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: gs - list of gene sets positively expressed in the cell type 
# @params: gs2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
}



############ User Input and Output ###############
# Input File
sc_RNAseq_dir = "single_cell_analysis_GSE205013/GSM6204109_P01/Data" 

# Input DB and study tissue (for sc_Type package)
db_ <- "ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# Project name
Prj = "GSE205013"

# Output directory
dir = getwd()

# Output File 1 cell-level meta
meta_dir = "GSM6204109_P01_meta_data_combined"

# Output File 2 cluster level annotatiion
cluster_annotation_dir = "GSM6204109_P01_cluster_annotation_percentage"



################ Load in seurat object ##################
# Load the sc-RNAseq dataset
pdac.data <- Read10X(data.dir = sc_RNAseq_dir)

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



################ scGate Annoatation ##################
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



################ SingleR annotation ##################
# Automated annotation using SingleR
reference_HPCA <- celldex::HumanPrimaryCellAtlasData()
singleR_result_HPCA <- SingleR(test = seu_obj@assays$RNA$data, ref = reference_HPCA, labels = reference_HPCA$label.main)
seu_obj$SingleR.labels.HPCA <- singleR_result_HPCA$labels
DimPlot(seu_obj, reduction = "umap", group.by = "SingleR.labels.HPCA", label = T, label.size = 3) + ggtitle("SingleR annotation HPCA")



################# sc_Type annotation #################
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seu_obj[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seu_obj[["RNA"]]$scale.data) else as.matrix(seu_obj[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. pbmc[["RNA"]]$counts or pbmc[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu_obj@meta.data[seu_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu_obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

seu_obj@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu_obj@meta.data$sctype_classification[seu_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seu_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')        



################# find the percentage of each cell types in each clusters #####################
# List of cell type columns to analyze
cell_type_columns <- c(
  "scGate_generic_Celltype", "scGate_HiTME_Celltype", "scGate_PBMC_Celltype",  
  "scGate_TME_broad_Celltype", "scGate_TME_HiRes_Celltype", "sctype_classification", "SingleR.labels.HPCA",  
  "scGate_CD4_TIL_Celltype", "scGate_CD8_TIL_Celltype"
)

# Initialize an empty list to store the results
cluster_percentage_list <- list()

# Loop through each cell type column
for (cell_type_col in cell_type_columns) {
  # Generate a df with cluster and cell type information
  cluster_celltype_df <- data.frame(
    cluster = seu_obj$RNA_snn_res.0.5, 
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
    summarise(!!paste0(cell_type_col, "_percentage") := paste0(cell_type, "(", round(percentage, 2), "%)", collapse = " "))
  
  # Store the result in the list
  cluster_percentage_list[[cell_type_col]] <- cluster_percentage_strings
}

# Merge all the results into a single data frame
cluster_percentage_strings_combined <- purrr::reduce(cluster_percentage_list, full_join, by = "cluster")

# View the final combined data frame
View(cluster_percentage_strings_combined)

# Generate cluster level annotation
write.csv(cluster_percentage_strings_combined, file = paste0(dir, "/", cluster_annotation_dir, ".csv"))

# Generate the combined metadata to a csv file
write.csv(seu_obj@meta.data, file = paste0(dir, "/", meta_dir, ".csv"))



# ##### Generate Umap #####
# # Metadata column to group
# group_col = "scGate_multi_generic"
# DimPlot(seu_obj, reduction = 'umap', group.by = group_col, label = T, label.size = 3) + ggtitle(group_col)
