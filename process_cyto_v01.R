## 21 June 2024 Foo Cheung 
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
library(htmltools)
library(future)
source("dev.cyto_functions_cite.R")
library("dplyr")
library("matrixStats")
library('tidyverse')
library(ggplot2)
library(cowplot)
library(patchwork)
library(harmony)
library(dsb,lib="./lib")
library(RSpectra)
library(ttservice, lib="./lib")
library(tidyseurat, lib="./lib")
library(paletteer, lib="./lib")
library(snakecase, lib="./lib")
library(janitor, lib="./lib")
library(ggprism, lib="./lib")
library(scCustomize, lib="./lib")
library(glmGamPoi)
library(R.utils) 
library(SeuratWrappers)
library(Azimuth) 
library(prismatic, lib="./lib")

# Function to process each dataset
process_dataset <- function(config_file) {
  # Load configuration from YAML file
  config <- yaml::read_yaml(config_file)
  
 # B1_US_data <- Read10X_h5(config$data_path)
 B1_US_data <- read_h5_files(config$lanes, config$rawdir, config$prefix)



 B1_US_SeuratObj <- create_seurat_objects(B1_US_data, config$lanes)
 B1_US_SeuratObj<-merge_seurat_objects(B1_US_SeuratObj,config$lanes)
 B1_US_SeuratObj <- add_metadata_and_filter(B1_US_SeuratObj,config$output_prefix)


DefaultAssay(B1_US_SeuratObj)<-"RNA"
B1_US_SeuratObj<-RunAzimuth(B1_US_SeuratObj, reference = "./pbmcref.SeuratData/azimuth/")
B1_US_SeuratObj[["RNA"]]<-split(B1_US_SeuratObj[["RNA"]], f = B1_US_SeuratObj$Lane)
obj<-B1_US_SeuratObj
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj,  dims=1:20,reduction = "pca", reduction.name = "umap.unintegrated")

obj$Lane<-as.factor(obj$Lane)


ggsave(paste(config$output_prefix,"_UMAP_L.pdf" , sep=""), width=length(obj$Lane  %>% unique()) + 1 * 3, height=length(obj$Lane  %>% unique()) + 1 * 3, plot=scCustomize::DimPlot_scCustom(obj, num_columns=1, reduction = "umap.unintegrated", group.by =  c("Lane", "predicted.celltype.l1" , "predicted.celltype.l2"), split.by="Lane", split_seurat = TRUE ))


ggsave(paste(config$output_prefix,"_UMAP.pdf" , sep=""), width=length(obj$Lane  %>% unique()) + 1 * 3, height=length(obj$Lane  %>% unique()) + 1 * 3, plot=scCustomize::DimPlot_scCustom(obj, num_columns=1, reduction = "umap.unintegrated", group.by =  c("Lane", "predicted.celltype.l1" , "predicted.celltype.l2"), split_seurat = TRUE ))




obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20)
obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, reduction.name = "harmony")

ggsave(paste(config$output_prefix,"_UMAP2_L.pdf" , sep=""), width=length(obj$Lane  %>% unique()) + 1 * 3, height=length(obj$Lane  %>% unique()) + 1 * 3, plot=scCustomize::DimPlot_scCustom(obj, num_columns= 1, reduction = "harmony", group.by =  c("Lane", "predicted.celltype.l1", "predicted.celltype.l2", "harmony_clusters"), split.by="Lane", split_seurat = TRUE ))


ggsave(paste(config$output_prefix,"_UMAP2.pdf" , sep=""), width=length(obj$Lane  %>% unique()) + 1 * 3, height=length(obj$Lane  %>% unique()) + 1 * 3, plot=scCustomize::DimPlot_scCustom(obj, num_columns= 1, reduction = "harmony", group.by =  c("Lane", "predicted.celltype.l1", "predicted.celltype.l2", "harmony_clusters"), split_seurat = TRUE ))



  saveRDS(obj, paste(config$output_prefix, "_ref.rds", sep = "")) 
}

process_dataset("config_1.yaml")
