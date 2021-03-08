
library(magrittr)
library(tidyverse)
library(Seurat)
library(readxl)
library(cowplot)
library(colorblindr)
library(viridis)
library(harmony)
setwd(rprojroot::find_rstudio_root_file())

theme_cowplot2 <- function(...) {
  theme_cowplot(font_size = 12, ...) %+replace%
    theme(strip.background = element_blank(),
          plot.background = element_blank())
}
theme_set(theme_cowplot2())

coi <- "Ovarian.cancer.super"
cell_sort <- "CD45-"
louvain_resolution <- 0.3
nF <- 1000
pMT_lower <- 1

# coi <- "T.super"
# cell_sort <- "CD45+"
# louvain_resolution <- 0.2
# nF <- 500
# pMT_lower <- 1

# coi <- "Myeloid.super"
# cell_sort <- "CD45-"
# louvain_resolution <- 0.2
# nF <- 500
# pMT_lower <- 1

### load all data ---------------------------------
seu_obj <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/outs_pre/", coi, "_seurat_", louvain_resolution, ".rds"))
seu_obj_sub <- subset(seu_obj, subset = nFeature_RNA > 1000 & percent.mt > 1)
set.seed(123)
seu_obj_sample <- subset(seu_obj_sub, cells = sample(colnames(seu_obj_sub), 10000))

helper_f <- function(x) ifelse(is.na(x), "", x)
helper_f2 <- function(x) select(unnest(enframe(x, "subtype", "gene"), cols = gene), gene, subtype)

clrs <- yaml::read_yaml("/home/uhlitzf/spectrum_tme/_data/small/signatures/hgsc_v6_colors.yaml") %>% 
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

meta_tbl <- read_excel("_data/small/MSK SPECTRUM - Single cell RNA-seq_v6.xlsx", sheet = 2) %>% 
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"))

## Discarding bad qc cancer cells
## Cancer cells with less than one percent mito reads and less than 1000 expressed features were discarded.

## Re-preprocessing of high qc cancer cells ---------------------------------

preprocess_wrapper <- . %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F) %>%
  RunHarmony(group.by.vars = "patient_id") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(res = 0.1) %>%
  FindClusters(res = 0.2) %>%
  FindClusters(res = 0.3) %>%
  RunUMAP(reduction = "pca", dims = 1:20, reduction.name = "umappca", reduction.key = "umappca_") %>%
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umapharmony", reduction.key = "umapharmony_")

## remove patient-specific clusters

find_markers_parallel <- function(seu_x, resolution) {

  Idents(seu_x) <- seu_x[[resolution]] %>% {setNames(.[,1], rownames(.))}
  marker_list <- parallel::mclapply(unique(Idents(seu_x)),
                                    function(x) FindMarkers(seu_x, ident.1 = x, only.pos = T),
                                    mc.cores = 20)

  marker_tbl <- marker_list %>%
    setNames(unique(Idents(seu_x))) %>%
    lapply(as_tibble, rownames = "gene") %>%
    bind_rows(.id = "cluster")

  return(marker_tbl)

}

## sample subset ---------------------------------

seu_obj_sample <- preprocess_wrapper(seu_obj_sample)
write_rds(seu_obj_sample, paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc_sample.rds"))
seu_obj_sample <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc_sample.rds"))

# Idents(seu_obj_sample) <- seu_obj_sample$RNA_snn_res.0.2
# marker_tbl_02_sample <- as_tibble(FindAllMarkers(seu_obj_sample, only.pos = T))
# write_tsv(marker_tbl_02_sample, paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc_markers_02_sample.tsv"))

Idents(seu_obj_sample) <- seu_obj_sample$RNA_snn_res.0.3
marker_tbl_03_sample <- as_tibble(FindAllMarkers(seu_obj_sample, only.pos = T))
write_tsv(marker_tbl_03_sample, paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc_markers_03_sample.tsv"))



## full subset ---------------------------------

seu_obj_sub <- preprocess_wrapper(seu_obj_sub)

## remove patient-specific clusters

patient_specific_clusters <- FetchData(seu_obj_sub, c("patient_id", "RNA_snn_res.0.3")) %>% 
  as_tibble() %>% 
  group_by(patient_id, RNA_snn_res.0.3) %>% 
  tally %>% 
  group_by(RNA_snn_res.0.3) %>% 
  mutate(nrel = n/sum(n),
         patient_specific = nrel > 0.5) %>% 
  filter(patient_specific) %>% 
  pull(RNA_snn_res.0.3) %>% 
  as.character() %>% 
  unique

seu_obj_sub

write_rds(seu_obj_sub, paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc.rds"))
seu_obj_sub <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc.rds"))

# Idents(seu_obj_sub) <- seu_obj_sub$RNA_snn_res.0.2
# marker_tbl_02 <- as_tibble(FindAllMarkers(seu_obj_sub, only.pos = T))
# write_tsv(marker_tbl_02, paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc_markers_02.tsv"))

Idents(seu_obj_sub) <- seu_obj_sub$RNA_snn_res.0.3
marker_tbl_03 <- as_tibble(FindAllMarkers(seu_obj_sub, only.pos = T))
write_tsv(marker_tbl_03, paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/", coi, "_highqc_markers_03.tsv"))


