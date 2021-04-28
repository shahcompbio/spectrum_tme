library(magrittr)
library(tidyverse)
library(readxl)
library(Seurat)
library(cowplot)

#setwd(rstudioapi::getActiveProject())
setwd("/home/uhlitzf/spectrum_tme/")
# setwd("/home/uhlitzf/spectrum_tme")
ncores <- 40

source("_src/global_vars.R")

# seu_merge <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/outs_pre/integrated_seurat.rds")
# seu_snn <- seu_merge@graphs$RNA_snn
# rm(seu_merge)
# write_rds(seu_snn, "/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/cohort_snn_matrix.rds")
seu_snn <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/cohort_snn_matrix.rds")
seu_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/outs_pre/cells.tsv") %>% 
  left_join(meta_tbl, by = "sample") %>% 
  filter(therapy == "pre-Rx") %>% 
  mutate(cell_type = str_replace_all(cell_type, "\\.", " "))

seu_snn <- seu_snn[seu_tbl$cell_id, seu_tbl$cell_id]
rownames(seu_snn) <- str_sub(rownames(seu_snn), 13, 15)

get_bins <- function(x, y) split(x, ceiling(seq_along(x)/y))
snn_bins <- get_bins(1:ncol(seu_snn), 20000)
snn_bins_multicore <- lapply(snn_bins, get_bins, 20000/40)

pid_freq <- select(seu_tbl, cell_id, cell_type, patient_id_short) %>% 
  group_by(cell_type, patient_id_short) %>% 
  tally %>% 
  group_by(cell_type) %>% 
  mutate(nrel = n/sum(n)) %>% 
  ungroup %>% 
  select(-n)

get_neighbours <- function(i) {
  seu_snn[,i] %>% 
    .[.!=0] %>% 
    sort(decreasing = T) %>% 
    enframe("patient_id_short2", "value")
}

get_pm_score <- function(cell_index_range) {
  
  snn_tbl <- lapply(cell_index_range, get_neighbours) %>% 
    setNames(colnames(seu_snn)[cell_index_range]) %>% 
    bind_rows(.id = "cell_id") %>% 
    mutate(patient_id_short = str_sub(cell_id, 13, 15))
  
  pm_score <- snn_tbl %>% 
    left_join(select(seu_tbl, cell_id, cell_type), by = "cell_id") %>% 
    left_join(pid_freq, by = c("patient_id_short", "cell_type")) %>% 
    group_by(cell_id) %>% 
    # slice(1:30) %>% 
    mutate(nexp = nrel * length(cell_id),
           nobs = sum(patient_id_short == patient_id_short2),
           score = log2(nobs/(nexp+1))) %>% 
    distinct(cell_id, cell_type, nexp, nobs, score)
  
  return(pm_score)
  
}

pb <- txtProgressBar(min = 0, max = length(snn_bins_multicore), style = 3)
pm_score_list <- list()
for (i in 1:length(snn_bins_multicore)) {
  now <- Sys.time()
  pm_score_list[[i]] <- bind_rows(parallel::mclapply(snn_bins_multicore[[i]], get_pm_score, mc.cores = ncores))
  setTxtProgressBar(pb, i)
  Sys.time() - now
}

pm_score_tbl <- bind_rows(pm_score_list)

write_tsv(pm_score_tbl, "/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/patient_mixing_scores.tsv")
pm_score_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/patient_mixing_scores.tsv")


