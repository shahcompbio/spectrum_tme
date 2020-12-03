library(magrittr)
library(tidyverse)
library(readxl)
library(Seurat)

setwd(rstudioapi::getActiveProject())
# setwd("/home/uhlitzf/spectrum_tme")
ncores <- 40

meta_tbl <- read_excel("_data/small/MSK SPECTRUM - Single cell RNA-seq_v6.xlsx", sheet = 2) %>% 
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-")) %>% 
  mutate(sample = isabl_id)

# seu_merge <- read_rds("/work/shah/isabl_data_lake/analyses/16/52/1652/cohort_merged.rdata")
# seu_snn <- seu_merge@graphs$RNA_snn
# write_rds(seu_snn, "/work/shah/uhlitzf/data/SPECTRUM/freeze/v5/cohort_snn_matrix.rds")
seu_snn <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v5/cohort_snn_matrix.rds")
seu_tbl <- read_tsv("/work/shah/isabl_data_lake/analyses/16/52/1652/cells.tsv") %>% 
  left_join(meta_tbl, by = "sample") %>% 
  filter(therapy == "pre-Rx") %>% 
  mutate(cell_type = str_replace_all(cell_type, "\\.", " "))

seu_snn <- seu_snn[seu_tbl$cell_id, seu_tbl$cell_id]
rownames(seu_snn) <- str_sub(rownames(seu_snn), 13, 15)

get_bins <- function(x, y) split(x, ceiling(seq_along(x)/y))
snn_bins <- get_bins(1:ncol(seu_snn), 20000)
snn_bins_multicore <- lapply(snn_bins, get_bins, 20000/40)

pid_freq <- select(seu_tbl, cell_id, cell_type, pid1) %>% 
  group_by(cell_type, pid1) %>% 
  tally %>% 
  group_by(cell_type) %>% 
  mutate(nrel = n/sum(n)) %>% 
  ungroup %>% 
  select(-n)

get_neighbours <- function(i) {
  seu_snn[,i] %>% 
    .[.!=0] %>% 
    sort(decreasing = T) %>% 
    enframe("pid2", "value")
}

get_pm_score <- function(cell_index_range) {
  
  snn_tbl <- lapply(cell_index_range, get_neighbours) %>% 
    setNames(colnames(seu_snn)[cell_index_range]) %>% 
    bind_rows(.id = "cell_id") %>% 
    mutate(pid1 = str_sub(cell_id, 13, 15))
  
  pm_score <- snn_tbl %>% 
    left_join(select(seu_tbl, cell_id, cell_type), by = "cell_id") %>% 
    left_join(pid_freq, by = c("pid1", "cell_type")) %>% 
    group_by(cell_id) %>% 
    # slice(1:30) %>% 
    mutate(nexp = nrel * length(cell_id),
           nobs = sum(pid1 == pid2),
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

write_tsv(pm_score_tbl, "/work/shah/uhlitzf/data/SPECTRUM/freeze/v5/patient_mixing_scores.tsv")


