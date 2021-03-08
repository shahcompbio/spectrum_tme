## RUN EXAMPLE
## Rscript /home/uhlitzf/spectrum_tme/_preprocessing/020_SPECTRUM_TCGA_signatures_helper.R --seu_obj_path /work/shah/isabl_data_lake/analyses/07/38/738/results/integrated_harmony_seurat.rdata --patient OV025 

library(magrittr)
library(tidyverse)

## command line argument parsing
cargs <- list()

cargs_user <- commandArgs(trailingOnly = T) %>% 
  paste0(collapse = " ") %>% 
  str_split("--") %>% 
  lapply(str_split, " ") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) set_names(list(x[-1]), x[1])) %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) x[x!=""]) %>% 
  .[-1]

for (i in 1:length(cargs_user)) {
  cargs[[names(cargs_user)[i]]] <- cargs_user[[names(cargs_user)[i]]]
}

knit_root_dir <- rprojroot::find_rstudio_root_file()
setwd(knit_root_dir)

run_sub <- function(seu_obj_path = cargs$seu_obj_path, patient) {
  rmarkdown::render(
    "_preprocessing/021_SPECTRUM_freeze_v6_TCGA_signatures_per_patient.Rmd", params = list(
      seu_obj_path = seu_obj_path
    ),
    output_file = paste0("020_SPECTRUM_v6_freeze_TCGA_signatures_", patient, ".html"),
    output_dir = "_html",
    knit_root_dir = knit_root_dir
  )
}

#cat(paste0("run_sub('", foo, "', '", str_remove_all(basename(foo), ".rds"), "')"), sep = "\n")

run_sub('/work/shah/isabl_data_lake/analyses/16/59/1659/patient.rdata', "002")
run_sub('/work/shah/isabl_data_lake/analyses/17/02/1702/patient.rdata', "003")
run_sub('/work/shah/isabl_data_lake/analyses/16/88/1688/patient.rdata', "007")
run_sub('/work/shah/isabl_data_lake/analyses/18/69/1869/patient.rdata', "008")
run_sub('/work/shah/isabl_data_lake/analyses/18/70/1870/patient.rdata', "009")
run_sub('/work/shah/isabl_data_lake/analyses/17/11/1711/patient.rdata', "014")
run_sub('/work/shah/isabl_data_lake/analyses/17/29/1729/patient.rdata', "022")
run_sub('/work/shah/isabl_data_lake/analyses/19/11/1911/patient.rdata', "024")
run_sub('/work/shah/isabl_data_lake/analyses/17/30/1730/patient.rdata', "025")
run_sub('/work/shah/isabl_data_lake/analyses/19/10/1910/patient.rdata', "026")
run_sub('/work/shah/isabl_data_lake/analyses/18/77/1877/patient.rdata', "028")
run_sub('/work/shah/isabl_data_lake/analyses/16/51/1651/patient.rdata', "031")
run_sub('/work/shah/isabl_data_lake/analyses/18/99/1899/patient.rdata', "036")
run_sub('/work/shah/isabl_data_lake/analyses/19/13/1913/patient.rdata', "037")
run_sub('/work/shah/isabl_data_lake/analyses/18/53/1853/patient.rdata', "041")
run_sub('/work/shah/isabl_data_lake/analyses/19/00/1900/patient.rdata', "042")
run_sub('/work/shah/isabl_data_lake/analyses/19/12/1912/patient.rdata', "045")
run_sub('/work/shah/isabl_data_lake/analyses/19/01/1901/patient.rdata', "049")
run_sub('/work/shah/isabl_data_lake/analyses/18/74/1874/patient.rdata', "050")
run_sub('/work/shah/isabl_data_lake/analyses/19/06/1906/patient.rdata', "051")
run_sub('/work/shah/isabl_data_lake/analyses/19/03/1903/patient.rdata', "052")
run_sub('/work/shah/isabl_data_lake/analyses/19/05/1905/patient.rdata', "053")
run_sub('/work/shah/isabl_data_lake/analyses/19/02/1902/patient.rdata', "054")
run_sub('/work/shah/isabl_data_lake/analyses/18/72/1872/patient.rdata', "065")
run_sub('/work/shah/isabl_data_lake/analyses/18/73/1873/patient.rdata', "067")
run_sub('/work/shah/isabl_data_lake/analyses/19/04/1904/patient.rdata', "068")
run_sub('/work/shah/isabl_data_lake/analyses/18/71/1871/patient.rdata', "070")
run_sub('/work/shah/isabl_data_lake/analyses/19/08/1908/patient.rdata', "071")
run_sub('/work/shah/isabl_data_lake/analyses/19/07/1907/patient.rdata', "075")
run_sub('/work/shah/isabl_data_lake/analyses/19/09/1909/patient.rdata', "077")
run_sub('/work/shah/isabl_data_lake/analyses/46/17/4617/patient.rdata', "080")
run_sub('/work/shah/isabl_data_lake/analyses/46/18/4618/patient.rdata', "081")
run_sub('/work/shah/isabl_data_lake/analyses/46/20/4620/patient.rdata', "082")
run_sub('/work/shah/isabl_data_lake/analyses/46/15/4615/patient.rdata', "083")
run_sub('/work/shah/isabl_data_lake/analyses/45/84/4584/patient.rdata', "085")
run_sub('/work/shah/isabl_data_lake/analyses/46/16/4616/patient.rdata', "089")
run_sub('/work/shah/isabl_data_lake/analyses/45/98/4598/patient.rdata', "090")
run_sub('/work/shah/isabl_data_lake/analyses/67/84/6784/patient.rdata', "105")
run_sub('/work/shah/isabl_data_lake/analyses/68/06/6806/patient.rdata', "107")
run_sub('/work/shah/isabl_data_lake/analyses/68/10/6810/patient.rdata', "110")
run_sub('/work/shah/isabl_data_lake/analyses/68/09/6809/patient.rdata', "112")
run_sub('/work/shah/isabl_data_lake/analyses/71/96/7196/patient.rdata', "115")
run_sub('/work/shah/isabl_data_lake/analyses/79/57/7957/patient.rdata', "116")
run_sub('/work/shah/isabl_data_lake/analyses/84/76/8476/patient.rdata', "118")



