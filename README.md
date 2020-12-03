# SPECTRUM TME scRNA figure repository

The repository includes all analysis necessary to arrive at the figures for the SPECTRUM TME paper. 

Scripts are separated into preprocessing scripts and plotting scripts. 

**STATUS**

- [_] OPEN
- [x] DONE

## Data dependencies

**Processed expression data** 

(File paths on juno) 

The output from Nick Ceglia's scRNA-seq pipeline serves as an input to the analysis. 

- [_] Cohort expression object: `/work/shah/isabl_data_lake/analyses/68/76/6876/RNASCP/outs/`
- [_] Cohort embedding table: `/work/shah/isabl_data_lake/analyses/68/76/6876/RNASCP/outs/cells.tsv`
- [_] Cell type objects: `/work/shah/isabl_data_lake/analyses/68/76/6876/RNASCP/outs/{}.rds`
- [_] Cell type object cluster markers: `/work/shah/isabl_data_lake/analyses/68/76/6876/RNASCP/outs/{}_markers.tsv`

**Meta data** 

(File paths in this repo) 

- [x] scRNA-seq sample sheet: 
    - [x] [google spreadsheet](https://docs.google.com/spreadsheets/d/1plhIL1rH2IuQ8b_komjAUHKKrnYPNDyhvNNRsTv74u8/edit?ts=5d406b84#gid=1078838729) owned by Ignacio
    - [x] downloaded copy, Dec 3rd, 2020: [`_data/small/MSK SPECTRUM - Single cell RNA-seq_v6.xlsx`](_data/small/MSK SPECTRUM - Single cell RNA-seq_v6.xlsx)
- [x] color code: [`_data/small/signatures/hgsc_v6_colors.yaml`](_data/small/signatures/hgsc_v6_colors.yaml)
- [x] mutational signature consensus labels: 

## Preprocessing scripts

**Cohort** 

- [_] compute patient specificity from SNN graph
    - [_] [`_preprocessing/010_SPECTRUM_freeze_v6_patient_mixing_multicore.R`](_preprocessing/010_SPECTRUM_freeze_v6_patient_mixing_multicore.R)
- [_] run consensusOV for TCGA subtype analysis 
    - [x] helper script to start runs: [`_preprocessing/020_SPECTRUM_freeze_v6_TCGA_signatures.R`](_preprocessing/020_SPECTRUM_freeze_v6_TCGA_signatures.R)
    - [x] parameterized markdown to compute consensusOV for a given `patient.rds` input file [`_preprocessing/021_SPECTRUM_freeze_v6_TCGA_signatures_per_patient.Rmd`](_preprocessing/021_SPECTRUM_freeze_v6_TCGA_signatures_per_patient.Rmd)
    - [_] summarized report of TCGA analysis [`_preprocessing/022_SPECTRUM_freeze_v6_TCGA_summary.Rmd`](_preprocessing/022_SPECTRUM_freeze_v6_TCGA_summary.Rmd)

- [_] output files and objects:
    - [_] consensusOV: `/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/consensusOV/SPECTRUM_freeze_v6_consensusOV.tsv`

**Cell type subsets**

- [_] common tasks
    - [_] annotate sub type clusters and filter out doublet clusters
    - [_] run progeny for pathway scoring
    - [_] score gene signature modules

- [_] individual scripts for common tasks:
    - [_] B.super [`_preprocessing/031_SPECTRUM_freeze_v6_prepro_B.super.Rmd`](_preprocessing/031_SPECTRUM_freeze_v6_prepro_B.super.Rmd)
    - [_] Fibroblast.super [`_preprocessing/032_SPECTRUM_freeze_v6_prepro_Fibroblast.super.Rmd`](_preprocessing/032_SPECTRUM_freeze_v6_prepro_Fibroblast.super.Rmd)
    - [_] T.super [`_preprocessing/033_SPECTRUM_freeze_v6_prepro_T.super.Rmd`](_preprocessing/033_SPECTRUM_freeze_v6_prepro_T.super.Rmd)
    - [_] Myeloid.super [`_preprocessing/034_SPECTRUM_freeze_v6_prepro_Myeloid.super.Rmd`](_preprocessing/034_SPECTRUM_freeze_v6_prepro_Myeloid.super.Rmd)
    - [_] Ovarian.cancer.cell.super [`_preprocessing/035_SPECTRUM_freeze_v6_prepro_Ovarian.cancer.cell.super.Rmd`](_preprocessing/035_SPECTRUM_freeze_v6_prepro_Ovarian.cancer.cell.super.Rmd)
    
- [_] cell type specific tasks
    - [_] Cancer cells: filter out low quality cells <1% mito reads
    - [_] T cells: subcluster CD8 T cells and run diffusion map on CD8 subset
    - [_] Myeloid cells: run diffusion map on macrophage subset

- [_] output files and objects:
    - [_] `B.super.annotated.rds`
    - [_] `Fibroblast.super.annotated.rds`
    - [_] `T.super.annotated.rds`
    - [_] `Myeloid.super.annotated.rds`
    - [_] `Ovarian.cancer.cell.super.annotated.rds`
    - [_] `CD8.T.cell.annotated.rds`
    - [_] `Macrophage.annotated.rds`


## Plotting scripts

**Cohort level**

[`110_cohort_plotting.Rmd`](110_cohort_plotting.Rmd)

- [_] UMAP embeddings
- [_] Patient specificity
- [_] TCGA subtypes
- [_] Composition barplots

**Cancer cell**

[`120_cancer_cell_plotting.Rmd`](120_cancer_cell_plotting.Rmd)

- [_] UMAP embeddings
- [_] Composition barplots
- [_] Marker heatmap

**T cell**

[`130_T_cell_plotting.Rmd`](130_T_cell_plotting.Rmd)

- [_] UMAP embeddings
- [_] Composition barplots
- [_] Marker heatmap

**Myeloid cell**

[`140_Myeloid_cell_plotting.Rmd`](140_Myeloid_cell_plotting.Rmd)

- [_] UMAP embeddings
- [_] Composition barplots
- [_] Marker heatmap

**CD8 T cell**

[`150_CD8_T_cell_plotting.Rmd`](150_CD8_T_cell_plotting.Rmd)

- [_] Diffusion maps
- [_] Module trajectories

**Macrophages**

[`160_macrophage_plotting.Rmd`](160_macrophage_plotting.Rmd)

- [_] Diffusion maps
- [_] Module trajectories



