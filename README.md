# SPECTRUM TME scRNA figure repository

The repository includes all analysis necessary to arrive at the figures for the SPECTRUM TME paper. 

Scripts are separated into preprocessing scripts and plotting scripts. 

**STATUS**

- [ ] OPEN
- [x] DONE

## Data dependencies

**Processed expression data** 

(File paths on juno) 

The output from Nick Ceglia's scRNA-seq pipeline serves as an input to the analysis. 

- [x] SPECTRUM v6 pre-Rx cohort expression object: `/work/shah/isabl_data_lake/analyses/68/75/6875/RNASCP/outs/`
- [x] SPECTRUM v6 pre-Rx cohort embedding table: `/work/shah/isabl_data_lake/analyses/68/75/6875/RNASCP/outs/cells.tsv`
- [x] v6 cell type objects: `/work/shah/isabl_data_lake/analyses/68/75/6875/RNASCP/outs/{cell_type}.rds`
- [x] v6 cell type cluster markers: `/work/shah/isabl_data_lake/analyses/68/76/6875/RNASCP/outs/{cell_type}_markers.tsv`

**Meta data** 

(File paths in this repo) 

- [x] scRNA-seq sample sheet: 
    - [x] [google spreadsheet](https://docs.google.com/spreadsheets/d/1plhIL1rH2IuQ8b_komjAUHKKrnYPNDyhvNNRsTv74u8/edit?ts=5d406b84#gid=1078838729) owned by Ignacio
    - [x] downloaded copy, Dec 3rd, 2020: [`_data/small/MSK\ SPECTRUM\ -\ Single\ cell\ RNA-seq_v6.xlsx`](_data/small/MSK\ SPECTRUM\ -\ Single\ cell\ RNA-seq_v6.xlsx)
- [x] color code: [`_data/small/signatures/hgsc_v6_colors.yaml`](_data/small/signatures/hgsc_v6_colors.yaml)
- [ ] mutational signature consensus labels: [`_data/small/mutational_signatures_summary_v3.tsv`](_data/small/mutational_signatures_summary_v3.tsv)

## Preprocessing scripts

**Cohort** 

- [ ] compute patient specificity from SNN graph
    - [ ] [`_preprocessing/010_SPECTRUM_freeze_v6_patient_mixing_multicore.R`](_preprocessing/010_SPECTRUM_freeze_v6_patient_mixing_multicore.R)
- [ ] run consensusOV for TCGA subtype analysis 
    - [x] helper script to start runs: [`_preprocessing/020_SPECTRUM_freeze_v6_TCGA_signatures.R`](_preprocessing/020_SPECTRUM_freeze_v6_TCGA_signatures.R)
    - [x] parameterized markdown to compute consensusOV for a given `patient.rds` input file [`_preprocessing/021_SPECTRUM_freeze_v6_TCGA_signatures_per_patient.Rmd`](_preprocessing/021_SPECTRUM_freeze_v6_TCGA_signatures_per_patient.Rmd)
    - [ ] summarized report of TCGA analysis [`_preprocessing/022_SPECTRUM_freeze_v6_TCGA_summary.Rmd`](_preprocessing/022_SPECTRUM_freeze_v6_TCGA_summary.Rmd)

- [ ] output files and objects:
    - [ ] patient specificity per cell
    - [x] TCGA consensusOV: `/work/shah/uhlitzf/data/SPECTRUM/freeze/v6/consensusOV/SPECTRUM_freeze_v6_consensusOV.tsv`

**Cell type subsets**

- [ ] common tasks
    - [ ] annotate sub type clusters and filter out doublet clusters
    - [ ] run progeny for pathway scoring
    - [ ] score gene signature modules

- [ ] individual scripts for common tasks:
    - [ ] B.super [`_preprocessing/031_SPECTRUM_freeze_v6_prepro_B.super.Rmd`](_preprocessing/031_SPECTRUM_freeze_v6_prepro_B.super.Rmd)
    - [ ] Fibroblast.super [`_preprocessing/032_SPECTRUM_freeze_v6_prepro_Fibroblast.super.Rmd`](_preprocessing/032_SPECTRUM_freeze_v6_prepro_Fibroblast.super.Rmd)
    - [ ] T.super [`_preprocessing/033_SPECTRUM_freeze_v6_prepro_T.super.Rmd`](_preprocessing/033_SPECTRUM_freeze_v6_prepro_T.super.Rmd)
    - [ ] Myeloid.super [`_preprocessing/034_SPECTRUM_freeze_v6_prepro_Myeloid.super.Rmd`](_preprocessing/034_SPECTRUM_freeze_v6_prepro_Myeloid.super.Rmd)
    - [x] Ovarian.cancer.cell.super [`_preprocessing/035_SPECTRUM_freeze_v6_prepro_Ovarian.cancer.cell.super.Rmd`](_preprocessing/035_SPECTRUM_freeze_v6_prepro_Ovarian.cancer.cell.super.Rmd)
    
- [ ] cell type specific tasks
    - [x] Cancer cells: filter out low quality cells <1% mito reads and recluster []()
    - [ ] T cells: subcluster CD8 T cells and run diffusion map on CD8 subset
    - [ ] Myeloid cells: run diffusion map on macrophage subset

- [ ] output files and objects:
    - [ ] `B.super.annotated.rds`
    - [ ] `Fibroblast.super.annotated.rds`
    - [ ] `T.super.annotated.rds`
    - [ ] `Myeloid.super.annotated.rds`
    - [x] `Ovarian.cancer.cell.super.annotated.rds`
    - [ ] `CD8.T.cell.annotated.rds`
    - [ ] `Macrophage.annotated.rds`


## Plotting scripts

**Cohort level**

[`110_cohort_plotting_umap.Rmd`](110_cohort_plotting_umap.Rmd)

- [x] UMAP embeddings
- [ ] Patient specificity
- [ ] TCGA subtypes

[`115_cohort_plotting_composition.Rmd`](110_cohort_plotting_composition.Rmd)

- [ ] Composition barplots

**Cancer cell**

[`120_cancer_cell_plotting.Rmd`](120_cancer_cell_plotting.Rmd)

- [x] UMAP embeddings
- [ ] Composition barplots

**T cell**

[`130_T_cell_plotting.Rmd`](130_T_cell_plotting.Rmd)

- [ ] UMAP embeddings
- [ ] Composition barplots
- [ ] Marker heatmap

**Myeloid cell**

[`140_Myeloid_cell_plotting.Rmd`](140_Myeloid_cell_plotting.Rmd)

- [ ] UMAP embeddings
- [ ] Composition barplots
- [ ] Marker heatmap

**CD8 T cell**

[`150_CD8_T_cell_plotting.Rmd`](150_CD8_T_cell_plotting.Rmd)

- [ ] Diffusion maps
- [ ] Module trajectories

**Macrophages**

[`160_macrophage_plotting.Rmd`](160_macrophage_plotting.Rmd)

- [ ] Diffusion maps
- [ ] Module trajectories



