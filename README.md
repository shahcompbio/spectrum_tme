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
- [x] v6 cell type objects (pipeline output, for annotated objects, see further below): `/work/shah/isabl_data_lake/analyses/68/75/6875/RNASCP/outs/{cell_type}.rds`
- [x] v6 cell type cluster markers: `/work/shah/isabl_data_lake/analyses/68/76/6875/RNASCP/outs/{cell_type}_markers.tsv`
- [x] SPECTRUM v7 pre-Rx cohort expression object: `/work/shah/isabl_data_lake/analyses/84/78/8478/RNASCP/outs/`
- [x] SPECTRUM v7 pre-Rx cohort embedding table: `/work/shah/isabl_data_lake/analyses/84/78/8478/RNASCP/outs/cells.tsv`
- [x] v7 cell type objects (pipeline output, for annotated objects, see further below): `/work/shah/isabl_data_lake/analyses/84/78/8478/RNASCP/outs/{cell_type}.rds`
- [x] v7 cell type cluster markers: `/work/shah/isabl_data_lake/analyses/84/78/8478/RNASCP/outs/{cell_type}_markers.tsv`

**Meta data** 

(File paths in this repo) 

- [x] scRNA-seq sample sheet: 
    - [x] [google spreadsheet](https://docs.google.com/spreadsheets/d/1plhIL1rH2IuQ8b_komjAUHKKrnYPNDyhvNNRsTv74u8/edit?ts=5d406b84#gid=1078838729) owned by Ignacio
    - [x] downloaded copy, Dec 3rd, 2020: [`_data/small/MSK\ SPECTRUM\ -\ Single\ cell\ RNA-seq_v6.xlsx`](_data/small/MSK\ SPECTRUM\ -\ Single\ cell\ RNA-seq_v6.xlsx)
- [x] color code: [`_data/small/signatures/hgsc_v7_colors.yaml`](_data/small/signatures/hgsc_v7_colors.yaml)
- [ ] mutational signature consensus labels: [`_data/small/mutational_signatures_summary_v3.tsv`](_data/small/mutational_signatures_summary_v3.tsv)

## Preprocessing scripts

**Cohort** 

- [ ] compute patient specificity from SNN graph
    - [ ] [`_preprocessing/010_SPECTRUM_freeze_v7_patient_mixing_multicore.R`](_preprocessing/010_SPECTRUM_freeze_v7_patient_mixing_multicore.R)
- [x] run consensusOV for TCGA subtype analysis 
    - [x] helper script to start runs: [`_preprocessing/020_SPECTRUM_freeze_v7_TCGA_signatures.R`](_preprocessing/020_SPECTRUM_freeze_v7_TCGA_signatures.R)
    - [x] parameterized markdown to compute consensusOV for a given `patient.rds` input file [`_preprocessing/021_SPECTRUM_freeze_v7_TCGA_signatures_per_patient.Rmd`](_preprocessing/021_SPECTRUM_freeze_v7_TCGA_signatures_per_patient.Rmd)
    - [x] summarized report of TCGA analysis [`_preprocessing/022_SPECTRUM_freeze_v7_TCGA_summary.html`](_preprocessing/022_SPECTRUM_freeze_v7_TCGA_summary.html)

- [ ] output files and objects:
    - [ ] patient specificity per cell
    - [x] TCGA consensusOV: `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/consensusOV/SPECTRUM_freeze_v7_consensusOV.tsv`

**Cell type subsets**

- [x] common tasks
    - [x] annotate sub type clusters and filter out doublet clusters
    - [x] run progeny for pathway scoring
    - [x] score gene signature modules

- [ ] individual scripts for common tasks:
    - [ ] B.super [`_preprocessing/031_SPECTRUM_freeze_v7_prepro_B.super.Rmd`](_preprocessing/031_SPECTRUM_freeze_v7_prepro_B.super.Rmd)
    - [x] Fibroblast.super [`_preprocessing/032_SPECTRUM_freeze_v7_prepro_Fibroblast.super.Rmd`](_preprocessing/032_SPECTRUM_freeze_v7_prepro_Fibroblast.super.Rmd)
    - [x] T.super [`_preprocessing/033_SPECTRUM_freeze_v7_prepro_T.super.Rmd`](_preprocessing/033_SPECTRUM_freeze_v7_prepro_T.super.Rmd)
    - [x] Myeloid.super [`_preprocessing/034_SPECTRUM_freeze_v7_prepro_Myeloid.super.Rmd`](_preprocessing/034_SPECTRUM_freeze_v7_prepro_Myeloid.super.Rmd)
    - [x] Ovarian.cancer.cell.super [`_preprocessing/035_SPECTRUM_freeze_v7_prepro_Ovarian.cancer.cell.super.Rmd`](_preprocessing/035_SPECTRUM_freeze_v7_prepro_Ovarian.cancer.cell.super.Rmd)
    - [ ] Endothelial.super [`_preprocessing/036_SPECTRUM_freeze_v7_prepro_Endothelial.super.Rmd`](_preprocessing/036_SPECTRUM_freeze_v7_prepro_Endothelial.super.Rmd)
    
- [x] cell type specific tasks
    - [x] Cancer cells: filter out low quality cells <1% mito reads and recluster 
    - [x] T cells: subcluster CD8 T cells and run diffusion map on CD8 subset
    - [x] Myeloid cells: run diffusion map on macrophage subset

## Output files and objects

- [ ] major cell types
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Ovarian.cancer.super_processed_filtered_annotated.rds`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/T.super_processed_filtered_annotated.rds`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Myeloid.super_processed_filtered_annotated.rds`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Fibroblast.super_processed_filtered_annotated.rds`
  - [ ] `B.super.annotated.rds`
  - [ ] `Endothelial.cell.super.annotated.rds`
- [x] subset objects 
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/CD8.T_processed_filtered.rds`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/DCs_processed.rds`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Macrophages_processed.rds`
- [ ] supplementary marker tables
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/Ovarian.cancer.super_marker_table_annotated.tsv`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/T.super_marker_table_annotated_full.tsv`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/Myeloid.super_marker_table_annotated.tsv`
  - [x] `/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/Fibroblast.super_marker_table_annotated.tsv`


## Plotting scripts

**Cohort (Fig 2)**

[`110_cohort_plotting_umap_v7.html`](110_cohort_plotting_umap_v7.html)

- [x] UMAP embeddings
- [x] Patient specificity
- [x] TCGA subtypes

[`115_cohort_plotting_composition_v7.html`](115_cohort_plotting_composition_v7.html)

- [x] Composition barplots scRNA
- [x] Composition barplots mpIF

**Cancer cell (Fig 3)**

[`120_cancer_cell_plotting_umap_v7.html`](120_cancer_cell_plotting_umap_v7.html)

- [x] UMAP embeddings
- [x] Marker heatmap
- [x] Composition barplots
- [x] Pathway heatmap and boxplots

**T cell (Fig 5)**

[`130_T_cell_plotting_umap_comp_heatmap_v7.html`](130_T_cell_plotting_umap_comp_heatmap_v7.html)

- [x] UMAP embeddings
- [x] Marker heatmap
- [x] Composition barplots

[`135_T_cell_CD8_diffusion_v7.html`](135_T_cell_CD8_diffusion_v7.html)

- [x] DC embeddings
- [x] Pseudotime trajectories
- [x] Module boxplots
- [x] Correlation scatter

**Myeloid cell (Fig 6)**

[`140_Myeloid_cell_plotting_umap_comp_heatmap_v7.html`](140_Myeloid_cell_plotting_umap_comp_heatmap_v7.html)

- [x] UMAP embeddings
- [x] Marker heatmap
- [x] Composition barplots

[`145_Macrophage_diffusion_v7.html`](145_Macrophage_diffusion_v7.html)

- [ ] UMAP marker gene embeddings
- [ ] Module score boxplots



