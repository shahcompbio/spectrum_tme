# SPECTRUM TME scRNA figure repository

The repository includes all analysis necessary to arrive at the figures for the SPECTRUM TME paper. 

Scripts are separated into preprocessing scripts and plotting scripts. 

## Data dependencies

**Processed expression data** 

The output from Nick Ceglia's scRNA-seq pipeline serves as an input to the analysis. 
- Cohort expression object: 
- Cohort embedding table: 
- Cell type objects: 
- Cell type object cluster markers: 

**Meta data** 

- scRNA-seq sample sheet (google spreadsheet)[https://docs.google.com/spreadsheets/d/1plhIL1rH2IuQ8b_komjAUHKKrnYPNDyhvNNRsTv74u8/edit?ts=5d406b84#gid=1078838729]
- color code
- mutational signatures

## Preprocessing scripts

**Cohort** 

- run consensusOV for TCGA subtype analysis
- compute patient specificity from SNN graph

**Cell type subsets**

- all cell types:
    - annotate sub type clusters and filter out doublet clusters
    - run progeny for pathway scoring
    - score gene signature modules
- Cancer cells: filter out low quality cells <1% mito reads
- T cells: subcluster CD8 T cells and run diffusion map ind CD8 subset
- Myeloid cells: run diffusion map on macrophage subset

## Plotting scripts

**Cohort level**

- UMAP embeddings
- Patient specificity
- TCGA subtypes
- Composition barplots

**Cancer cell**

- UMAP embeddings
- Composition barplots
- Marker heatmap

**T cell**

- UMAP embeddings
- Composition barplots
- Marker heatmap
- Diffusion maps
- Module trajectories

**Myeloid cell**

- UMAP embeddings
- Composition barplots
- Marker heatmap
- Diffusion maps
- Module trajectories





