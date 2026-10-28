# Single-Cell RNA-seq Analysis of Lung Tissue

This project uses R and Seurat to analyse single-cell RNA sequencing data from the NML1 lung sample in the GSE132771 dataset.

The analysis explores gene expression patterns, assesses cell quality, identifies cell clusters, and examines marker genes to develop provisional cell-type annotations.

## Dataset

- **Dataset accession:** GSE132771
- **Sample analysed:** NML1
- **Data type:** Single-cell RNA sequencing (scRNA-seq)
- **Input format:** 10x Genomics gene-expression matrix
- **Source:** [NCBI Gene Expression Omnibus (GEO), GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771)


The analysis uses the raw gene-expression matrix for the NML1 sample. The dataset's original publication and metadata should be consulted when interpreting the biological findings.

## Software and Packages

The analysis was performed in R using RStudio.

Required R packages:

- `Seurat` for single-cell RNA-seq analysis
- `Matrix` for handling sparse gene-expression matrices
- `ggplot2` for visualisation
- `dplyr` for data manipulation

### Software Environment

The analysis was run using the following environment:

- **R:** 4.6.1
- **Operating system:** Windows 10 (64-bit)
- **Seurat:** 5.6.0
- **SeuratObject:** 5.4.0
- **Matrix:** 1.7-5
- **ggplot2:** 4.0.3
- **dplyr:** 1.2.1

Package versions are reported from the R session used for this analysis.


## Project Structure

```text
Single Cell RNA/
├── Single Cell RNA.Rproj
├── single_cell_analysis.R
├── README.md
├── GSE132771_RAW/
│   └── NML1/
├── seu_nml1_filtered.rds
└── figures/
    ├── quality_control_violin_plots.png
    ├── cluster_qc_violin.png
    ├── umap_clusters.png
    ├── selected_marker_dotplot.png
    ├── cluster_annotations.csv
    ├── significant_cluster_markers.csv
    ├── top10_markers_per_cluster.csv
    └── marker_expression_by_cluster.csv
```

The project structure above describes the principal analysis files and outputs. Additional files may be present in the project directory.

## Analysis Steps

The analysis is organised into eight main steps:

1. **Load and inspect the data:** Import the raw gene-expression matrix and examine its dimensions.
2. **Create a Seurat object:** Store gene counts and cell-level metadata in a Seurat object.
3. **Quality control:** Examine total counts, detected genes, and mitochondrial gene percentages.
4. **Filter cells:** Apply initial quality-control thresholds to remove cells that do not meet the selected criteria.
5. **Normalisation and variable features:** Normalise gene expression and identify highly variable genes.
6. **Dimensionality reduction and clustering:** Scale the data, perform principal component analysis (PCA), construct a neighbourhood graph, identify clusters, and visualise them using UMAP.
7. **Marker analysis and provisional annotation:** Identify cluster marker genes, compare selected marker expression patterns, and assign provisional cell-type labels.
8. **Documentation and reproducibility:** Document the methods, results, limitations, and steps required to reproduce the analysis.

## Quality Control and Filtering

The initial dataset contained 3,669 cells.

The following initial filtering criteria were applied:

- At least 200 detected genes per cell (`nFeature_RNA >= 200`).
- Mitochondrial gene percentage below 10% (`percent.mt < 10`).

After filtering, 3,627 cells remained. A total of 42 cells were excluded.

These thresholds were used for this exploratory analysis and have not been established as universally optimal. Alternative thresholds and additional quality-control methods may be evaluated in future analyses.

## Normalisation, Dimensionality Reduction and Clustering

Gene expression was normalised using Seurat's `LogNormalize` method with a scale factor of 10,000.

The analysis identified 2,000 highly variable genes using the `vst` selection method.

Principal component analysis (PCA) was performed on the variable genes. The first 20 principal components were used for neighbourhood graph construction, clustering, and UMAP visualisation.

Clustering was performed using the Louvain algorithm with a resolution of 0.5, producing 18 clusters in this analysis.

UMAP was used to visualise the relationships between cells in two dimensions. The resulting clusters are computational groupings and should not be treated as definitive cell identities without further validation.

## Marker Gene Analysis and Cell-Type Annotation

Cluster marker genes were identified using Seurat's `FindAllMarkers()` function with the following settings:

- `only.pos = TRUE`
- `min.pct = 0.25`
- `logfc.threshold = 0.25`

Markers were selected based on their expression differences between each cluster and the remaining cells. Significant marker results were retained using an adjusted p-value threshold of 0.05.

Selected marker genes were visualised using a dot plot, and their expression patterns were examined to support provisional cell-type annotations.

The current annotations include epithelial, endothelial, stromal, smooth muscle, myeloid, mast, plasma, and proliferating cell populations.

**Unresolved cluster:** Cluster 13 expresses both endothelial-associated and macrophage-associated markers. Its identity remains uncertain and requires further investigation.

The annotations are preliminary interpretations of one sample and should be validated using additional marker genes, appropriate reference datasets, and, where possible, independent samples.

## Results and Output Files

The `figures/` directory contains the principal visualisations and tables generated during the analysis:

- `quality_control_violin_plots.png`: Initial cell-level quality-control distributions.
- `cluster_qc_violin.png`: Quality-control metrics across clusters.
- `umap_clusters.png`: UMAP visualisation of the 18 clusters.
- `selected_marker_dotplot.png`: Expression patterns of selected marker genes across clusters.
- `cluster_annotations.csv`: Provisional cell-type annotations.
- `significant_cluster_markers.csv`: Marker genes meeting the selected significance criteria.
- `top10_markers_per_cluster.csv`: Selected top marker genes for each cluster.
- `marker_expression_by_cluster.csv`: Average expression of selected marker genes by cluster.

## Limitations

This analysis examines a single sample, NML1. The results describe patterns in that sample and cannot establish how consistently the observed populations occur across patients.

Cell-type annotations are provisional. Additional validation, doublet assessment, contamination checks, and comparison with suitable reference datasets may be needed.

The clustering resolution, number of principal components, and filtering thresholds are analytical choices that can affect the results.

Marker-gene significance at the cell level does not, by itself, establish patient-level statistical significance or biological generalisability.

## Reproducibility

The R script `single_cell_analysis.R` contains the analysis code. The project also includes a saved Seurat object, `seu_nml1_filtered.rds`, for resuming the analysis without repeating the initial steps.

To reproduce the analysis from the beginning, use the original dataset files and run the script in the intended order from the project directory. Required packages must be installed, and the expected input-file structure must be preserved.

## How to Run the Analysis

1. Install R and RStudio.
2. Install the required R packages if they are not already installed.
3. Open `Single Cell RNA.Rproj` in RStudio.
4. Confirm that the `GSE132771_RAW/NML1/` directory contains the required input files.
5. Open `single_cell_analysis.R`.
6. Run the script sequentially from the beginning.

The script uses relative file paths, so run it from the project directory. The original dataset files must be available for a complete run.

## Data Availability

The dataset is publicly available through the [NCBI Gene Expression Omnibus (GEO), accession GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771).

## Project Status

This is an exploratory single-sample analysis. Cluster annotations are provisional, and additional validation is required before drawing firm biological conclusions.