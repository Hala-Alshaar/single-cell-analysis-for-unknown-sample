# Single-Cell RNA-seq Analysis Pipeline

## Overview

This repository contains a complete single-cell RNA sequencing (scRNA-seq) analysis workflow implemented in **R** using the **Seurat** ecosystem. The pipeline performs quality control, normalization, dimensionality reduction, clustering, cell cycle analysis, cell type annotation, and tissue-of-origin inference for human scRNA-seq datasets.

The workflow is designed to provide an end-to-end analysis from raw count matrices to biologically interpretable cell populations.

---

## Features

* Quality control and filtering of low-quality cells
* Mitochondrial content assessment
* Data normalization and scaling
* Identification of highly variable genes
* Cell cycle scoring (G1, S, G2/M phases)
* Principal Component Analysis (PCA)
* UMAP visualization
* Graph-based clustering
* Cluster marker gene identification
* Cell type annotation using SingleR and Human Primary Cell Atlas
* Additional annotation using CellMarkerAccordion
* Cell composition analysis
* Tissue-of-origin inference based on cell-type distribution and marker genes

---

## Workflow

### 1. Data Loading

Create a Seurat object from raw count matrices and verify species-specific gene annotations.

### 2. Quality Control

Calculate mitochondrial gene percentages and remove cells with:

* Low gene counts
* Excessively high gene counts
* High mitochondrial content

### 3. Data Preprocessing

* Normalize expression values
* Identify highly variable genes
* Scale expression data

### 4. Cell Cycle Analysis

Assign cells to G1, S, or G2/M phases using canonical cell-cycle markers.

### 5. Dimensionality Reduction

* Perform PCA
* Select informative principal components
* Generate UMAP embeddings

### 6. Clustering

Identify transcriptionally similar cell populations using graph-based clustering.

### 7. Marker Gene Identification

Detect cluster-specific marker genes and visualize them using heatmaps.

### 8. Cell Type Annotation

Annotate clusters using:

* SingleR
* Human Primary Cell Atlas reference
* CellMarkerAccordion database

### 9. Biological Interpretation

Evaluate cell-type composition and infer the likely tissue of origin using canonical marker genes and annotation results.

---

## Required Packages

```r
Seurat
SingleR
celldex
SingleCellExperiment
dplyr
ggplot2
patchwork
cellmarkeraccordion
```

---

## Input Data

The pipeline expects a raw count matrix with:

* Genes as rows
* Cells as columns

Example:

```r
load("YKG.rds")
```

---

## Output

The pipeline generates:

* Filtered Seurat object
* UMAP visualizations
* Cluster assignments
* Marker gene tables
* Cell-type annotations
* Annotation comparison tables
* Tissue-of-origin interpretation

---

## Biological Interpretation

The workflow is particularly useful for identifying immune and stromal populations in human tissues. By combining reference-based annotation with marker-gene validation, it provides robust characterization of cellular heterogeneity within single-cell datasets.

---

## Citation

If you use this workflow in your research, please cite:

* Seurat
* SingleR
* Human Primary Cell Atlas
* CellMarkerAccordion

and the corresponding publications associated with these tools.

