**RNAseK scRNA-seq Analysis Pipeline**

This repository contains an R-based workflow to analyze single-cell RNA sequencing data, focusing on the role of RNAseK in EMT trajectories and melanocyte differentiation.

**Key Features**

Data Integration: Combines WT and KO datasets using the Seurat pipeline.
Cell Annotation: Identifies populations including iPSCs, Migratory NCCs, and Melanoblasts.
Trajectory Inference: Maps differentiation paths and pseudotime using Monocle3.
Iterative DE: Employs a 25-iteration bootstrapping method to identify high-confidence differentially expressed genes.
EMT and Stemness Scoring: Calculates E/M ratios and transcriptomic diversity via CytoTRACE.
Functional Enrichment: Performs GO term analysis (BP, MF, CC) using clusterProfiler.

**Requirements**

The pipeline requires the following R libraries:

Seurat
monocle3
tidyverse
clusterProfiler
CytoTRACE
patchwork

**Workflow Summary**

QC and Filtering: Standardizes data based on mitochondrial gene percentage and feature density.
Normalization and Integration: Harmonizes Control and Knockout conditions.
Clustering and UMAP: Visualizes cell populations and condition-specific shifts.
Trajectory Analysis: Models the transition from epithelial progenitors to neural crest lineages.
Statistical Validation: Uses iterative sampling to ensure gene marker consistency.
Visualization: Generates publication-ready volcano plots, heatmaps, and EMT shift arrow plots.

**Biological Context**

This study investigates how the loss of RNAseK affects the progression from iPSCs to neural crest-derived melanocytes, specifically identifying transcriptomic stalling or the emergence of aberrant intermediate states.

Author: Logeshwar V R
Institution: CSIR–Institute of Genomics and Integrative Biology (IGIB)
