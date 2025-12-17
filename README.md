# τ-Axis: Immune–Metabolic Demand–Capacity Framework

This repository contains the analysis code and supporting resources for the **τ-Axis framework**, a systems biology approach that quantifies the relationship between immune activation and cellular energy-producing capacity.

The τ-Axis integrates pathway-level transcriptomic modules to distinguish **metabolically compensated inflammation** from **energy-constrained inflammatory states**, with application to **autism spectrum disorder (ASD)** and **acute sepsis**.

## Conceptual Overview

Biological systems face fluctuating immune and metabolic demands. While acute inflammatory states (e.g. sepsis) can upregulate glycolysis and substrate-level phosphorylation to meet energy demand, chronic immune activation may persist in the absence of adequate metabolic compensation.

The τ-Axis captures this balance by integrating:
- Immune signaling modules (e.g. TNFα/NF-κB, IL-6/STAT3, IFNγ, IL-10, IL-4/Th2)
- Metabolic capacity modules (e.g. glycolysis, cytosolic and mitochondrial substrate-level phosphorylation, FAO)
- Hypoxia and stress adaptation (HIF-1α–related signaling)

## Repository Contents

- `data/`  
  Curated and preprocessed transcriptomic datasets used in the analyses.

- `modules/`  
  Gene sets defining immune and metabolic pathway modules.

- `scripts/`  
  R scripts for module scoring, τ-Axis computation, statistical analysis, and visualization.

- `figures/`  
  Code used to generate publication-quality figures.

## Methods Summary

For each sample, pathway module scores are computed as the mean expression of member genes.  
The τ-Axis metric is defined as a weighted integration of metabolic and immune module scores, capturing **demand–capacity coupling or mismatch** at the systems level.

All analyses are fully reproducible using the provided scripts.

## Requirements

- R (≥ 4.2)
- tidyverse
- limma / edgeR
- ggplot2
- patchwork

Additional package requirements are specified within individual scripts.

## Reproducibility

All figures and tables reported in the associated manuscript can be regenerated from this repository. 
Data Availability and Ethics

This project exclusively analyzes publicly available, de-identified datasets. No new human or animal data were generated for this study.

The primary datasets used include:

GSE185263 – Whole-blood RNA-seq data from adults with community-acquired sepsis and healthy controls

Autism spectrum disorder cohort (GSE18123)

All datasets were obtained from the NCBI Gene Expression Omnibus (GEO) under their respective data-use policies. Original studies received appropriate ethical approval and informed consent as reported by the dataset authors.

This repository provides analysis code, gene module definitions, and visualization pipelines. 

## Citation

If you use this code or framework, please cite:

> Dervishi A. *et al.*  
> **The τ-Axis: Quantifying Immune–Metabolic Demand–Capacity Mismatch Across Disease States.**  
> (Manuscript under review)

## License

This repository is released for academic and non-commercial use.  
Please contact the author for reuse beyond these terms.
<img width="936" height="1290" alt="image" src="https://github.com/user-attachments/assets/2d979371-5809-4054-beaa-2cc6d2a6ed0d" />
