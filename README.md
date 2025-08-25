# Aging transcriptome is Associated with Poor Outcome Post-Surgical Removal of Early-Stage Lung Cancer

This repository contains the analysis code used in our study of transcriptomic differences between younger and older patients with stage I lung adenocarcinoma (LUAD). 
The project focuses on lower-airway gene expression programs associated with age and prognosis after surgical resection.

## Contents
- **/scripts** – R scripts for data preprocessing, differential expression, co-expression network analysis (MEGENA), and pathway enrichment.
- **/figures** – Code for generating key plots used in the manuscript.
- **/data** – Metadata 

## Key Analyses
- Differential gene expression between younger vs. older patients and between recurrence vs. no recurrence among each age group. 
- Identification of recurrence-associated transcriptional modules  
- Pathway annotation and cross-cohort consistency checks (NYU, TCGA, TRACERx)

## Requirements
- R (≥4.2) with packages: `DESeq2`, `edgeR`, `MEGENA`, `ggplot2`, `ComplexHeatmap`
