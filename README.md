# RNA-Sequence-Data-Analysis-uing-R
An RNA-seq analysis project investigating host transcriptional responses to respiratory virus infection (SARS-CoV-2) using R.

## Overview

RNA-Seq data analysis was performed using R, focusing on a real-world dataset published by Blanco-Melo et al. in Cell in 2020. In this project, normal human bronchial epithelial cells (NHBE), were infected with respiratory viruses, and the resulting transcriptomic profiles were compared to understand the distinct host responses triggered by the virus.

## Objectives

- Explore RNA-seq count data
- Identify differentially expressed genes
- Compare gene-expression patterns across experimental conditions
- Visualize biological patterns associated with viral infection
- Interpret the resulting host transcriptional responses

## Dataset

The analysis is based on data from ;

Accession: GSE147507

>Blanco-Melo, Daniel et al. “Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19.” Cell vol. 181,5 (2020): 1036-1045.e9. doi:10.1016/j.cell.2020.04.026

## Workflow

Analysis of a real-world RNA sequencing RNA-seq count data.
 Import and preprocess the dataset.
 Perform quality control.
 Normalize the expression values.
 Conduct differential expression analysis.
 Create visualization plots.


## Repository Structure

├── README.md
├── Data # Raw or processed count data
├── RNA Seq data analysis.R 
├This project explores transcriptional responses in  normal human bronchial── data/
└── results


## Requirements

- R 4.6.1
- RStudio
install BiocManager and gprofiler2 by using the command install.packages()
follow scripts <- RNA Seq data analysis.R

## Key Results

Blanco-Melo et al. highlighted in their publication that "SARS-CoV-2 infection induces low IFN-I and -III levels with a moderate ISG response" and report a strong chemokine expression. I simplified and normalized the data in such a way that I can perform the basic analyses also on a laptop with reduced CPU power. Throughout this exercise, I explored how the authors arrived at these conclusions.

## Visualization Examples
<img width="4527" height="3937" alt="sampleHeatmap_NHBE" src="https://github.com/user-attachments/assets/a211485e-c399-43de-aaa3-efde5c84c4d8" />
<img width="4330" height="3937" alt="PCAplot_NHBE" src="https://github.com/user-attachments/assets/67b5c28d-0553-4133-9807-bd1c97babdf5" />
<img width="480" height="480" alt="Volcano PLot_NHBE" src="https://github.com/user-attachments/assets/0d493166-9ca1-44d0-8da3-cdda0d5f7365" />
[gseaNHBEup.xlsx](https://github.com/user-attachments/files/35281695/gseaNHBEup.xlsx)
[gseaNHBEdown.xlsx](https://github.com/user-attachments/files/35281681/gseaNHBEdown.xlsx)
<img width="866" height="621" alt="BarPlot_Upregulated" src="https://github.com/user-attachments/assets/d4f93023-50e5-44e2-9015-6d03bb36f5c3" />
<img width="866" height="621" alt="BarPLot_downregulated" src="https://github.com/user-attachments/assets/591fae3c-7068-41d4-8928-70dd9b0a45a7" />


> sampleHeatmap_NHBE.png
> PCAplot_NHBE.png
> VolcanoPlot_NHBE.png
> BarPlot_downregulated.png
> Barplot_Upregulated.png


## Citation

Blanco-Melo, Daniel et al. “Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19.” Cell vol. 181,5 (2020): 1036-1045.e9. doi:10.1016/j.cell.2020.04.026
