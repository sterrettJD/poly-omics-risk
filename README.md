# Poly-omics risk scores

## Project Description
Integrating multi-omics data to a polygenic risk score framework to model inflammatory bowel disease.

What started as an [Interdisciplinary Quantitative Biology](https://www.colorado.edu/certificate/iqbiology/) team rotation project in Spring 2022 has spiraled into a larger project and conversation around multi-omics modeling. 

This repository supports the following publication:

>Poly-omic risk scores predict inflammatory bowel disease diagnosis. CH Arehart,\* JD Sterrett,\* RL Garris, RE Quispe-Pilco, CR Gignoux, LM Evans,\^ MA Stanislawski.\^ *mSystems* in press (2023). doi: 10.1128/msystems.00677-23.

*: join first authors    
^: joint senior authors

## Repository Contents
This repository contains directories for each -omics data type, with a script in each for quality control, transformation, feature selection, model evaluation, and assessing variable importance.

Additionally, `AMON` contains scripts used for integrating the metagenomics data with metabolomics data to assess potential microbial sources of LASSO-selected compounds.

Scripts integrating risk scores and feature importance from multiple -omics datatypes are found in this main directory.

`train_test_split.R` contains steps for the train-test split and plots supporting the rationalization for performing it as such.

## Data Source
- Human Microbiome Project 2
  - [Data](https://ibdmdb.org/tunnel/public/summary.html)
  - [Outcome manuscript](https://www.nature.com/articles/s41586-019-1237-9)
  - Data types
    - Metagenomics
    - Viromics
    - Metabolomics
    - Metatranscriptomics