# TSG-twohit

## Overview

This repository contains scripts used for analyses and figure generation in the study:

**Genomic and Evolutionary Determinants of Two-hit Frequencies in Tumor Suppressor Genes**

Each Python script is named after the subplot it generates, enabling direct mapping between code and manuscript figures.


## Repository Structure

├── scripts/ # Python and R scripts for analysis and plotting
<br>├── environment.yml # Python environment
<br>├── renv.lock # R environment

## Data Availability

All input data required to run the scripts can be downloaded from:
https://drive.ncbs.res.in/index.php/s/3Byt2dpsDoYmfaB

After downloading, place the files in a directory named 'data'.

## Environment Setup

### Python (Conda)

```bash
conda env create -f environment.yml
conda activate tsgtwohit
```
### R

```R
install.packages("renv")
renv::restore()
```

## Citation

Genomic and Evolutionary Determinants of Two-hit Frequencies in Tumor Suppressor Genes
Nivedita Mukherjee, Radhakrishnan Sabarinathan
bioRxiv 2026.02.17.706175; doi: https://doi.org/10.64898/2026.02.17.706175

