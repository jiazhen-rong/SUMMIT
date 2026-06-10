# SUMMIT Overview
This GitHub contains the beta version of the Spatial Mitochondrial Variant Test - **SUMMIT** package.

All the analysis codes are within the [Analysis](https://github.com/jiazhen-rong/SUMMIT/tree/main/Analysis) folder. 

All the results in our paper was based on SUMMIT v1.0.0. 

For Single-cell enriched mitochondrial variant related test and filtering, please view our previous package [scmtVT](https://github.com/jiazhen-rong/scmtVT).

## System Requirements
### Hardware requirements
```SUMMIT``` package requires a standard computer with enough RAM.

### Software requirements
#### OS Requirements
The software was tested on:

 - MacOS Sequoia (15.6.1)
 - Rocky Linux 8.10 (Green Obsidian)
   
#### R Dependencies
```SUMMIT``` depends on the following R packages:

``` r
Matrix, ggplot2, dplyr, tidyr, cowplot, distances, igraph, biglm,
pheatmap, ComplexHeatmap, circlize 
```

## Installation
The installation shall take < 5 minutes.

To install the package, please use the following commands:
``` r
install.packages("devtools")
devtools::install_github("jiazhen-rong/SUMMIT") # install
library(SUMMIT) # load
```
or directly copy from git:
``` linux
git clone https://github.com/jiazhen-rong/SUMMIT.git
```

## Tutorial

 - [SUMMIT on Barrett's Esophagus Sample with Localized Clones](https://jiazhen-rong.github.io/SUMMIT/Tutorial_a200_s6.html) (~10 minutes)
 - [Mixture Differentially Expressed Genes Analysis on Barrett's Esophagus Sample with Localized Clones](https://jiazhen-rong.github.io/SUMMIT/Tutorial_a200_s6_mixtureDEG.html) (~2 hours)

## Citations
If you used the package in your research, please cite:

***Mitochondrial clone tracing within spatially intact human tissues** <br/>
Sydney A Bracht^*, Jiazhen Rong^*, Rodrigo A Gier, Maureen DeMarshall, Hailey Golden, Diya Dhakal, Jayne C McDevitt, Feiyan Mo, Emma E Furth, Alexandra Strauss Starling, Amanda B Muir, Gary W Falk, Bryson W Katona, Ben Z Stanger, Nancy R Zhang, Sydney M Shaffer <br/>
bioRxiv 2025.07.11.664452; doi: https://doi.org/10.1101/2025.07.11.664452* <br/>

## Contributors/Maintainers
Sydney Bracht & Jiazhen Rong

## License
Apache License 2.0
