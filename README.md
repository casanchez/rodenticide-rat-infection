[![DOI](https://zenodo.org/badge/387547164.svg)](https://zenodo.org/badge/latestdoi/387547164)

A repository for: 
================

Murray MH, Sánchez CA. 2021. Urban rat exposure to anticoagulant rodenticides
and zoonotic infection risk. Biology Letters. doi: 10.1098/rsbl.2021.0311


**This repository has 5 R scripts used for this analysis:** The first script
(1_calculateAgeSMI.R) loads the data and calculates two new variables: age 
(based on weight) and scaled mass index (SMI), a measure of body condition. 
The second script (2_analysis.R) runs models to assess how variance in infection 
status is associated with multiple predictor variables. The third script 
(3_Figure2.R) generates Figure 2 of the manuscript. Figure 1 was generated 
separately in QGIS. The other scripts contain libraries and functions necessary
to run the rest of the code, and are sourced when needed.
