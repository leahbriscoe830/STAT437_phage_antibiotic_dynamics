# Modelling the bacteriophage-antibiotic synergistic effect on Enterococcus faecalis

## About
This repo modifies the model created by Leclerc et al. to model the synergistic effect of bacteriophages and antibiotics on Enterococcus faecalis. This project was created by Zoe Berge, Shreya Bhagi, and Leah Briscoe for STAT437 Quantitative Bioinformatics at Loyola University Chicago.

## Model

The core model code is in the "/Model" folder. This folder contains an R script called efaecalis_model_final.R encoding a function for the model and the code for ef_fig1. This relies on the `deSolve` package and must be run before any analysis scripts. 

## Parameters

The "/Parameters" folder contains an RDS object containing the parameters defined in the model script. To adjust the parameters, change them in the model script and rerun the model to save the new parameters.

## Analysis

The "/Analysis" folder contains various scripts used to create the figures shown in the paper. They are named to indicate the figure they generate.

## Figures

The "/Figures" folder contains the figures created by the different scripts. Each figure is named according to its position in the paper (e.g. "ef_fig1" is Figure 1 in the paper).

## References 
Leclerc, Q. J., Lindsay, J. A., & Knight, G. M. (2022). Modelling the synergistic effect of bacteriophage and antibiotics on bacteria: Killers and drivers of resistance evolution. PLoS computational biology, 18(11), e1010746. https://doi.org/10.1371/journal.pcbi.1010746 

Moryl, M., Szychowska, P., Dziąg, J., Różalski, A., & Torzewska, A. (2024). The Combination of Phage Therapy and β-Lactam Antibiotics for the Effective Treatment of Enterococcus faecalis Infections. International journal of molecular sciences, 26(1), 11. https://doi.org/10.3390/ijms26010011 

## License

This work is distributed under the MIT license (see LICENSE file).
