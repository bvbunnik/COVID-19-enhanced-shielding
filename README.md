# COVID-19 Enhanced Shielding Model

Repository for R and C++ code used to model enhanced shielding measures on a simulated COVID-19 outbreak.

## Running the Code
### R 
R code should run without further modification of the code (after installing pre-requisite packages). Plots can be saved by altering the working directory in the `Cairo` functions at the bottom of the script.  

### C++

## Navigating the Repository 
### R
R code can be found in the `Enhanced Shielding` folder and is organised according to the code used to plot figures found in the manuscript. Current code is available for:
* Figure 2 - Baseline Analysis using "Central" Parameters
* Figure 3 - Sensitivity analysis for the ramp-up and ramp-down period for lockdown release
* Figure 4 - Sensitivity analysis of the lockdown trigger day and the Phase 2 R_e
* Figure 5 - Sensitivity analysis of the zeta parameter and Phase 1 R_e
* Figure S2 & 3 - Differing population structure

# C++

## Programs and Packages Used
COVID-19 modelling code was implemented using R and C++ (3.6.2). ODEs were solved using the "desolve" (1.27.1) package in R. 
Plotting in R was carried out using the "ggplot" package (3.3.0) in R.   

## Acknowledgements 
All analysis were performed by members and affiliates of Epigroup, University of Edinburgh: 
https://www.wiki.ed.ac.uk/display/Epigroup/Epigroup+-+Home

