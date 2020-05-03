# COVID-19 Enhanced Shielding Model

Repository for R and C++ code used to model enhanced shielding measures on a simulated COVID-19 outbreak.

## Running the Code
### R 
R code should run without further modification (after installing pre-requisite packages). 
The chosen directory where plots are saved can be altered by changing the 'setwd()' function at the top of each script. Plots are automatically saved in the chosen working directory if each script is run. 

Code to run analysis of the population structure 'Enhanced_Shielding_PopStruct_FigS3_5.R' (Figure S3 & 5) require the corresponding '.csv' files to be present in the User's chosen working directory.
The name of '.csv' correspond to each population structure analysed and represent simulations implemented in C++ and imported into R for graphical output.

### C++
**WIP**

## Navigating the Repository 
### R
R code can be found in the `Enhanced Shielding` folder and is organised according to the code used to plot figures found in the manuscript. **WORK IN PROGRESS** Current code is available for:
* Figure 2 - Baseline Analysis using "Central" Parameters
	* `Enhanced_Shielding_Baseline_Fig2.R`
* Figure 3 - Sensitivity analysis for the ramp-up and ramp-down period for lockdown release
	* `Enhanced_Shielding_Phase3_Fig3.R`
* Figure 4 - Sensitivity analysis of the lockdown trigger day and the Phase 2 R_e
	* `Enhanced_Shielding_TrigDay_Phase2_Fig4.R`
* Figure 5 - Sensitivity analysis of the zeta parameter and Phase 1 R_e
	* `Enhanced_Shielding_Zeta_Phase1R0_Fig5.R`
* Figure 8 - Sensitivity analysis of Phase 2 lockdown duration
	* `Enhanced_Shielding_LockdownDur_Fig8.R`
* Figure 9 - Analysis of Vulnerable population enhanced shielding % compliance (Phase 4)
	* `Enhanced_Shielding_Compliance_Fig9.R`
* Figure 10 - Sensitivity analysis of all Phase 4 Betas 
	* `Enhanced_Shielding_Phase4_Fig10.R`
* Figure 12 - Analysis of efficacy of Shielders % testing
	* `Enhanced_Shielding_ShieldTest_Fig12.R`
* Figure S2 - Sensitivity analysis SIS model w/ baseline parameters
	* `Enhanced_Shielding_SIS_FigS2.R`
* Figure S3 & 5 - Sensitivity analysis differing population structure
	* `Enhanced_Shielding_PopStruct_FigS3_5.R`

### C++
**WIP**

## Programs and Packages Used
COVID-19 modelling code was implemented using R and C++ (3.6.2). Main Packages used in R include: ODEs were solved using the `desolve` (1.27.1) package in R. 
Plotting in R was carried out using the `ggplot` package (3.3.0) in R. Dataframe manipulation was performed using `dplyr` (0.8.3). Finalised outputting was performed
using `ggpubr` (0.2.4) and `Cairo` (1.5-10) packages. 

## Acknowledgements 
All analysis were performed by members and affiliates of Epigroup, University of Edinburgh: 
https://www.wiki.ed.ac.uk/display/Epigroup/Epigroup+-+Home
