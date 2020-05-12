# COVID-19 Enhanced Shielding Model

Repository for R and C++ code used to model enhanced shielding measures on a simulated COVID-19 outbreak.

## Running the Code
### R 
R code should run without further modification (after installing pre-requisite packages). 
The chosen directory where plots are saved can be altered by changing the `setwd()` function at the top of each script. Plots are automatically saved in the chosen working directory if each script is run. 

Code to run analysis of the population structure `Enhanced_Shielding_PopStruct_FigS3_5.R` (Figure S3 & 5) require the corresponding `.csv` files to be present in the User's chosen working directory.
The name of `.csv` correspond to each population structure analysed and represent simulations implemented in C++ and imported into R for graphical output.

### C++
C++ code can be found in main.cpp. This file can be compiled using g++. Implementation makes use of boost libraries (odeint).

## Navigating the Repository 
### R
R code can be found in the `Enhanced Shielding` folder and is organised according to the code used to plot figures found in the manuscript. **WORK IN PROGRESS** Current code is available for:
* Figure 2 & S1 - Baseline Analysis using "Central" Parameters and Recovered Fraction Plot
	* `Enhanced_Shielding_Baseline_Fig2_S1.R`
* Figure 3 - Sensitivity Analysis for Phase4 R_e, Zeta (1/Duration of Immunity), Phase 1 R_e and % Compliance of Vulnerable Population 
	* `Enhanced_Shielding_Phase4_Zeta_Phase1_Compliance_Fig3.R`
* Figure S2 - Sensitivity analysis of Phase 2 lockdown duration
	* `Enhanced_Shielding_LockdownDur_S2.R`
* Figure S3 - Sensitivity Analysis for Phase3 R_e Duration
	* `Enhanced_Shielding_Phase3_FigS3.R`
* Figure S4 & S5 - Sensitivity analysis of the lockdown trigger day and the Phase 2 R_e
	* `Enhanced_Shielding_TrigDay_Phase2_FigS4_S5.R`
* Figure S6 & 7 - Sensitivity analysis differing population structure and vulnerable/shielders ratio
	* `Enhanced_Shielding_PopStruct_FigS6_7.R`
* Figure S8 - SIS Baseline Plot
	* `Enhanced_Shielding_SIS_FigS8.R`
* Figure S9 - Analysis of efficacy of Shielders % testing
	* `Enhanced_Shielding_ShieldTest_FigS9.R`


## Programs and Packages Used
COVID-19 modelling code was implemented using R (3.6.2) and C++ independently. ODEs were solved using the `desolve` (1.27.1) package in R and `odeint` in C++. Plotting in R was carried out using the `ggplot` package (3.3.0). Dataframe manipulation was performed using `reshape2` (1.4.4). Finalised plot output was performed using `ggpubr` (0.2.4) and `Cairo` (1.5-10) packages. 

## Acknowledgements 
All analysis were performed by members and affiliates of Epigroup, University of Edinburgh: 
https://www.wiki.ed.ac.uk/display/Epigroup/Epigroup+-+Home
