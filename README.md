# Tree_Growth_Death

This repository contains all background code and much of the data files behind 
the manuscript *Tree growth, mortality, and climatic water deficit in west-coast states, USA*
by Jeremiah Groom ([Groom Analytics LLC](https://www.groomanalytics.com/)) and 
Bryce Frank (USDA Forest Service).  It also contains files for the supplemental material and 
the code for its accompanying [Shiny dashboard](https://groomanalyticsllc.shinyapps.io/FIA_Growth_Mortality/).   

The manuscript describes the estimation of growth and mortality rates for 23 tree
species using FIA plot revisit data.  Estimates are provided for species overall,
by state, and by domains defined by plot initial visit and change in climatic water deficit (CWD). 

## Structure

The main directory contains the draft manuscript, draft supplemental materials, associated
RMD files, code for the Shiny app (`app.r`), the R Project file (`Tree_Growth_Death.Rproj`), and 
a file of settings and connections (`Global.R`).    

The RMD files for the manuscript and supplemental materials  source a number of outputs from other code 
so I won't describe them here.  

If you wish to trace the analysis path, let's head to the folder `Analysis Code`.

## Analysis Code

This folder contains a number of R scripts, only some of which are still useful.   

*Data_Prep_Analysis.R*:  This is a very good place to start.  The code is fairly well documented. 
It will source `Global.R` and `Functions.R`.  It will try to load code unless code doesn't exist,
in which case it will source files.  Specifically, it sources `FIA_SQL_compile.R` to grab 
downloaded FIA SQLite databases and process a bit, and then `FIA_data_distillation.R` to finish 
the data processing (well, enough so that files can be saved and worked from).  It brings in 
climate data, LATLON data, and others to round things out.  The code will produce several things: 

* A mortality table that provides a breakdown of mortality causes by CWD domain
* Estimates overall and by state for tree species growth and mortality (sources `Overall_Mort_Est.R`)
* CWD CWD domain estimates for growth and mortality
* Provides output files specifically for the manuscript via sourcing `Manuscript_information.R`.
* Produces some plots used in the manuscript.

If you wish to fully run this code you will need to download the FIA SQLite databases for Oregon, California, and 
Washington and point the code to your saved location in the `Global.R` script.

*Nplot_Climvar_analysis*: If the code were cleaner, this script would have been run as part of the previous 
file.  A number of important things happen, along with data exploration.  This is where we determine how many
species to analyze (minimum of 400 plots), and where we examine the relationships between CWD and delta CWD. 
This output is mentioned in the manuscript and displayed in the Supplemental Material document.

*terraclimate_download.R*: Here's how the TerraClimate data can be obtained via R.  

## Data
Data summaries and intermediary products are placed here, along with some outside files (e.g., `Occ_OriginalVisit2.csv`
 was originally from the Groom and Monleon 2023 analysis).  
 
## Results
 This folder contains the subfolders `Other_Results` and `Quantile_Results`.  The folders are populated as the code is run.  


