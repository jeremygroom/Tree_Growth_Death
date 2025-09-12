## Global conditions for all scripts to use ##


##### ---- Libraries ---- #####
library(readr)
library(viridis)
library(cowplot)
library(patchwork)
library(readxl)
library(tidyverse)
library(here)         # File path management
library(furrr)
library(parallel)
# Generation of maps
library(sf)
library(ggspatial)
library(maps)



##### ---- Constants ---- #####

# Code controls:
CALC.FIREPROP <- TRUE   # Set this to TRUE if wanting to calculate or recalculate
                          # the proportion of dead trees killed by fire by quantile.
RUN.STATES <- TRUE # Should Data_Prep_Analysis.R find mortality estimates by state and overall?

RUN.SUMMARY <- TRUE # Should Data_Prep_Analysis.R calculate 

# File locations:
DATA.LOC <- "Data/"
CODE.LOC <- "Analysis Code/"
RESULTS.LOC <- "Results/"
RESULTS1.LOC <- "Results/Quantile_Only/"  # For ANALYSIS.PATHWAY 1
RESULTS.OTHER <- "Results/Other_Results/"  # For summaries, maps, etc.



# Assuming WA/OR/CA zipped databases are in the same location (outside of GitHub clone)
SQL.LOC <- "G:/My Drive/Consulting Practice/Contracts/ODF_FIA_2024/sqlDB_WA_CA/"
CLIMATE.LOC <- "G:/My Drive/Consulting Practice/Contracts/ODF_FIA_2024/TerraClimate data/"

# Which analysis to do??
# Pathway 1 = Mortality/Growth analysis by 9 quantiles, 4 climate variables.
# Pathway 2 = Pathway 1 run for each of 2 tree size classes (< / >= DBH cutpoint)
# Pathway 3 = Mortality/Growth analysis by 9 quantiles * 7 site class categories,
#                4 climate variables.
ANALYSIS.PATHWAY <- 1   # 1, 2, or 3
DBH.CUTPOINT <- 12      # DBH inches. Examine trees above & below this amount.

### Select spp for analysis
#SEL.SPP <- c("X11", "X122", "X117", "X202", "X242", "X805", "X81", "X818", "X93") # Original list of 9 species I decided to focus upon
#SEL.SPP <- c("X15", "X17", "X19", "X93", "X108", "X117", "X351", "X361", "X631", "X818") # List of 10 species with overall mortality > 0.15 and number of plots > 500
SEL.SPP <- c("X11", "X15", "X17", "X19", "X64", "X73", "X81", "X93", "X108", "X116", "X117", "X119", "X122", "X202", "X242", "X263", "X264", "X312", "X351", "X361", "X631", "X805", "X818") # Species w/ > 400 plots
CLIM.VAR <- c("aet", "cwd", "vpdmin", "vpdmax", "temp", "precip")  # Climate variables under examination.
CLIM.VAR.USE <- "cwd" # "aet" "pet"
CLIM.SUMMARY <- "summer_results"# "annual_results"       Select annual or summer climate summaries.

ANALYSIS.TYPE <- c("grow", "mort")


VAR.DELTA.BOUNDARIES <- tibble(clim.var = CLIM.VAR, max.min = c(3, 5, 0.5, 0.5, 0.5, 100)) # For setting absolute +/- boundaries for the listed variables

USE.QUANT.PROBS <- TRUE # TRUE = quantiles (defined in QUANT.PROBS) used to define initial climate variable breaks.
THREE.CHANGE.CATEGORIES <- FALSE # TRUE = Delta climate variable (e.g., change in CWD) has 3 categories.  FALSE = 2 categories.
QUANT.PROBS <- c(0.25, 0.75)  # First-visit climate variable values: Broken between these values into three parts.
INIT.CLIM.BREAKS <- c(75, 125) # if USE.QUANT.PROBS are FALSE, these are the domain limits for the initial climate variable values.

if(THREE.CHANGE.CATEGORIES) {
DOMAIN.LEVELS <- c("LiHc", "MiHc", "HiHc", "LiMc", "MiMc", "HiMc", "LiLc", "MiLc", "HiLc")
} else {
DOMAIN.LEVELS <- c("AL", "AM", "AH", "BL", "BM", "BH") # Above/Below threshold and Low/Medium/High
}
N.PLOT.LIM <- 10  # Limit to the minimum number of plots for which an estimate will be calculated
BS.N <- 1000    # Bootstrap iteration number


## Shiny app settings
SHINYAPP.IN.USE <- FALSE
SHINY.FONTSIZE <- 14

# Species number as number
sel.spp <- as.numeric(gsub("X", "", SEL.SPP))

## Stratum info - not sure we need this. 'orig' already has W_h, # of strata = length(unique(W_h))
#strat <- read_csv(paste0(DATA.LOC, "strat_info052120.csv"), show_col_types = FALSE) %>%
#  mutate(W_h = P1POINTCNT/p1pntcnt_eu)             # W_h is the stratum weight
#strat2 <- strat %>% select(STRATUM, P1POINTCNT, W_h)  # Reducing the number of columns to those we need 
#strat.num <- dim(strat2)[1]


n_domain <- length(DOMAIN.LEVELS)
domain.level.table <- tibble(var.deltvar = factor(DOMAIN.LEVELS, level = DOMAIN.LEVELS), q.num = 1:length(DOMAIN.LEVELS))

# From the previous analysis (Groom and Monleon 2023)
orig <- read_csv(paste0(DATA.LOC, "Occ_OriginalVisit2.csv")) %>%
  dplyr::select(STATECD, PLOT_FIADB, State_Plot, STRATUM, W_h, propfor, propnonfor, propwater, propmiss, intensification, starts_with("X"))

# Number of possible species, indexes, names list
# Data frame of plots and species considered, and plot weights

spp.list <- as.numeric(gsub("X", "", names(orig)[(ncol(orig) - (length(grep("X", names(orig))) - 1)):ncol(orig)]))
spp.id <- names(orig)[(ncol(orig) - (length(grep("X", names(orig))) - 1)):ncol(orig)]

spp.names <- readxl::read_xlsx(paste0(DATA.LOC, "FullSppNames.xlsx")) %>% dplyr::filter(SPCD %in% spp.list)


## Items for parallel computing
n.cores <- round(detectCores() * 0.75, 0) # Number of cores, package "parallel".  Using 75% of the available cores, rounded.


   # State layers
states <- map_data("state")   # Loaded with ggplot2
west_df <- subset(states, region == "california" | region == "oregon" | region == "washington")

