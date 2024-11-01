## Global conditions for all scripts to use ##


##### ---- Libraries ---- #####
library(tidyverse)
library(readr)
library(viridis)
library(cowplot)
library(patchwork)
library(readxl)

library(furrr)
library(parallel)
library(RSQLite) # For obtaining SQLite FIA databases
library(tictoc)  # For development, timing routines.



##### ---- Constants ---- #####

# File locations:
DATA.LOC <- "Data/"
CODE.LOC <- "Analysis Code/"
RESULTS.LOC <- "Results/"
RESULTS1.LOC <- "Results/Quantile_Only/"  # For ANALYSIS.PATHWAY 1
RESULTS2.LOC <- "Results/Size_Class/"  # For ANALYSIS.PATHWAY 2
RESULTS3.LOC <- "Results/Site_Class/"  # For ANALYSIS.PATHWAY 3



# Assuming WA/OR/CA zipped databases are in the same location (outside of GitHub clone)
SQL.LOC <- "G:/My Drive/Consulting Practice/Contracts/ODF_FIA_2024/sqlDB_WA_CA/"

# Which analysis to do??
# Pathway 1 = Mortality/Growth analysis by 9 quantiles, 4 climate variables.
# Pathway 2 = Pathway 1 run for each of 2 tree size classes (< / >= DBH cutpoint)
# Pathway 3 = Mortality/Growth analysis by 9 quantiles * 7 site class categories,
#                4 climate variables.
ANALYSIS.PATHWAY <- 2   # 1, 2, or 3
DBH.CUTPOINT <- 12      # DBH inches. Examine trees above & below this amount.

### Select spp for analysis
SEL.SPP <- c("X11", "X122", "X117", "X202", "X242", "X805", "X81", "X818", "X93") # List of 9 species I decided to focus upon
CLIM.VAR <- c("vpdmin", "vpdmax", "temp", "precip")  # Climate variables under examination.
ANALYSIS.TYPE <- c("grow", "mort")

#QUANT.DELTA.PROBS.POS <- 0.33 # Delta climate variable: positive values broken between this quantile. 
#QUANT.DELTA.PROBS.NEG <- 0.33 # Delta climate variable: negative values broken between this quantile.
VAR.DELTA.BOUNDARIES <- tibble(clim.var = CLIM.VAR, max.min = c(0.5, 0.5, 0.5, 100)) # For setting absolute +/- boundaries for the listed variables


QUANT.PROBS <- c(0.25, 0.75)  # First-visit climate variable values: Broken between these values into three parts.

QUANT.LEVELS <- c("LiHc", "MiHc", "HiHc", "LiMc", "MiMc", "HiMc", "LiLc", "MiLc", "HiLc")
N.PLOT.LIM <- 10  # Limit to the minimum number of plots for which an estimate will be calculated
BS.N <- 200    # Bootstrap iteration number
 

# Species number as number
sel.spp <- as.numeric(gsub("X", "", SEL.SPP))

## Stratum info - not sure we need this. 'orig' already has W_h, # of strata = length(unique(W_h))
strat <- read_csv(paste0(DATA.LOC, "strat_info052120.csv"), show_col_types = FALSE) %>%
  mutate(W_h = P1POINTCNT/p1pntcnt_eu)             # W_h is the stratum weight
strat2 <- strat %>% select(STRATUM, P1POINTCNT, W_h)  # Reducing the number of columns to those we need 
strat.num <- dim(strat2)[1]


n_quant <- length(QUANT.LEVELS)
quant.level.table <- tibble(var.deltvar = factor(QUANT.LEVELS, level = QUANT.LEVELS), q.num = 1:length(QUANT.LEVELS))

# From the previous analysis (Groom and Monleon 2023)
orig <- read_csv(paste0(DATA.LOC, "Occ_OriginalVisit2.csv")) %>%
  dplyr::select(STATECD, PLOT_FIADB, State_Plot, STRATUM, W_h, propfor, propnonfor, propwater, propmiss, intensification, starts_with("X"))

# Number of possible species, indexes, names list
# Data frame of plots and species considered, and plot weights

spp.list <- as.numeric(gsub("X", "", names(orig)[(ncol(orig) - (length(grep("X", names(orig))) - 1)):ncol(orig)]))
spp.id <- names(orig)[(ncol(orig) - (length(grep("X", names(orig))) - 1)):ncol(orig)]

spp.names <- read_xlsx(paste0(DATA.LOC, "FullSppNames.xlsx")) %>% filter(SPCD %in% spp.list)

## Items for parallel computing
n.cores <- round(detectCores() * 0.75, 0) # Number of cores, package "parallel".  Using 75% of the available cores, rounded.

