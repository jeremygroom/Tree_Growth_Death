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



##### ---- Constants ---- #####
DATA.LOC <- "Data/"
CODE.LOC <- "Analysis Code/"
RESULTS.LOC <- "Results/"

# Assuming WA/OR/CA zipped databases are in the same location
SQL.LOC <- "G:/My Drive/Consulting Practice/Contracts/ODF_FIA_2024/sqlDB_WA_CA/"

### Select spp for analysis
#SEL.SPP <- "X202" # Douglas fir, 8600 + 
#SEL.SPP <- "X98" # Sitka spruce
#SEL.SPP <- "X108" # Lodgepole pine
#SEL.SPP <- "X312" # Bigleaf maple
# SEL.SPP <- "X242" # Western Red Cedar
SEL.SPP <- c("X11", "X122", "X117", "X202", "X242", "X805", "X81", "X818", "X93") # List of 9 species I decided to focus upon
VAR1 <- "pre.vpdmin"
VAR.DELT <- "delt.vpdmin"
QUANT.PROBS <- c(0.25, 0.75)
QUANT.LEVELS <- c("LiHc", "MiHc", "HiHc", "LiMc", "MiMc", "HiMc", "LiLc", "MiLc", "HiLc")
N.PLOT.LIM <- 10  # Limit to the minimum number of plots for which an estimate will be calculated
BS.N <- 1000    # Bootstrap iteration number
 

## Stratum info - not sure we need this. 'orig' already has W_h, # of strata = length(unique(W_h))
strat <- read_csv(paste0(DATA.LOC, "strat_info052120.csv"), show_col_types = FALSE) %>%
  mutate(W_h = P1POINTCNT/p1pntcnt_eu)             # W_h is the stratum weight
strat2 <- strat %>% select(STRATUM, P1POINTCNT, W_h)  # Reducing the number of columns to those we need 
strat.num <- dim(strat2)[1]


n_quant <- length(QUANT.LEVELS)

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

