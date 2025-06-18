#### ----------------------------------------------------------------- ####
###  Creation of data set from FIA downloads for CA, WA, OR
#### ----------------------------------------------------------------- ####



# This code may be sourced from Data_Prep_Analysis.R and will likely produce
#  an RDS file. Its purpuse is to directly construct a dataset from the 
#  FIA data that includes REMPER, AGENTCD, intensification plots, etc.

PLOT <- fia.tables$plot_vals
COND <- fia.tables$cond_vals
TREE <- fia.tables$tree_vals
SCCM <- fia.tables$sccm_vals %>% dplyr::select(-CREATED_DATE) # From the table SUBP_COND_CHNG_MTRX, tracking changes in portions of subplots sampled.
POP_ESTN_UNIT <- fia.tables$peu_vals
POP_EVAL <- fia.tables$pev_vals
POP_EVAL_TYP <- fia.tables$pet_vals
POP_PLOT_STRATUM_ASSGN <- fia.tables$ppsa_vals 
POP_STRATUM <- fia.tables$popstrm_vals

rm(fia.tables)



data <- PLOT %>%
  ## Add a PLT_CN column for easy joining
  mutate(PLT_CN = CN) %>%
  ## Join COND & TREE
  left_join(COND, by = c('PLT_CN', 'STATECD', 'PLOT')) %>%
  left_join(TREE, by = c('PLT_CN', 'PLOT', 'CONDID', 'STATECD')) %>%  
  ## Population tables
  left_join(POP_PLOT_STRATUM_ASSGN, by = c('PLT_CN')) %>% # Links PLT_CN to STRATUM_CN
  left_join(POP_STRATUM, by = c('STRATUM_CN' = 'CN', 'STATECD')) %>% # Links strata info to PLOT, ultimately
  left_join(POP_ESTN_UNIT, by = c('ESTN_UNIT_CN' = 'CN', 'STATECD')) %>%
  left_join(POP_EVAL, by = c('EVAL_CN' = 'CN', 'STATECD')) %>%
  mutate(stratum = as.numeric(paste0(STRATUMCD, 0, ESTN_UNIT, 0, STATECD))) 


# Detour: Obtaining all plots for use in determining plot-level weights (as opposed to area weights)
plotN <- data %>% mutate(puid = paste(PLOT, STATECD, COUNTYCD, sep = "_")) %>%
  filter(INTENSITY == "1", is.na(STRATUMCD) == FALSE) %>%
  dplyr::select(puid, stratum, STATECD, ESTN_UNIT, STRATUMCD, INTENSITY, EVALID, END_INVYR) %>% 
  distinct()
  
# Back to the tree-plot dataset: The joining of TREE to PLOT and COND introduced NAs into DIA2
data <- data %>% 
  filter(!is.na(STRATUM_CN), !is.na(DIA2))            
  


# Adding variable tAdj, which is the adjustment factor for partially unsampled plots.  Here the process treats 
#  all trees in the subplot and macroplot greater in DBH than the macro breakpoints as belonging to the macro plots.  
#  If the trees initially had DIA measured below the breakpoints they were only measured in the subplots.  Thus,
#  they remain as subplot trees.  Conversely, trees larger than the breakpoints in the subplots are treated as belonging 
#  to the macroplot the whole time.
data <- data %>% 
  mutate(
    ## TREE
    treeType = case_when(
      MACRO_BREAKPOINT_DIA == 30 & PREVDIA < 30 ~ "subplot",
      MACRO_BREAKPOINT_DIA == 24 & PREVDIA < 24 ~ "subplot",
      ## DIA is greater than macro breakpoint diameter
      MACRO_BREAKPOINT_DIA == 30 & PREVDIA >=30 ~ "macro",
      MACRO_BREAKPOINT_DIA == 24 & PREVDIA >=24 ~ "macro",
      .default = NA
    ),
    tAdj = case_when(
      treeType == "macro" ~ ADJ_FACTOR_MACR,
      treeType == "subplot" ~ ADJ_FACTOR_SUBP,
      .default = NA
    ),
    # Dead trees without DBH values the 2nd visit are not given TPA. Setting these based on PREVDIA.
    TPA_UNADJ2 = case_when(
      !is.na(TPA_UNADJ) ~ TPA_UNADJ, 
      is.na(TPA_UNADJ) & treeType == "subplot" ~ 6.018046,
      is.na(TPA_UNADJ) & treeType == "macro"  ~ 0.999188,
      .default = NA
    ),
    puid = paste(PLOT, STATECD, COUNTYCD, sep = "_")
  ) %>%
  # Removing non-forest-condition CONDID
  filter(COND_STATUS_CD == 1) %>%
  # Will want to not filter by INTENSITY and check the EVALID exclusion when we tighten up the analysis strata.
  filter(INTENSITY == "1", EVALID != 61603)      # This is a CA EVALID with an END_INVYR earlier (2016) than those used for OR and WA.    

#dat.1.n <- nrow(data)


## There are some datasets with three visits.  We want the last revisit set. Removes 200 plots.
dat.revisit <- data %>% dplyr::select(puid, MEASYEAR) %>% distinct() %>%
  group_by(puid) %>%
  summarize(last.yr = max(MEASYEAR))

data <- data %>% left_join(dat.revisit, by = "puid") %>%
  filter(MEASYEAR == last.yr) %>%
  dplyr::select(-last.yr) 

#dat.2.n <- nrow(data)  

# Looking for area adjustments.  Combining PLOT, SCCM, and COND, reducing to the 
#  measurement years of interest, and then adding the subplot (1-4) values for subplots (small diameter)
#  and macroplots. The goal is to get the overall plot proportion for sub/macro plots that are available as forest.
pltA <- PLOT %>% mutate("PLT_CN" = CN) %>%
  dplyr::select(-CN) %>%
  left_join(SCCM, by = c("PLT_CN", "STATECD")) %>%
  left_join(COND, by = c("PLT_CN", "STATECD", "CONDID", "PLOT")) %>%
  mutate(puid = paste(PLOT, STATECD, COUNTYCD, sep = "_"))

data.yrs <- data %>% dplyr::select(puid, MEASYEAR) %>%
  group_by(puid) %>%
  reframe(data.yrs = max(MEASYEAR))

pltA2 <- pltA %>% left_join(data.yrs, by = "puid") %>%
  filter(MEASYEAR == data.yrs, COND_STATUS_CD == 1, is.na(CONDPROP_UNADJ) == FALSE) %>% 
  dplyr::select(puid, SUBP, SUBPTYP, SUBPTYP_PROP_CHNG) %>%
  group_by(puid, SUBP, SUBPTYP) %>%
  reframe(sum.prop = sum(SUBPTYP_PROP_CHNG) * 0.25) %>% 
  group_by(puid, SUBPTYP) %>%
  reframe(tot.sum.prop = sum(sum.prop)) %>%
  pivot_wider(names_from = "SUBPTYP", values_from = "tot.sum.prop") %>%
  rename("sub.prop" = `1`, "macro.prop" = `3`)


data <- left_join(data, pltA2, by = 'puid')

# Tests to see if there are NA values for sub.prop or macro.prop where we need them.
# --> Turns out there are no conditions where these are true. 
#data %>% filter(is.na(sub.prop), treeType == "subplot")
#data %>% filter(is.na(macro.prop), treeType == "macro")


### Key bit: figuring out how many trees each record represents. We'll multiply the 
#  tree by the macro/sub plot area (macroplots total 1 acre, about) and divide
#  by the appropriate proportion of macroplot or subplot available.  
data <- data  %>% 
  mutate(ntree.rep = ifelse(treeType == "macro", TPA_UNADJ2 / macro.prop,
                            TPA_UNADJ2 / sub.prop))



## Weights ------------------------------------

#length(unique(data$P1PNTCNT_EU))
total.p1 <- sum(unique(data$P1PNTCNT_EU))

# Create weights:
strata.w <- data %>%
  group_by(STATECD, ESTN_UNIT, STRATUMCD, P1POINTCNT, P1PNTCNT_EU) %>%   
  reframe(w =first(P1POINTCNT) / total.p1)

# Verify that the strata points summed equal the EU point count for the state.  
#wa <- strata.w %>% filter(STATECD == 53)
#sum(wa$P1POINTCNT) / sum(unique(wa$P1PNTCNT_EU))  # 0.99949

#or <- strata.w %>% filter(STATECD == 41)
#sum(or$P1POINTCNT) / sum(unique(or$P1PNTCNT_EU))  # 1.000

#ca <- strata.w %>% filter(STATECD == 6)
#sum(ca$P1POINTCNT) / sum(unique(ca$P1PNTCNT_EU))  # 0.9929


tree.plt.data <- left_join(data, strata.w, by = c("STATECD", "ESTN_UNIT", "STRATUMCD", "P1POINTCNT", "P1PNTCNT_EU"))
#dim(data %>% filter(is.na(w)))  # Check to verify all trees have weights
# nrow(data) == dat.2.n   # Verify no trees lost

# Plot weights: making sure the same EVALID is used in both data sets.
plotN2 <- plotN %>% filter(is.na(END_INVYR) == FALSE, EVALID %in% unique(tree.plt.data$EVALID), STRATUMCD %in% tree.plt.data$STRATUMCD)
write_csv(plotN2, file.path(DATA.LOC, "N_Plots.csv"))  # For plot-level weights.  <2MB, not zipping

# The main data set. Big. Writing then zipping.
write_csv(tree.plt.data, file.path(DATA.LOC, "Distilled_Tree_Data.csv"))
zip(zipfile = file.path(DATA.LOC, "Distilled_Tree_Data.zip"), files = file.path(DATA.LOC, "Distilled_Tree_Data.csv"),
    flags = '-r9Xj' )  # This weird bit keeps the parent directory from being included in the zip folder

# Deleting written CSV file (~200 MB), cleaning up after the data zip process. 
if (file.exists(file.path(DATA.LOC, "Distilled_Tree_Data.csv")) == TRUE) {
  file.remove(file.path(DATA.LOC, "Distilled_Tree_Data.csv"))
}


