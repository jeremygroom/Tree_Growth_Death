####### ----------------------------------------------------------#######
###           Data preparation and analysis code                   ###
####### ----------------------------------------------------------#######


# The purpose of this code is to allow the user to have control over the desired
#  analysis pathway.  There are six different analysis pathways, three each 
#  for tree mortality and basal area growth.  
#  1) Analysis of four climate variables in 9 climate variable/change variable
#     quadrants.
#  2) Analysis by size class. The user selects the DBH breakpoint and the 
#    TREE file is separated into two files. Each is analyzed for the 9 quadrants.
#       - The diameter at the first visit is used for filtering
#  3) Analysis by site class (SITECLCD, 7).  The 9 climate quadrants are each 
#   subdivided into 7 site classes.  This is a change in plot interpretation,
#    not TREE file.  There become effectively 63 quadrant conditions.  
#
#   The data preparation is combined with the analysis code because of the 
#   variety of analysis types.  The data preparation code does not separately
#   create output files for analysis code to operate upon because the file size
#   becomes unwieldy. This process will keep the Git repo size substantially
#   smaller.
#


#### 1) Loading constants, libraries, and functions----------------------------

source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))


#### 2) Base data prep --------------------------------------------------------
# These are the base data sets that will be transformed by the analysis types.



### ==> Load the tree/plot data ===================================================

# Working backwards: is the tree/plot data set ready to go? If so, let's unzip the CSV file and move ahead. 
#   If not, do the SQLite databases for CA/OR/WA need data extracted?  Those datasets are then processed by the 
# code to produce the tree/plot dataset.  
if(file.exists(file.path(DATA.LOC, "Distilled_Tree_Data.zip"))) {
  tree.plt.data <- read_csv(unzip(paste0(DATA.LOC, "Distilled_Tree_Data.zip"), "Distilled_Tree_Data.csv")) # This unzips the folder in the parent directory.
  file.remove("Distilled_Tree_Data.csv")
} else {
  ## Some info (harvest, fire deaths for trees) was absent in the range-shift analysis.  We need to introduce it here. 
  ##  The following code relies on downloaded and zipped SQLite FIA data from Oregon, Washington, and CA.  See the Global.R file:
  #    the zip files should be placed in the same folder, defined by SQL.LOC.  The following code unpacks, reduces, combines, and saves
  #   the PLOT, TREE, and COND files for each state as a zipped RDS file.  
  if (file.exists(paste0(DATA.LOC, "Addl_PlotTreeInfo.zip")) == FALSE) {
    source(paste0(CODE.LOC, "FIA_SQL_compile.R"))  # Data filtering for TREE occurs here as well. 
  }
  
  fia.tables <- read_rds(unzip(paste0(DATA.LOC, "Addl_PlotTreeInfo.zip"), "Addl_PlotTreeInfo.rds")) # This unzips the folder in the parent directory.
  
  # Deleting written RDS file (~200+ MB)
  if (file.exists("Addl_PlotTreeInfo.rds") == TRUE) {
    file.remove("Addl_PlotTreeInfo.rds")
  }
  
  # Now hopping over to run this code which will save the tree/plot zipped CSV file.
  source(paste0(CODE.LOC, "FIA_data_distillation.R"))
  
}

### ==> Load the climate data ===================================================
if (file.exists(paste0(DATA.LOC, "Climate_plot_results.rds")) == FALSE) {
  source(paste0(CODE.LOC, "Climate_AET_Preparation.R"))  # Data filtering for TREE occurs here as well. 
}

climate.data <- read_rds(paste0(DATA.LOC, "Climate_plot_results.rds")) # This unzips the folder in the parent directory.

### ==> Load the latitude/longitude data ===================================================
# Bringing in Lat/Lon from earlier study. Needed to recreate State_Plot (above) for the join.  
latlon <- read_csv(paste0(DATA.LOC, "PlotLatLon.csv")) %>% select(-n)

# Selecting summer AET, adding the variable State_Plot from the Groom/Vicente analysis 
# to enable joining of lat/lon info.
climate.use <- climate.data$summer_results %>% select (puid, pre_mean, difference) %>%
  rowwise() %>%
  mutate(plot = strsplit(puid, "_")[[1]][1], 
         state = strsplit(puid, "_")[[1]][2],
         State_Plot = as.numeric(paste0(plot, state))) %>%
  left_join(latlon, by = "State_Plot") %>%
  dplyr::select(-plot, -state, -State_Plot)

data.use <- left_join(tree.plt.data, climate.use, by = "puid")
dat.na <- data.use %>% filter(is.na(difference)) # 7 plots for which we don't have AET data
data.use <- anti_join(data.use, dat.na)



# ANALYSIS.PATHWAY 2 has two tree data sets.  We reduce 'data.use' to include
#  only those selected trees. DBH.CUTPOINT is set in Global.R

if (ANALYSIS.PATHWAY == 2) {
  
  data.use.S <- data.use %>% filter(PREVDIA < DBH.CUTPOINT)
  data.use.L <- data.use %>% filter(PREVDIA >= DBH.CUTPOINT)
  
}

# For the small/large tree analysis, the small tree dataset runs through the code
#  like ANALYSIS.PATHWAY 1 or 3.  Large trees are treated separately (see below).
treedat.use <- if(ANALYSIS.PATHWAY != 2) data.use else data.use.S

n_quant <- if(ANALYSIS.PATHWAY != 3) length(QUANT.LEVELS) else 7

PlotDat <- treedat.use %>% dplyr::select(STATECD, puid, ESTN_UNIT, STRATUMCD, w) %>% 
  distinct() %>%
  mutate(stratum = as.numeric(paste0(STRATUMCD, 0, ESTN_UNIT, 0, STATECD))) 


#pd.check <- PlotDat %>% select(State_Plot) %>% distinct() # Yes, same number of rows


#### 3) Specialized data prep --------------------------------------------------------
## First the data are summarized by species and climate variable (clim.mort.resp.fcn, 
#   clim.growth.resp.fcn).  Then the data are combined into arrays for 
#    analysis (parse.tree.clim.fcn).
tree.mort.dat <- map(sel.spp, clim.mort.resp.fcn, clim.var = "aet", treedat.sel = treedat.use, clim.dat = climate.use) 

# Create lists of tree growth and number of trees that grew in each plot
tree.grow.dat <- map(sel.spp, clim.growth.resp.fcn, clim.var = "aet", treedat.sel = treedat.use, clim.dat = climate.use) 

# Analysis pathway # 2 additionally creates these for larger trees
if(ANALYSIS.PATHWAY == 2) {
  
  tree.mort.dat.aet2 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "aet", treedat.sel = data.use.L) 
  tree.grow.dat.aet2 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "vpdmin", treedat.sel = data.use.L) 
  
}


# Combining data into arrays and such, preparing for analysis.
arrays.aet.mort1 <- parse.tree.clim.fcn(tree.mort.dat, "aet", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP, clim.dat = climate.use)

arrays.aet.grow1 <- parse.tree.clim.fcn(tree.grow.dat, "aet", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP, clim.dat = climate.use)


if(ANALYSIS.PATHWAY == 2) {
  arrays.aet.mort2 <- parse.tree.clim.fcn(tree.mort.dat.aet2, "aet", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
  arrays.aet.grow2 <- parse.tree.clim.fcn(tree.grow.dat.aet2, "aet", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
}


## Fire mortality: what proportion died from fire?  This info is only generated for pathway 1 and mortality.
if (ANALYSIS.PATHWAY == 1 & CALC.FIREPROP == TRUE) {
  fire.frac.table.fcn(tablename = "Fire_Prop_aet_mort.csv", tableloc = paste0(RESULTS1.LOC, "Mort_figs_aet/"), treedat = treedat.use, parseddat = arrays.aet.mort1)  
}



#### 4) Analysis and plotting  --------------------------------------------------------

# If desired, the analysis will obtain mortality estimates for the selected species
#  by state and overall. These values can be used to evaluate species for inclusion in the
#  main analysis
if(RUN.STATES == TRUE) {
  source(paste0(CODE.LOC, "Overall_Mort_Est.R"))
}


# The remaining code breaks up the analyses by type (mortality, growth) and climate variable (AET).
# If ANALYSIS.PATHWAY is 2 then the process is run twice, once for each size class. 
#  The plotting functions are located at the end of the process, with each ANALYSIS.PATHWAY
#  having a distinct function.

dataset.number <- if(ANALYSIS.PATHWAY == 2) 2 else 1 # How many times to run through.


tic() 
# For ANALYSIS.PATHWAY == 1, one climate variable, and 16 dedicated cores, the analysis takes about 8 minutes.

# Subscript i cycles through the climate variables. 
#for (i in 1:length(CLIM.VAR)){
i <- 1  # One climate variable at this time, 'aet'
for(j in 1:length(ANALYSIS.TYPE)) {
  
  # Setting up data set outputs, one for ANALYSIS.PATHWAY = 1, 3,
  #  two for ANALYSIS.PATHWAY 2.
  
  results <- list(results1 = NULL, results2 = NULL)
  
  for(k in 1:dataset.number) {
    
    
    var1 <- "pre_mean" #paste0("pre.", clim.var)
    var.delt <- "difference"  #paste0("delt.", clim.var)
    
    # Need climate names for files and axes.
    clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE)
    
    var.label <- clim.names$label[clim.names$values == var1]
    #var.filename <- clim.names$filename[clim.names$values == var1]
    var.delt.label <- clim.names$label[clim.names$values == var.delt]
    
    # Obtaining the data to work with: Grabbing the list objects from above
    mort.grow.dat <- get(paste0("arrays.", CLIM.VAR[i], ".", ANALYSIS.TYPE[j], k)) 
    
    vals_dat <- mort.grow.dat$vals_dat
    all_dat <- mort.grow.dat$all_dat
    quant.array <- mort.grow.dat$quant.array
    state.array <- mort.grow.dat$state.array
    state.n <- mort.grow.dat$state.n
    quant.matrix <- mort.grow.dat$quant.matrix
    centroid.array <- mort.grow.dat$centroid.array
    quant.n <- mort.grow.dat$quant.n 
    quant.lims <- mort.grow.dat$quant.lims
    quant.lims.delta <- mort.grow.dat$quant.lims.delt
    
    
    
    ## First, adjusting the species list
    #spp.id <- paste0("X", spp.list)     # Can use SEL.SPP
    spp.list <- as.numeric(gsub("X", "", SEL.SPP))
    
    strat.num <- length(unique(vals_dat$stratum))

    n.spp <- length(SEL.SPP)
    
    
    # Strata for bootstrap resampling procedure
    strata <- unique(vals_dat$stratum)
    strata.num <- vals_dat %>% dplyr::select(stratum, puid) %>%
      ungroup() %>% 
      mutate(row.id = row_number())
    
    
    plan(multisession, workers = n.cores) # Setting up parallel computing. See Global.R for n.cores.
    
    
    ## Quantile estimates for mortality or growth, across select species.
    quant.index <- 1:n_quant
    all.quants <- quant.index %>% 
      map(\(x) q_mort.grow.fcn(q1 = x, vals.dat = vals_dat, all.dat = all_dat, # Alternate mapping formulation, only advantage is clarity.
                               array.name = quant.array, category.n = quant.n,
                               selected.spp = SEL.SPP) )  
    quant.table <- bind_rows(all.quants) %>% arrange(Species, Quantile)
    
    results[[k]] <- list(quant.table = quant.table, var.label = var.label, var.delt.label = var.delt.label, 
                         quant.matrix = quant.matrix, quant.lims = quant.lims, quant.n = quant.n) 
    
  }   # end k 
  
  # Analysis Pathway 1: no subgroups of size or site class
  if (ANALYSIS.PATHWAY == 1) {
    # Plotting paired plots of mortality/growth by quantile and a scatterplot of plot distribution by quantiles.
    map(SEL.SPP, pair.plts.fcn, var.filename = "aet", quant.table = results$results1$quant.table, var.label = results$results1$var.label, 
        var.delt.label = results$results1$var.delt.label, var1 = var1, var.delt = var.delt,
        quant.matrix = results$results1$quant.matrix, quant.lims = quant.lims, quant.n = quant.n)
  }
  
  if (ANALYSIS.PATHWAY == 2) {
    quant.table.1.2 <- bind_rows(results$results1$quant.table, results$results2$quant.table, .id = "smaller_larger")  
    quant.matrix.1 <- results$results1$quant.matrix
    quant.matrix.2 <- results$results2$quant.matrix
    quant.lims.1 <- results$results1$quant.lims
    quant.lims.2 <- results$results2$quant.lims
    quant.n.1 <- results$results1$quant.n 
    quant.n.2 <- results$results2$quant.n # used for plotting scatterplot; just need to know which quantiles > 0
    
    map(SEL.SPP, pair.plts2.fcn, quant.table = quant.table.1.2, var.label = results$results1$var.label, 
        var1 = var1, var.delt = var.delt, var.delt.label = results$results1$var.delt.label, 
        quant.matrix.1 = quant.matrix.1, quant.matrix.2 = quant.matrix.2, 
        quant.lims.1 = quant.lims.1, quant.lims.2 = quant.lims.2,
        quant.n.1 = quant.n.1,
        quant.n.2 = quant.n.2)
    
  }
  
  if (ANALYSIS.PATHWAY == 3) {
    map(SEL.SPP, pair.plts3.fcn, quant.table = results$results1$quant.table, var.label = results$results1$var.label, 
        var.delt.label = results$results1$var.delt.label, var1 = var1, var.delt = var.delt,
        quant.matrix = quant.matrix, quant.n = quant.n)
  }
  
  
}   # end j
#}   # end i

toc()











