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

# Climate variables and reduction to a single table
tmp <- read_csv(paste0(DATA.LOC, "tmp.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.temp = post.temp - pre.temp) 
precip <- read_csv(paste0(DATA.LOC, "precip.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.precip = post.precip - pre.precip)
vpdmin <- read_csv(paste0(DATA.LOC, "vpdmin.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.vpdmin = post.vpdmin - pre.vpdmin)
vpdmax <- read_csv(paste0(DATA.LOC, "vpdmax.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.vpdmax = post.vpdmax - pre.vpdmax)


# clim.dat = all climate data
clim.dat <- list(tmp, precip, vpdmin, vpdmax) %>% reduce(left_join, by = c("STATECD", "PLOT_FIADB")) %>% 
  mutate(State_Plot = paste0(PLOT_FIADB, STATECD),
         State_Plot = as.numeric(State_Plot)) %>% 
  dplyr::select(-STATECD, -PLOT_FIADB) %>% 
  left_join(latlon, by = 'State_Plot')  # latlon loaded in Global.R .  For mapping.



##### Tree and plot info  ####
## Some info (harvest, fire deaths for trees) was absent in the range-shift analysis.  We need to introduce it here. 
##  The following code relies on downloaded and zipped SQLite FIA data from Oregon, Washington, and CA.  See the Global.R file:
#    the zip files should be placed in the same folder, defined by SQL.LOC.  The following code unpacks, reduces, combines, and saves
#   the PLOT, TREE, and COND files for each state as a zipped RDS file.  
if (file.exists(paste0(DATA.LOC, "Addl_PlotTreeInfo.zip")) == FALSE) {
  source(paste0(CODE.LOC, "FIA_SQL_compile.R"))
}

fia.tables <- read_rds(unzip(paste0(DATA.LOC, "Addl_PlotTreeInfo.zip"), "Addl_PlotTreeInfo.rds")) # This unzips the folder in the parent directory.

# Deleting written RDS file (~80 MB)
if (file.exists("Addl_PlotTreeInfo.rds") == TRUE) {
  file.remove("Addl_PlotTreeInfo.rds")
}

## Cleaned tree data, removing trees beyond 24 feet (only subplot data)
Tree1 <- readr::read_csv(unzip(paste0(DATA.LOC, "Cleaned_Trees_2019v2.zip"), "Cleaned_Trees_2019v2.csv")) %>% filter(DIST <= 24, SDN != 3)
# SDN = 1 is alive-alive, SDN = 2 = alive-dead.  We are removing all ingrowth.

# Deleting extracted csv file (~125 MB)
if (file.exists("Cleaned_Trees_2019v2.csv") == TRUE) {
  file.remove("Cleaned_Trees_2019v2.csv")
}
# Obtaining the (2024) current TREE data for OR/WA/CA to grab the AGENTCD values
fia.tree <- fia.tables$tree_vals %>% filter(INVYR > 2010 & INVYR < 2020) %>%
  select(-INVYR) %>%                      # Need to leave in to double-match trees to harvest events.
  mutate(CN = as.numeric(CN), 
         join.test = 1)

# Obtaining plot SITECLCD values. Note that the preparation of this file happened outside of this repository.
# (Bryce Frank's GitHub repo: https://github.com/bfrank5/site_class)
fia.cond <- read_rds(paste0(DATA.LOC, "plot_siteclcd.RDS")) %>%
  mutate(PLT_CN = as.numeric(PLT_CN)) %>%
  select(-MAICF_plot, -sum_p) %>%
  left_join(Tree1 %>% filter(INVYR < 2011) %>%  # Need to attach STATECD and PLOT to Bryce's file.
              select(PLT_CN, PLOT_FIADB, STATECD) %>%    # Just grabbing 3 columns from Tree1.
              distinct(), by = "PLT_CN") %>%
  select(-PLT_CN)

#fia.cond %>% filter(is.na(PLOT_FIADB))  # Check: any missing State_Plot vals? 

Tree1.1 <- Tree1 %>% left_join(fia.tree, by = c("STATECD", "PLOT_FIADB" = "PLOT", "TREE", "SUBP")) 

# Check: Did the join work as intended? Which trees were harvested, using STATUSCD and AGENTCD?
#sum(Tree1.1$join.test, na.rm = TRUE) # All fia.tree entries were successfully joined with Tree1 trees.
#table(Tree1.1$STATUSCD, Tree1.1$AGENTCD)   # 32073 trees burned (12609), were harvested (15260), or both (4204)
#tree.a80 <- Tree1.1 %>% filter(AGENTCD == 80, STATUSCD != 3)   # 4,204 trees = AGENTCD == 80, STATUSCD == 2
#tree.s3 <- Tree1.1 %>% filter(AGENTCD != 80, STATUSCD == 3)   # No trees were coded as harvest without AGENTCD == 80

tree.cut <- Tree1.1 %>% filter(AGENTCD == 80, STATUSCD == 3) %>% select(State_Plot_Subp_Tree)   # Extracting only trees that were cut. Burned remain.

# Removing cut/cleared trees from the tree list.
Tree1.2 <- Tree1.1 %>% anti_join(tree.cut, by = "State_Plot_Subp_Tree") %>% # First removing all harvested trees.
  filter(STATUSCD != 0, !RECONCILECD %in% 1:3)      # Removing trees that shrank below 5" DIA,
                                                    #  ingrowth, throughgrowth, missed alive/dead

# Removing remaining single tree records
single.trees <- Tree1.2 %>% group_by(SUBP, TREE, STATECD, PLOT_FIADB) %>%
  summarize(n = n()) %>%
  filter(n == 1)

Tree1.2 <- Tree1.2 %>% anti_join(single.trees, by = c("SUBP", "TREE", "STATECD", "PLOT_FIADB"))

# ANALYSIS.PATHWAY 2 has two tree data sets.  We need to pare the first visit 
#  trees by the DBH.CUTPOINT value and then reduce the Tree1.2 data set to include
#  only those selected trees

if (ANALYSIS.PATHWAY == 2) {
  smaller.trees <- Tree1.2 %>% filter(INVYR > 2010, PREVDIA < DBH.CUTPOINT) %>%
    select(SUBP, TREE, STATECD, PLOT_FIADB) %>%
    mutate(smaller = 1)
  
  # Using "S" and "L" for Smaller and Larger than DBH CUTPOINT
  Tree1.2S <- Tree1.2 %>% 
    left_join(smaller.trees, by = c("SUBP", "TREE", "STATECD", "PLOT_FIADB")) %>%
    filter(smaller == 1) %>%
    select(-smaller)
  
  larger.trees <- Tree1.2 %>% filter(INVYR > 2010, PREVDIA >= DBH.CUTPOINT) %>%
    select(SUBP, TREE, STATECD, PLOT_FIADB) %>%
    mutate(larger = 1)
  
  Tree1.2L <- Tree1.2 %>% 
    left_join(larger.trees, by = c("SUBP", "TREE", "STATECD", "PLOT_FIADB")) %>%
    filter(larger == 1) %>%
    select(-larger)

  # nrow(Tree1.2) - (nrow(Tree1.2S) + nrow(Tree1.2L))   # Check to make sure no rows are missed.
  # 
}

PlotDat <- orig[, 1:5] %>% 
  left_join(fia.cond, by = c("PLOT_FIADB", "STATECD")) # Attaching SITECLCD for site productivity class code.

## REMOVE/UPDATE THIS ONCE I RECEIVE CHANGES FROM BRYCE ##
remove.plot.fiadb <- fia.cond %>% filter(is.na(SITECLCD_plot) == TRUE)
PlotDat <- PlotDat %>% anti_join(remove.plot.fiadb %>% select(PLOT_FIADB, STATECD), by = c("PLOT_FIADB", "STATECD"))
#nrow(remove.plot.fiadb) == nrow(orig) - nrow(PlotDat)  # Verify that the correct number of rows were removed. Should be "TRUE".

Tree1.3 <- left_join(Tree1.2, PlotDat %>% select(-STRATUM, -W_h), by = c("STATECD", "PLOT_FIADB", "State_Plot"))


#  All three analysis pathways create these files. The first batch of eight tree...dat.. 
#   objects will work for ANALYSIS.PATHWAYs 1 through 3.  "smaller" tells the function how to 
#   handle pathway #2 for the smaller
# Create a list of the number of dead trees of each species in each plot.

treedat.use <- if(ANALYSIS.PATHWAY == 2) Tree1.2S else Tree1.3

n_quant <- if(ANALYSIS.PATHWAY != 3) length(QUANT.LEVELS) else 7


#### 3) Specialized data prep --------------------------------------------------------
## First the data are summarized by species and climate variable (clim.mort.resp.fcn, 
#   clim.growth.resp.fcn).  Then the data are combined into arrays for 
#    analysis (parse.tree.clim.fcn).

tree.mort.dat.vpdmin1 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "vpdmin", treedat.sel = treedat.use) 
tree.mort.dat.vpdmax1 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "vpdmax", treedat.sel = treedat.use)
tree.mort.dat.temp1 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "temp", treedat.sel = treedat.use)
tree.mort.dat.precip1 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "precip", treedat.sel = treedat.use)

# Create lists of tree growth and number of trees that grew in each plot
tree.grow.dat.vpdmin1 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "vpdmin", treedat.sel = treedat.use) 
tree.grow.dat.vpdmax1 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "vpdmax", treedat.sel = treedat.use)
tree.grow.dat.temp1 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "temp", treedat.sel = treedat.use)
tree.grow.dat.precip1 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "precip", treedat.sel = treedat.use)

# Analysis pathway # 2 additionally creates these for larger trees
if(ANALYSIS.PATHWAY == 2) {
  
  # Create a list of the number of dead trees of each species in each plot.
  tree.mort.dat.vpdmin2 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "vpdmin", treedat.sel = Tree1.2L) 
  tree.mort.dat.vpdmax2 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "vpdmax", treedat.sel = Tree1.2L)
  tree.mort.dat.temp2 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "temp", treedat.sel = Tree1.2L)
  tree.mort.dat.precip2 <- map(sel.spp, clim.mort.resp.fcn, clim.var = "precip", treedat.sel = Tree1.2L)
  
  
  # Create lists of tree growth and number of trees that grew in each plot
  tree.grow.dat.vpdmin2 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "vpdmin", treedat.sel = Tree1.2L) 
  tree.grow.dat.vpdmax2 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "vpdmax", treedat.sel = Tree1.2L)
  tree.grow.dat.temp2 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "temp", treedat.sel = Tree1.2L)
  tree.grow.dat.precip2 <- map(sel.spp, clim.growth.resp.fcn, clim.var = "precip", treedat.sel = Tree1.2L)
}

# Combining data into arrays and such, preparing for analysis.
arrays.vpdmin.mort1 <- parse.tree.clim.fcn(tree.mort.dat.vpdmin1, "vpdmin", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
arrays.vpdmax.mort1 <- parse.tree.clim.fcn(tree.mort.dat.vpdmax1, "vpdmax", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
arrays.temp.mort1 <- parse.tree.clim.fcn(tree.mort.dat.temp1, "temp", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
arrays.precip.mort1 <- parse.tree.clim.fcn(tree.mort.dat.precip1, "precip", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)


arrays.vpdmin.grow1 <- parse.tree.clim.fcn(tree.grow.dat.vpdmin1, "vpdmin", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
arrays.vpdmax.grow1 <- parse.tree.clim.fcn(tree.grow.dat.vpdmax1, "vpdmax", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
arrays.temp.grow1 <- parse.tree.clim.fcn(tree.grow.dat.temp1, "temp", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
arrays.precip.grow1 <- parse.tree.clim.fcn(tree.grow.dat.precip1, "precip", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)

if(ANALYSIS.PATHWAY == 2) {
  arrays.vpdmin.mort2 <- parse.tree.clim.fcn(tree.mort.dat.vpdmin2, "vpdmin", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
  arrays.vpdmax.mort2 <- parse.tree.clim.fcn(tree.mort.dat.vpdmax2, "vpdmax", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
  arrays.temp.mort2 <- parse.tree.clim.fcn(tree.mort.dat.temp2, "temp", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
  arrays.precip.mort2 <- parse.tree.clim.fcn(tree.mort.dat.precip2, "precip", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP)
  
  
  arrays.vpdmin.grow2 <- parse.tree.clim.fcn(tree.grow.dat.vpdmin2, "vpdmin", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
  arrays.vpdmax.grow2 <- parse.tree.clim.fcn(tree.grow.dat.vpdmax2, "vpdmax", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
  arrays.temp.grow2 <- parse.tree.clim.fcn(tree.grow.dat.temp2, "temp", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
  arrays.precip.grow2 <- parse.tree.clim.fcn(tree.grow.dat.precip2, "precip", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP)
}

## Fire mortality: what proportion died from fire?  This info is only generated for pathway 1 and mortality.
if (ANALYSIS.PATHWAY == 1 & CALC.FIREPROP == TRUE) {
  fire.frac.table.fcn(tablename = "Fire_Prop_vpdmin_mort.csv", tableloc = paste0(RESULTS1.LOC, "Mort_figs_vpdmin/"), treedat = treedat.use, parseddat = arrays.vpdmin.mort1)  
  fire.frac.table.fcn(tablename = "Fire_Prop_vpdmax_mort.csv", tableloc = paste0(RESULTS1.LOC, "Mort_figs_vpdmax/"), treedat = treedat.use, parseddat = arrays.vpdmax.mort1)  
  fire.frac.table.fcn(tablename = "Fire_Prop_temp_mort.csv", tableloc = paste0(RESULTS1.LOC, "Mort_figs_temp/"), treedat = treedat.use, parseddat = arrays.temp.mort1)  
  fire.frac.table.fcn(tablename = "Fire_Prop_precip_mort.csv", tableloc = paste0(RESULTS1.LOC, "Mort_figs_precip/"), treedat = treedat.use, parseddat = arrays.precip.mort1)  
}



#### 4) Analysis and plotting  --------------------------------------------------------

# If desired, the analysis will obtain mortality estimates for each of the 50 species
#  by state and overall. These values can be used to evaluate species for inclusion in the
#  main analysis
if(RUN.STATES == TRUE) {
  source(paste0(CODE.LOC, "Overall_Mort_Est.R"))
}


# The remaining code breaks up the analyses by type (mortality, growth) and climate variables (4).
# If ANALYSIS.PATHWAY is 2 then the process is run twice, once for each size class. 
#  The plotting functions are located at the end of the process, with each ANALYSIS.PATHWAY
#  having a distinct function.

dataset.number <- if(ANALYSIS.PATHWAY == 2) 2 else 1 # How many times to run through.

# For ANALYSIS.PATHWAY = 1 and ITER = 200, the new machine takes: 23.5 minutes
# For ANALYSIS.PATHWAY = 1 and ITER = 1000, the new machine takes: 42 minutes
# ANALYSIS.PATHWAY = 2 and ITER = 200 takes the new machine 28 minutes

tic()
for (i in 1:length(CLIM.VAR)){
  for(j in 1:length(ANALYSIS.TYPE)) {
    
    # Setting up data set outputs, one for ANALYSIS.PATHWAY = 1, 3,
    #  two for ANALYSIS.PATHWAY 2.
    
    results <- list(results1 = NULL, results2 = NULL)
    
    for(k in 1:dataset.number) {
      
      
      var1 <- paste0("pre.", CLIM.VAR[i])
      var.delt <- paste0("delt.", CLIM.VAR[i])
      
      # Need climate names for files and axes.
      clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE)
      
      var.label <- clim.names$label[clim.names$values == var1]
      var.filename <- clim.names$filename[clim.names$values == var1]
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
      
      
      # Questions about the plot data
      # 1) Are there intensification plots
      #sum(orig$intensification)  # Yes, 6327
      # 2) How many forested plots?
      #length(orig$propfor[orig$propfor > 0])  # 18672 
      # 3) Do our trees exclusively occur in forested plots?
      #orig2 <- orig %>% mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
      #  left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
      #orig2$trees <- apply(orig2[, 10:59], 1, sum)
      #summary(orig2$propfor[orig2$trees > 0])
      #orig2[orig2$propfor == 0 & orig2$trees > 0,] # 79 plots
      # It is evident that the proportion of forest/non-forest/missed areas of plots is not necessarily constant between the two visits. 
      # 4) We have climate data for 37,612 plots and plot data for 38656 plots.  Are we missing climate data for any plots with our species of interest?
      #dim(clim.dat) ; dim(orig)   # 37612 13 , 38656 59
      #orig2 <- orig %>% mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
      #left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
      #length(orig2$pre.temp[is.na(orig2$pre.temp) == FALSE])  # 37612; all climate variable plot values had matches
      #x1 <- orig2[is.na(orig2$pre.temp) == TRUE,  ]
      #sum(apply(x1[, 10:59], 1, sum))   # result = zero.  None of these plots contain our trees of interest.
      
      
      
      n.spp <- length(SEL.SPP)
      
      
      # Strata for bootstrap resampling procedure
      strata <- unique(vals_dat$STRATUM)
      strata.num <- data.frame(val = 1:nrow(vals_dat), stratum = vals_dat$STRATUM) 
      
      
      plan(multisession, workers = n.cores) # Setting up parallel computing
      
      
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
    map(SEL.SPP, pair.plts.fcn, quant.table = results$results1$quant.table, var.label = results$results1$var.label, 
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
}   # end i

toc()











