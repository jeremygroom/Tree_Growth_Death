

source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))


##### ---- Code ---- #####

# DATA PREP

# Climate variables and reduction to a single table
tmp <- read_csv(paste0(DATA.LOC, "tmp.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.temp = post.temp - pre.temp) 
precip <- read_csv(paste0(DATA.LOC, "precip.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.precip = post.precip - pre.precip)
vpdmin <- read_csv(paste0(DATA.LOC, "vpdmin.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.vpdmin = post.vpdmin - pre.vpdmin)
vpdmax <- read_csv(paste0(DATA.LOC, "vpdmax.20.10.csv")) %>% dplyr::select(-INVYR) %>% mutate(delt.vpdmax = post.vpdmax - pre.vpdmax)


# clim.dat = all climate data
clim.dat <- list(tmp, precip, vpdmin, vpdmax) %>% reduce(left_join, by = c("STATECD", "PLOT_FIADB")) %>% 
  mutate(State_Plot = paste0(PLOT_FIADB, STATECD),
         State_Plot = as.numeric(State_Plot)) %>% 
  dplyr::select(-STATECD, -PLOT_FIADB)



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

fia.tree <- fia.tables$tree_vals %>% filter(INVYR > 2010 & INVYR < 2020) %>%
  select(-INVYR) %>%                      # Need to leave in to double-match trees to harvest events.
  mutate(CN = as.numeric(CN), 
         join.test = 1)

fia.cond <- fia.tables$cond_vals %>% filter(INVYR < 2011) %>% 
  select(STATECD, PLOT, SITECLCD) %>%
  filter(is.na(SITECLCD) == FALSE) %>% # May as well remove these to avoid multiple join values where SITECLCD = NA and an integer for a single PLT_CN.
  group_by(STATECD, PLOT) %>%
  reframe(SITECLCD = first(SITECLCD))   # This probably needs to be changed (or will change with Bryce's CLCD values) as I'm arbitrarily using the first value when more than one values exist.
  

Tree1.1 <- Tree1 %>% left_join(fia.tree, by = c("STATECD", "PLOT_FIADB" = "PLOT", "TREE", "SUBP")) # CN (Better if this is left out, as we want )

# Check: Did the join work as intended? Which trees were harvested, using STATUSCD and AGENTCD?
#sum(Tree1.1$join.test, na.rm = TRUE) # All fia.tree entries were successfully joined with Tree1 trees.
#table(Tree1.1$STATUSCD, Tree1.1$AGENTCD)   # 32073 trees burned (12609), were harvested (15260), or both (4204)
#tree.a80 <- Tree1.1 %>% filter(AGENTCD == 80, STATUSCD != 3)   # 4,204 trees = AGENTCD == 80, STATUSCD == 2
#tree.s3 <- Tree1.1 %>% filter(AGENTCD != 80, STATUSCD == 3)   # No trees were coded as harvest without AGENTCD == 80

# Extracting the trees which burned or were cut or cleared.  63,472
#tree.cut.fire <- Tree1.1 %>% filter(AGENTCD %in% c(30, 80)) %>% select(State_Plot_Subp_Tree)  # extracting trees that were cut and/or burned
tree.cut <- Tree1.1 %>% filter(AGENTCD == 80) %>% select(State_Plot_Subp_Tree)   # Extracting only trees that were cut. Burned remain.

# Removing cut/cleared trees from the tree list.
#Tree1.2 <- Tree1.1 %>% anti_join(tree.cut.fire, by = "State_Plot_Subp_Tree")
Tree1.2 <- Tree1.1 %>% anti_join(tree.cut, by = "State_Plot_Subp_Tree")



PlotDat <- orig[, 1:5] %>% 
  left_join(fia.cond, by = c("PLOT_FIADB" = "PLOT", "STATECD")) # Attaching SITECLCD for site productivity class code.


dat.cent.start <- tibble(var.deltvar = factor(QUANT.LEVELS, level = QUANT.LEVELS))

# Create a list of the number of dead trees of each species in each plot.
tree.mort.dat.vpdmin <- map(spp.list, clim.mort.resp.fcn, clim.var = "vpdmin") 
tree.mort.dat.vpdmax <- map(spp.list, clim.mort.resp.fcn, clim.var = "vpdmax")
tree.mort.dat.temp <- map(spp.list, clim.mort.resp.fcn, clim.var = "temp")
tree.mort.dat.precip <- map(spp.list, clim.mort.resp.fcn, clim.var = "precip")


# Create lists of tree growth and number of trees that grew in each plot
tree.grow.dat.vpdmin <- map(spp.list, clim.growth.resp.fcn, clim.var = "vpdmin") 
tree.grow.dat.vpdmax <- map(spp.list, clim.growth.resp.fcn, clim.var = "vpdmax")
tree.grow.dat.temp <- map(spp.list, clim.growth.resp.fcn, clim.var = "temp")
tree.grow.dat.precip <- map(spp.list, clim.growth.resp.fcn, clim.var = "precip")


parse.tree.clim.fcn <- function(tree.dat, clim.var, analysis.type, resp.dat, tot.dat) {  # tree.dat = prepped list of data, clim.var = which climate variable name, analysis.type = "mort" or "grow", resp.dat = response 
                                                                                        #   data set (i.e., "died.out" or "growth.val"), tot.dat = total trees data set (i.e., "all.trees" or "growth.n.trees")
  # Makes a list of all of the response tables and all-tree tables for the next step to operate on.
  extract.resp <- map(tree.dat , ~.[[resp.dat]]) 
  vals_dat <- reduce(extract.resp, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM","SITECLCD_plot")) # Combining the list into a single tibble.
  
  extract.all <- map(tree.dat , ~.[[tot.dat]]) 
  all_dat <- reduce(extract.all, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM", "SITECLCD_plot")) 
  
  
  # Now extracting the quantile category information for each species 
  extract.quant <- map(tree.dat , ~.[["quant.out"]]) 
  # Creating a 3D array out of the quantile info (plots x quantiles x spp)
  quant.array <- array(unlist(extract.quant), dim = c(nrow(extract.quant[[1]]), ncol(extract.quant[[1]]), length(extract.quant)), 
                       dimnames = list(c(1:nrow(extract.quant[[1]])),
                                       QUANT.LEVELS, 
                                       SEL.SPP))
  
  # Matrix of quantile values for each species' plots joined with climate data
  # The values for plots with trees of a species are given the ID number for their respective quantile (e.g., i = 9 has plot values of 0 or 9).
  # Note that this is just for the construction of the quantile matrix (values 0 through 9), not quant.array.
  q2 <- quant.array
  for (i in 1:n_quant) {
    q2[, i, ] <- q2[, i, ] * i
  }
  quant.matrix <- apply(q2, c(1, 3), sum) %>%  # collapsing the quantile values into a single matrix
    bind_cols(vals_dat %>% dplyr::select(State_Plot)) %>% # It does not matter whether vals_dat or all_dat, just want State_Plot
    left_join(clim.dat, by =  "State_Plot")
  
  # Number of plots in each quantile
  quant.n <- map(1:length(extract.quant), function(x) apply(quant.array[, , x], 2, sum))
  quant.n <- bind_cols(quant.n)
  names(quant.n) <- SEL.SPP
  
  extract.quant.lims <- map(tree.dat, ~.[["quant.lims"]]) 
  names(extract.quant.lims) <- SEL.SPP
  
  extract.delt.quant.lims <- map(tree.dat, ~.[["quant.lims.delta"]]) 
  names(extract.delt.quant.lims) <- SEL.SPP
  
  
  # Now extracting the centroid information for each quantile/species)
  extract.centroid <- map(tree.dat, ~.[["cent.loc"]]) 
  centroid.array <- array(unlist(extract.centroid), dim = c(nrow(extract.centroid[[1]]), ncol(extract.centroid[[1]]), length(extract.centroid)), 
                          dimnames = list(QUANT.LEVELS, 
                                          c("var.deltvar", "end.x", "end.y"),
                                          SEL.SPP))
  
  
  ### State level ###
  state.matrix <- quant.matrix[, 1:length(SEL.SPP)]
  state.matrix <- (state.matrix / state.matrix)  # obtain 1s and nan's
  state.matrix <- apply(state.matrix, 2, function(x) {ifelse(is.nan(x) == TRUE, 0, x)}) # Change the nan to 0
  all.array <- state.matrix
  state.matrix <- state.matrix %>% cbind(vals_dat[, 1])
  
  state.list <- unique(vals_dat$STATECD) # 6 = California, 41 = Oregon, 53 = Washington
  
  state.array.fcn <- function(sel.state, q.matrix) {
    q.matrix2 <- q.matrix %>% mutate(STATECD = ifelse(STATECD == sel.state, 1, 0))
    q.matrix2[, 1:length(SEL.SPP)] <- q.matrix2[, 1:length(SEL.SPP)] * q.matrix2[, ncol(q.matrix2)]
    q.matrix2 <- q.matrix2[, 1:length(SEL.SPP)]
    return(q.matrix2)
  }
  
  extract.state <- map(state.list, state.array.fcn, state.matrix)
  # Creating a 3D array out of the quantile info (plots x quantiles x spp)
  state.array <- array(unlist(extract.state), dim = c(nrow(extract.state[[1]]), ncol(extract.state[[1]]), length(state.list)),
                       dimnames = list(c(1:nrow(extract.state[[1]])),
                                       SEL.SPP,
                                       state.list))
  state.array <- aperm(state.array, c(1, 3, 2))  # Changing the order of the dimensions to match those for quantiles. 
  
  # Obtaining the number of plots per state.
  state.n.extract <- map(1:3, function(x) apply(extract.state[[x]], 2, sum)) #sum plots by state.
  state.n <- reduce(state.n.extract, bind_rows) %>%  # Combine state data into a tibble.
    mutate(state = state.list)
  
  ## We now have the data at the level needed to conduct the mortality analysis.  We need to select a species and a quantile and find the 
  #  estimated mortality rate for that quantile.  The quantile array (quant.array) offers a mask to remove all values other than the quantile 
  #  plots of interest for the species. 
  
  
  mort.grow.dat <- list(
    spp.list = tree.dat[[1]]$spp.list,
    spp.id = tree.dat[[1]]$spp.id,
    vals_dat = vals_dat, 
    all_dat = all_dat, 
    quant.array = quant.array, 
    state.array = state.array,
    all.array = all.array,
    quant.matrix = quant.matrix,
    centroid.array = centroid.array,
    quant.n = quant.n,
    state.n = state.n,
    quant.lims = extract.quant.lims,
    quant.lims.delt = extract.delt.quant.lims,
    tree.dat = Tree1.2)

## --- NOTE: This function used to save each file.  That's mothballed.  
#  rds.name <- paste0(DATA.LOC, analysis.type, "_dat_", clim.var, ".rds")
#  write_rds(mort.grow.dat, rds.name)
# zip(zipfile = paste0(DATA.LOC, analysis.type, "_dat_", clim.var, ".zip"), files = rds.name)
  
  # Deleting written RDS file (~150 MB)
#  if (file.exists(rds.name) == TRUE) {
#    file.remove(rds.name)
# }
}

parse.tree.clim.fcn(tree.mort.dat.vpdmin, "vpdmin", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees")
parse.tree.clim.fcn(tree.mort.dat.vpdmax, "vpdmax", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees")
parse.tree.clim.fcn(tree.mort.dat.temp, "temp", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees")
parse.tree.clim.fcn(tree.mort.dat.precip, "precip", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees")


parse.tree.clim.fcn(tree.grow.dat.vpdmin, "vpdmin", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees")
parse.tree.clim.fcn(tree.grow.dat.vpdmax, "vpdmax", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees")
parse.tree.clim.fcn(tree.grow.dat.temp, "temp", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees")
parse.tree.clim.fcn(tree.grow.dat.precip, "precip", analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees")



