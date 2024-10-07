

source("Global.R")


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

Tree1.1 <- Tree1 %>% left_join(fia.tree, by = c("STATECD", "PLOT_FIADB" = "PLOT", "TREE", "SUBP"))  # CN (Better if this is left out, as we want )

# Check: Did the join work as intended? Which trees were harvested, using STATUSCD and AGENTCD?
#sum(Tree1.1$join.test, na.rm = TRUE) # All fia.tree entries were successfully joined with Tree1 trees.
#table(Tree1.1$STATUSCD, Tree1.1$AGENTCD)   # 32073 trees burned (12609), were harvested (15260), or both (4204)
#tree.a80 <- Tree1.1 %>% filter(AGENTCD == 80, STATUSCD != 3)   # 4,204 trees = AGENTCD == 80, STATUSCD == 2
#tree.s3 <- Tree1.1 %>% filter(AGENTCD != 80, STATUSCD == 3)   # No trees were coded as harvest without AGENTCD == 80

# Extracting the trees which burned or were cut or cleared.  63,472
tree.cut.fire <- Tree1.1 %>% filter(AGENTCD %in% c(30, 80)) %>% select(State_Plot_Subp_Tree)

# Removing cut/cleared trees from the tree list.
Tree1.2 <- Tree1.1 %>% anti_join(tree.cut.fire, by = "State_Plot_Subp_Tree")



PlotDat <- orig[, 1:5]
dat.cent.start <- tibble(var.deltvar = factor(QUANT.LEVELS, level = QUANT.LEVELS))



## Function to determine the number of trees that died in a plot and their climate quantiles.
clim.tree.resp.fcn <- function(spcd, clim.var) {
  var1 <- paste0("pre.", clim.var)
  var.delt <- paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  Tree2 <- Tree1.2 %>% filter(SPCD == spcd) %>%
    dplyr::select(State_Plot, SUBP, TREE, SDN)
  
  tree.plots <- clim.dat %>% filter(State_Plot %in% Tree2$State_Plot) %>% 
    dplyr::select(State_Plot, all_of(var1), all_of(var.delt))

  # Data quantiles
  quant.lims <- apply(tree.plots[, c(which(names(tree.plots) == var1), which(names(tree.plots) == var.delt))], 2, quantile, probs = QUANT.PROBS)
  
  # Applying quantiles to a particular species
  tree.plots2 <- tree.plots %>% dplyr::select(State_Plot, all_of(c(var1, var.delt))) %>% 
    mutate(var.quants = ifelse(get(var1) > quant.lims[2, 1], 1, ifelse(get(var1) < quant.lims[1, 1], -1, 0)),
           deltvar.quants = ifelse(get(var.delt) > quant.lims[2, 2], 1, ifelse(get(var.delt) < quant.lims[1, 2], -1, 0)),
           vq = case_match(var.quants, -1 ~ "Li", 0 ~ "Mi", 1 ~ "Hi"),
           dvq = case_match(deltvar.quants, -1 ~ "Lc", 0 ~ "Mc", 1 ~ "Hc"),
           var.deltvar = factor(paste0(vq, dvq), levels = c("LiHc", "MiHc", "HiHc", "LiMc", "MiMc", "HiMc", "LiLc", "MiLc", "HiLc")))
  
    
  # For centroid calcluation
  dat.cent <- tree.plots2 %>% group_by(var.deltvar) %>% # Quartile centers, good for finding centroids in each quadrant
    reframe(end.x = mean(get(var1)),        # Referred to as "end" because of distance measurements from the main centroid to each of the outer (end) centroids.
            end.y = mean(get(var.delt))) %>% 
    right_join(dat.cent.start, by = "var.deltvar") %>%
    arrange(var.deltvar)
  
  # Joining tree data with quantile assignment, providing plot-level summary of the number of trees, number that died, and percent that died.
  Tree3 <- left_join(Tree2, tree.plots2, by = "State_Plot") %>%
    group_by(State_Plot, var.deltvar, get(var1), get(var.delt)) %>%  # This will assign the names `get(var1)` and `get(var.delt)`, which are cleaned up below
    reframe(n.trees = n(),
            n.died = length(SDN[SDN == 2]),
            pct.died = n.died / n.trees) %>%
    rename(!!var1 := `get(var1)`,
           !!var.delt := `get(var.delt)`) 
           
  # Separate out a data set for those that died
  died.out <- left_join(PlotDat, Tree3 %>% dplyr::select(State_Plot, n.died), by = "State_Plot") %>%
    rename(!!return.col := n.died)
  died.out[is.na(died.out)] <- 0
  
  # Data set for all trees, alive and died
  all.trees <- left_join(PlotDat, Tree3 %>% dplyr::select(State_Plot, n.trees), by = "State_Plot") %>%
    rename(!!return.col := n.trees)
  all.trees[is.na(all.trees)] <- 0

  ## Quantile masks for each species
  quant.divide.fcn <- function(quant.level) {
    quant.name <- paste0(quant.level, ".",  spcd)
    quant.died.out <- left_join(PlotDat, tree.plots2 %>% 
                                  dplyr::select(State_Plot, var.deltvar) %>%
                                  filter(var.deltvar == quant.level), by = "State_Plot") %>%
      mutate(var.deltvar = ifelse(is.na(var.deltvar), 0, 1)) %>%
      rename(!!quant.name := var.deltvar) 
  }
  
  quant.out1 <- map(QUANT.LEVELS, quant.divide.fcn)
  quant.out2 <- reduce(quant.out1, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM")) %>%
    dplyr::select(-c("STATECD", "PLOT_FIADB", "State_Plot", "STRATUM", "W_h"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, died.out = died.out, all.trees = all.trees, quant.out = quant.out2, quant.lims = quant.lims, cent.loc = dat.cent)
  
  return(out.list)
  
}


# Create a list of the number of dead trees of each species in each plot.
tree.dat <- map(spp.list, clim.tree.resp.fcn, clim.var = "vpdmin") 

# Makes a list of all of the "(all)died.out" tables for the next step to operate on.
extract.diedout <- map(tree.dat , ~.[["died.out"]]) 
ado_vals <- reduce(extract.diedout, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM")) # Combining the list into a single tibble.

extract.all <- map(tree.dat , ~.[["all.trees"]]) 
all_vals <- reduce(extract.all, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM")) 

# Now extracting the quantile category information for each species 
extract.quant <- map(tree.dat , ~.[["quant.out"]]) 
  # Creating a 3D array out of the quantile info (plots x quantiles x spp)
quant.array <- array(unlist(extract.quant), dim = c(nrow(extract.quant[[1]]), ncol(extract.quant[[1]]), length(extract.quant)), 
                     dimnames = list(c(1:nrow(extract.quant[[1]])),
                                     QUANT.LEVELS, 
                                     spp.id))

 # Matrix of quantile values for each species' plots joined with climate data
# The values for plots with trees of a species are given the ID number for their respective quantile (e.g., i = 9 has plot values of 0 or 9).
# Note that this is just for the construction of the quantile matrix (values 0 through 9), not quant.array.
q2 <- quant.array
for (i in 1:n_quant) {
  q2[, i, ] <- q2[, i, ] * i
}
quant.matrix <- apply(q2, c(1, 3), sum) %>%  # collapsing the quantile values into a single matrix
  bind_cols(ado_vals %>% dplyr::select(State_Plot)) %>% 
  left_join(clim.dat, by =  "State_Plot")

  # Number of plots in each quantile
quant.n <- map(1:length(extract.quant), function(x) apply(quant.array[, , x], 2, sum))
quant.n <- bind_cols(quant.n)
names(quant.n) <- spp.id

extract.quant.lims <- map(tree.dat, ~.[["quant.lims"]]) 
names(extract.quant.lims) <- spp.id

# Now extracting the centroid information for each quantile/species)
extract.centroid <- map(tree.dat, ~.[["cent.loc"]]) 
centroid.array <- array(unlist(extract.centroid), dim = c(nrow(extract.centroid[[1]]), ncol(extract.centroid[[1]]), length(extract.centroid)), 
                     dimnames = list(QUANT.LEVELS, 
                                       c("var.deltvar", "end.x", "end.y"),
                                       spp.id))


### State level ###
state.matrix <- quant.matrix[, 1:length(spp.list)]
state.matrix <- (state.matrix / state.matrix)
state.matrix <- apply(state.matrix, 2, function(x) {ifelse(is.nan(x) == TRUE, 0, x)})
all.array <- state.matrix
state.matrix <- state.matrix %>% cbind(ado_vals[, 1])

state.list <- unique(ado_vals$STATECD) # 6 = California, 41 = Oregon, 53 = Washington

state.array.fcn <- function(sel.state, q.matrix) {
  q.matrix2 <- q.matrix %>% mutate(STATECD = ifelse(STATECD == sel.state, 1, 0))
  q.matrix2[, 1:length(spp.list)] <- q.matrix2[, 1:length(spp.list)] * q.matrix2[, ncol(q.matrix2)]
  q.matrix2 <- q.matrix2[, 1:length(spp.list)]
  return(q.matrix2)
}

extract.state <- map(state.list, state.array.fcn, state.matrix)
# Creating a 3D array out of the quantile info (plots x quantiles x spp)
state.array <- array(unlist(extract.state), dim = c(nrow(extract.state[[1]]), ncol(extract.state[[1]]), length(state.list)),
                     dimnames = list(c(1:nrow(extract.state[[1]])),
                                     spp.id,
                                     state.list))
state.array <- aperm(state.array, c(1, 3, 2))  # Changing the order of the dimensions to match those for quantiles. 

# Obtaining the number of plots per state.
state.n.extract <- map(1:3, function(x) apply(extract.state[[x]], 2, sum)) #sum plots by state.
state.n <- reduce(state.n.extract, bind_rows) %>%  # Combine state data into a tibble.
  mutate(state = state.list)

## We now have the data at the level needed to conduct the mortality analysis.  We need to select a species and a quantile and find the 
#  estimated mortality rate for that quantile.  The quantile array (quant.array) offers a mask to remove all values other than the quantile 
#  plots of interest for the species. 


dieout.dat <- list(
                   spp.list = tree.dat[[1]]$spp.list,
                   spp.id = tree.dat[[1]]$spp.id,
                   ado_vals = ado_vals, 
                   all_vals = all_vals, 
                   quant.array = quant.array, 
                   state.array = state.array,
                   all.array = all.array,
                   quant.matrix = quant.matrix,
                   centroid.array = centroid.array,
                   quant.n = quant.n,
                   state.n = state.n,
                   quant.lims = extract.quant.lims,
                   tree.dat = Tree1.2)
write_rds(dieout.dat, paste0(DATA.LOC, "dieout_dat.rds"))

zip(zipfile = paste0(DATA.LOC, "dieout_dat.zip"), files = paste0(DATA.LOC, "dieout_dat.rds"))

# Deleting written RDS file (~150 MB)
if (file.exists(paste0(DATA.LOC, "dieout_dat.rds")) == TRUE) {
  file.remove(paste0(DATA.LOC, "dieout_dat.rds"))
}

