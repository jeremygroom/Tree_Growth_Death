
##### ---- Functions ---- #####

### DATA PREP FUNCTIONS - MORTALITY ####


## Quantile masks for each species. Used in clim.tree.resp.fcn (called in Data_Prep.R), operates on a single quantile level and species 
quant.divide.fcn <- function(quant.level, spcd2, tree.plots2, PlotDat) {
  quant.name <- paste0(quant.level, ".",  spcd2)
  quant.died.out <- left_join(PlotDat, tree.plots2 %>% 
                                dplyr::select(puid, var.deltvar) %>% # var.deltvar defined in clim.tree.resp.fcn, "LiHc" etc.
                                filter(var.deltvar == quant.level), by = "puid") %>%
    mutate(var.deltvar = ifelse(is.na(var.deltvar), 0, 1)) %>%
    rename(!!quant.name := var.deltvar) 
}



#### Data prep function for mortality --------------------------------------------------

## Function to determine the number of trees that died in a plot and their climate quantiles, used in Data_Prep.R.
clim.mort.resp.fcn <- function(spcd, clim.var, treedat.sel, clim.dat) {
  var1 <- "pre_mean" #paste0("pre.", clim.var)
  var.delt <- "difference"  #paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
# For alive/dead tree assessments
  Tree2 <- treedat.sel %>% 
    filter(SPCD == spcd) %>%
    dplyr::select(puid, REMPER, TREE, STATUSCD, ntree.rep) %>% distinct()
  
  tree.plots <- clim.dat %>% filter(puid %in% Tree2$puid) %>% 
    dplyr::select(puid, all_of(var1), all_of(var.delt))
  
  # Data quantiles
  # Setting the delta-var amounts to be fixed values.  See Global for setting VAR.DELTA.BOUNDARIES.
  quant.lims.delta.pos <- VAR.DELTA.BOUNDARIES$max.min[VAR.DELTA.BOUNDARIES$clim.var == clim.var] #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- -quant.lims.delta.pos #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
  quant.lims.delta <- c(quant.lims.delta.neg, quant.lims.delta.pos)
  quant.lims <- quantile(get(var1, tree.plots), probs = QUANT.PROBS)
  
  # Applying quantiles to a particular species
  if(ANALYSIS.PATHWAY != 3) {
    tree.plots2 <- tree.plots %>% dplyr::select(puid, all_of(c(var1, var.delt))) %>% 
      mutate(var.quants = ifelse(get(var1) > quant.lims[2], 1, ifelse(get(var1) < quant.lims[1], -1, 0)),
             deltvar.quants = ifelse(get(var.delt) > quant.lims.delta[2], 1, ifelse(get(var.delt) < quant.lims.delta[1], -1, 0)),
             vq = case_match(var.quants, -1 ~ "Li", 0 ~ "Mi", 1 ~ "Hi"),
             dvq = case_match(deltvar.quants, -1 ~ "Lc", 0 ~ "Mc", 1 ~ "Hc"),
             var.deltvar = factor(paste0(vq, dvq), levels = c("LiHc", "MiHc", "HiHc", "LiMc", "MiMc", "HiMc", "LiLc", "MiLc", "HiLc")))
    
    # For centroid calcluation
    dat.cent <- tree.plots2 %>% group_by(var.deltvar) %>% # Quartile centers, good for finding centroids in each quadrant
      reframe(end.x = mean(get(var1)),        # Referred to as "end" because of distance measurements from the main centroid to each of the outer (end) centroids.
              end.y = mean(get(var.delt))) %>% 
      right_join(quant.level.table %>% select(-q.num), by = "var.deltvar") %>%
      arrange(var.deltvar)
    
  } else {
    # For Analysis Pathway 3.  A different set-up - not using the 9 quantiles, but in effect 7 (for Site Classes). The plots will be much simpler.
    tree.plots2 <- tree.plots %>% left_join(PlotDat %>% select(puid, SITECLCD_plot), by = "puid") %>% # Adding on the Site Class variable
      filter(is.na(SITECLCD_plot) == FALSE) %>%   # Removing any NA's in case they snuck in there.
      mutate(var.deltvar = factor(paste0("Class.", SITECLCD_plot)))
    
    dat.cent <- NULL # Don't need those.
    
    QUANT.LEVELS <- levels(tree.plots2$var.deltvar) # Renaming the QUANT.LEVELS away from "LiHc" and such
  }
  
  
  
  # Joining tree data with quantile assignment, providing plot-level summary of the number of trees, number that died, and percent that died.
  
  # Mortality data
  Tree3 <- left_join(Tree2, tree.plots2, by = "puid") %>%
    group_by(puid) %>% #, var.deltvar, get(var1), get(var.delt)) %>%  # This will assign the names `get(var1)` and `get(var.delt)`, which are cleaned up below
    reframe(n.trees = sum(ntree.rep),
            n.died.yr = sum(ntree.rep[STATUSCD == 2])/mean(REMPER),
            pct.died.yr = n.died.yr / n.trees)
  
  # Separate out a data set for those that died
  died.out <- left_join(PlotDat, Tree3 %>% dplyr::select(puid, n.died.yr), by = "puid") %>%
    rename(!!return.col := n.died.yr)
  died.out[is.na(died.out)] <- 0
  
  # Data set for all trees, alive and died
  all.trees <- left_join(PlotDat, Tree3 %>% dplyr::select(puid, n.trees), by = "puid") %>%
    rename(!!return.col := n.trees)
  all.trees[is.na(all.trees)] <- 0
  
  
  quant.out1 <- map(QUANT.LEVELS, quant.divide.fcn, spcd2 = spcd, tree.plots2 = tree.plots2, PlotDat = PlotDat)
  quant.out2 <- reduce(quant.out1, left_join, by = c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts")) %>%
    dplyr::select(-c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, died.out = died.out, all.trees = all.trees, 
                   quant.out = quant.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}


#### Data prep function for growth  --------------------------------------------------

## Function to determine the number of trees that died in a plot and their climate quantiles, used in Data_Prep.R.
clim.growth.resp.fcn <- function(spcd, clim.var, treedat.sel, clim.dat) {
  var1 <- "pre_mean" #paste0("pre.", clim.var)
  var.delt <- "difference"  #paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  Tree.spp <- treedat.sel %>% filter(SPCD == spcd, STATUSCD == 1)
  
  # For basal area growth/tree
  TreeG <- Tree.spp %>%
    # Amount of tree growth per year since the previous visit, multiplied by the number of trees per acre the tree represents. 
    mutate(BA_Growth = ntree.rep * (pi * (DIA / 2)^2 - pi * (PREVDIA / 2)^2)/REMPER) %>%
    group_by(puid) %>%
    #Goal: sum(change in BA/REMPER)/Total # trees
    reframe(n_g.trees = sum(ntree.rep),
            total.growth = sum(BA_Growth))

  tree.plots <- Tree.spp %>% 
    dplyr::select(puid, all_of(var1), all_of(var.delt)) %>% distinct()
  
  # Data quantiles
  quant.lims.delta.pos <- VAR.DELTA.BOUNDARIES$max.min[VAR.DELTA.BOUNDARIES$clim.var == clim.var] #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- -quant.lims.delta.pos #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
  quant.lims.delta <- c(quant.lims.delta.neg, quant.lims.delta.pos)
  quant.lims <- quantile(get(var1, tree.plots), probs = QUANT.PROBS)
  
  # Applying quantiles to a particular species
  if(ANALYSIS.PATHWAY != 3) {
    tree.plots2 <- tree.plots %>% dplyr::select(puid, all_of(c(var1, var.delt))) %>% 
      mutate(var.quants = ifelse(get(var1) > quant.lims[2], 1, ifelse(get(var1) < quant.lims[1], -1, 0)),
             deltvar.quants = ifelse(get(var.delt) > quant.lims.delta[2], 1, ifelse(get(var.delt) < quant.lims.delta[1], -1, 0)),
             vq = case_match(var.quants, -1 ~ "Li", 0 ~ "Mi", 1 ~ "Hi"),
             dvq = case_match(deltvar.quants, -1 ~ "Lc", 0 ~ "Mc", 1 ~ "Hc"),
             var.deltvar = factor(paste0(vq, dvq), levels = c("LiHc", "MiHc", "HiHc", "LiMc", "MiMc", "HiMc", "LiLc", "MiLc", "HiLc")))
    
    # For centroid calcluation
    dat.cent <- tree.plots2 %>% group_by(var.deltvar) %>% # Quartile centers, good for finding centroids in each quadrant
      reframe(end.x = mean(get(var1)),        # Referred to as "end" because of distance measurements from the main centroid to each of the outer (end) centroids.
              end.y = mean(get(var.delt))) %>% 
      right_join(quant.level.table %>% select(-q.num), by = "var.deltvar") %>%
      arrange(var.deltvar)
    
  } else {
    # For Analysis Pathway 3.  A different set-up - not using the 9 quantiles, but in effect 7 (for Site Classes). The plots will be much simpler.
    tree.plots2 <- tree.plots %>% left_join(PlotDat %>% select(puid, SITECLCD_plot), by = "puid") %>% # Adding on the Site Class variable
      filter(is.na(SITECLCD_plot) == FALSE) %>%   # Removing any NA's in case they snuck in there.
      mutate(var.deltvar = factor(paste0("Class.", SITECLCD_plot)))
    
    dat.cent <- NULL
    
    QUANT.LEVELS <- levels(tree.plots2$var.deltvar) # Renaming the QUANT.LEVELS away from "LiHc" and such
  }
  
  
  
  # Joining tree data with quantile assignment, providing plot-level summary of the number of trees, number that died, and percent that died.
  # Growth data  
  # Overall data:
  Growth.dat <- left_join(PlotDat, TreeG, by = "puid") 
  Growth.dat[is.na(Growth.dat)] <- 0
  
  # Cubic inches of growth by plot
  growth.val <- Growth.dat %>% 
    dplyr::select(-n_g.trees) %>%                 
    rename(!!return.col := total.growth) 
  
  # Number of 'growth trees' by plot
  growth.n.trees <- Growth.dat %>% 
    dplyr::select(-total.growth) %>%                 
    rename(!!return.col := n_g.trees) 
  
  
  quant.out1 <- map(QUANT.LEVELS, quant.divide.fcn, spcd2 = spcd, tree.plots2 = tree.plots2, PlotDat = PlotDat)
  quant.out2 <- reduce(quant.out1, left_join, by = c("STATECD", "puid", "n_h.plts", "ESTN_UNIT", "STRATUMCD", "w", "stratum")) %>%
    dplyr::select(-c("STATECD", "puid", "n_h.plts", "ESTN_UNIT", "STRATUMCD", "w", "stratum"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, growth.val = growth.val, growth.n.trees = growth.n.trees,
                   quant.out = quant.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}


#### Preparing mortality/growth output for analysis -------------------------------


## This is a large function that takes species-specific data and creates arrays to assist with analyses.
parse.tree.clim.fcn <- function(tree.dat, clim.var, analysis.type, resp.dat, tot.dat, selected.spp, clim.dat) {  
  # tree.dat = prepped list of data, clim.var = which climate variable name, 
  #  analysis.type = "mort" or "grow", resp.dat = response data set (i.e., 
  # "died.out" or "growth.val"), tot.dat = total trees data set (i.e., "all.trees" or "growth.n.trees")
  # Makes a list of all of the response tables and all-tree tables for the next step to operate on.
  
  extract.resp <- map(tree.dat , ~.[[resp.dat]]) 
  vals_dat <- reduce(extract.resp, left_join, by = c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts")) # Combining the list into a single tibble.
  #vals_dat <- reduce(extract.resp, left_join, by = c("STATECD", "PLOT_FIADB", "puid", "W_h", "STRATUM","SITECLCD_plot")) # Combining the list into a single tibble.
  
  extract.all <- map(tree.dat , ~.[[tot.dat]]) 
  all_dat <- reduce(extract.all, left_join, by = c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts")) 
  
  
  # Now extracting the quantile category information for each species 
  extract.quant <- map(tree.dat , ~.[["quant.out"]]) 
  # Creating a 3D array out of the quantile info (plots x quantiles x spp)
  quant.array <- array(unlist(extract.quant), dim = c(nrow(extract.quant[[1]]), ncol(extract.quant[[1]]), length(extract.quant)), 
                       dimnames = list(c(1:nrow(extract.quant[[1]])),
                                       QUANT.LEVELS, 
                                       selected.spp))
  
  # Matrix of quantile values for each species' plots joined with climate data
  # The values for plots with trees of a species are given the ID number for their respective quantile (e.g., i = 9 has plot values of 0 or 9).
  # Note that this is just for the construction of the quantile matrix (values 0 through 9), not quant.array.
  q2 <- quant.array
  for (i in 1:n_quant) {
    q2[, i, ] <- q2[, i, ] * i
  }
  quant.matrix <- apply(q2, c(1, 3), sum) %>%  # collapsing the quantile values into a single matrix
    bind_cols(vals_dat %>% dplyr::select(puid)) %>% # It does not matter whether vals_dat or all_dat, just want puid
    left_join(clim.dat, by =  "puid")
  
  # Number of plots in each quantile
  quant.n <- map(1:length(extract.quant), function(x) apply(quant.array[, , x], 2, sum))
  quant.n <- bind_cols(quant.n)
  names(quant.n) <- selected.spp
  
  extract.quant.lims <- map(tree.dat, ~.[["quant.lims"]]) 
  names(extract.quant.lims) <- selected.spp
  
  extract.delt.quant.lims <- map(tree.dat, ~.[["quant.lims.delta"]]) 
  names(extract.delt.quant.lims) <- selected.spp
  
  
  # Now extracting the centroid information for each quantile/species)
  if(ANALYSIS.PATHWAY != 3) {
    extract.centroid <- map(tree.dat, ~.[["cent.loc"]]) 
    centroid.array <- array(unlist(extract.centroid), dim = c(nrow(extract.centroid[[1]]), ncol(extract.centroid[[1]]), length(extract.centroid)), 
                            dimnames = list(QUANT.LEVELS, 
                                            c("var.deltvar", "end.x", "end.y"),
                                            selected.spp))
  } else {
    extract.centroid <- centroid.array <-  NULL
  }
  
  ### State level ###
  state.matrix <- quant.matrix[, 1:length(selected.spp)]
  state.matrix <- (state.matrix / state.matrix)  # obtain 1s and nan's
  state.matrix <- apply(state.matrix, 2, function(x) {ifelse(is.nan(x) == TRUE, 0, x)}) # Change the nan to 0
  state.matrix <- state.matrix %>% cbind(vals_dat[, 1])
  state.matrix <- data.frame(state.matrix)
  colnames(state.matrix)[ncol(state.matrix)] <- "STATECD"
  
  state.list <- unique(vals_dat$STATECD) # 6 = California, 41 = Oregon, 53 = Washington
  
  state.array.fcn <- function(sel.state, q.matrix) {
    q.matrix2 <- q.matrix %>% mutate(STATECD = ifelse(STATECD == sel.state, 1, 0))
    q.matrix2[, 1:length(selected.spp)] <- q.matrix2[, 1:length(selected.spp)] * q.matrix2[, ncol(q.matrix2)]
    q.matrix2 <- q.matrix2[, 1:length(selected.spp)]
    return(q.matrix2)
  }
  
  extract.state <- map(state.list, state.array.fcn, state.matrix)
  # Creating a 3D array out of the quantile info (plots x quantiles x spp)
  state.array <- array(unlist(extract.state), dim = c(nrow(extract.state[[1]]), ncol(extract.state[[1]]), length(state.list)),
                       dimnames = list(c(1:nrow(extract.state[[1]])),
                                       selected.spp,
                                       state.list))
  state.array <- aperm(state.array, c(1, 3, 2))  # Changing the order of the dimensions to match those for quantiles. 
  
  # Obtaining the number of plots per state.
  state.n.extract <- map(1:3, function(x) apply(extract.state[[x]], 2, sum)) #sum plots by state.
  state.n <- reduce(state.n.extract, bind_rows) %>%  # Combine state data into a tibble.
    mutate(state = state.list)
  
  
  
  ### All data level ###
  all.matrix <- quant.matrix[, 1:length(selected.spp)] 
  all.matrix <- apply(all.matrix, 2, function(x) ifelse(x > 0, 1, 0)) %>% 
    data.frame()
  all.array <- abind(all.matrix, all.matrix, along = 3)
  all.array <- aperm(all.array, c(1, 3, 2))  # Changing the order of the dimensions to match those for quantiles. 
  
  
  all.n <- apply(all.array[, 1, ], 2, sum)
  all.n <- tibble(n = all.n, labs = selected.spp) %>% 
    pivot_wider(names_from = "labs", values_from = "n") %>%
    mutate(all = 1)
  
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
    all.n = all.n,
    quant.lims = extract.quant.lims,
    quant.lims.delt = extract.delt.quant.lims)
  
  ## --- NOTE: This function used to save each file.  That's mothballed. --- ## 
  #  rds.name <- paste0(DATA.LOC, analysis.type, "_dat_", clim.var, ".rds")
  #  write_rds(mort.grow.dat, rds.name)
  # zip(zipfile = paste0(DATA.LOC, analysis.type, "_dat_", clim.var, ".zip"), files = rds.name)
  
  # Deleting written RDS file (~150 MB)
  #  if (file.exists(rds.name) == TRUE) {
  #    file.remove(rds.name)
  # }
}




###  Side project: determining the proportion of dead trees that died from fire-------------


 # This function is run by species in the fire.frac.table.fcn.  
fire.frac.dat.fcn <- function(tspp, treedat.sel, parsed.dat) {  # Tree species "X123", 
                                           # the tree data for analysis, and the files output
                                           # from parse.tree.clim.fcn.  
  x <- treedat.sel %>% filter(SPCD %in% as.numeric(gsub("X", "", tspp))) %>%
    select(puid, STATUSCD, AGENTCD, ntree.rep) %>%
    mutate(death_fire = ifelse(AGENTCD == 30, 1, 0),
           death_insect = ifelse(AGENTCD == 10, 1, 0),
           death_disease = ifelse(AGENTCD == 20, 1, 0),
           death_animal = ifelse(AGENTCD == 40, 1, 0),
           death_weather = ifelse(AGENTCD == 50, 1, 0),
           death_vegetation = ifelse(AGENTCD == 60, 1, 0),
           death_unknown = ifelse(AGENTCD == 70, 1, 0),
           death_silviculture = ifelse(AGENTCD == 80, 1, 0),
           death_all = ifelse(STATUSCD == 2, 1, 0),
           n_all = 1) %>%
    group_by(puid) %>%
    reframe(n_dead.fire = sum(death_fire * ntree.rep, na.rm = TRUE),
            n_dead.insect = sum(death_insect * ntree.rep, na.rm = TRUE),
            n_dead.disease = sum(death_disease * ntree.rep, na.rm = TRUE),
            n_dead.animal = sum(death_animal * ntree.rep, na.rm = TRUE),
            n_dead.weather = sum(death_weather * ntree.rep, na.rm = TRUE),
            n_dead.vegetation = sum(death_vegetation * ntree.rep, na.rm = TRUE),
            n_dead.unknown = sum(death_unknown * ntree.rep, na.rm = TRUE),
            n_dead.silviculture = sum(death_silviculture * ntree.rep, na.rm = TRUE),
            n_dead.all = sum(death_all * ntree.rep, na.rm = TRUE),
            n_all = sum(n_all * ntree.rep, na.rm = TRUE)) 
  
  y <- parsed.dat$quant.matrix %>% dplyr::select(all_of(tspp), puid) %>%
    left_join(x, by = "puid") %>%
    rename("Quants" := (!!tspp)) %>%
    filter(Quants > 0) %>%
    group_by(Quants) %>%
    reframe(tot = sum(n_all, na.rm = TRUE),
            tot.dead = sum(n_dead.all, na.rm = TRUE),
            tot.burn = sum(n_dead.fire, na.rm = TRUE),
            tot.insect = sum(n_dead.insect, na.rm = TRUE),
            tot.disease = sum(n_dead.disease, na.rm = TRUE),
            tot.animal = sum(n_dead.animal, na.rm = TRUE),
            tot.weather = sum(n_dead.weather, na.rm = TRUE),
            tot.vegetation = sum(n_dead.vegetation, na.rm = TRUE),
            tot.unknown = sum(n_dead.unknown, na.rm = TRUE),
            tot.silviculture = sum(n_dead.silviculture, na.rm = TRUE),
            ) %>%
    mutate(frac.died = tot/tot.dead,
           frac.burn = tot.burn/tot.dead,
           frac.insect = tot.insect/tot.dead,
           frac.disease = tot.disease/tot.dead,
           frac.animal = tot.animal/tot.dead,
           frac.weather = tot.weather/tot.dead,
           frac.vegetation = tot.vegetation/tot.dead,
           frac.unknown = tot.unknown/tot.dead,
           frac.silviculture = tot.silviculture/tot.dead,
           spp = tspp)
  
  return(y)
}

fire.frac.table.fcn <- function(tablename, tableloc, treedat, parseddat) {
  fire.frac.dat <- map(SEL.SPP, fire.frac.dat.fcn, treedat.sel  = treedat, parsed.dat = parseddat)
  fire.frac.dat2 <- rbindlist(fire.frac.dat) %>% data.frame()   # data.table function
  fwrite(fire.frac.dat2, file = paste0(tableloc, tablename))    # data.table csv write function
}



### ------ Mortality and growth estimation -------  ####

#  This function is nested within quant.est.spp.fcn, which is in turn nested within bs.fcn and finally q_mort.grow.fcn.  It is called in DeadTree_Analysis
# Find the estimates for mu_hat (the ratio estimator), Y_hat (numerator, thing we're interested in), and 
#    X_Hat (denominator, density of things like total trees of a species). See Equation 6 and related.
mean.q.fcn <- function(dat_tot, dat_vals, spp.sel1) {     #dat_tot = data frame of species total for a quantile (denominator). dat_vals = response data frame (numerator) for quantile.
  
  Zt <- Yv <- R <- 0 #  rep(0,length(sppname))    # Mean difference in response for a given species: Zt = Mean value for total trees, 
  #Yv = mean value for/of trees of interest (mean # that died), R = ratio of the two.
  #for (i in 1:length(sppname)) {                         # i indexes species
  col.spp <- grep(paste0("nd.", spp.sel1), colnames(dat_tot))
  col.name <- paste0("nd.", spp.sel1)
  
  # Old slow code.  For understanding the data.table usage.
#  Zt_i <- Yv_i <- rep(0, length(strata) )  # Weighted values for each strata
#  for (h in 1:strat.num) {                  # h indexes strata
#    w_h <- unique(dat_tot$w[dat_tot$stratum == strata[h]])
#    n_h <- unique(dat_tot$n_h.plts[dat_tot$stratum == strata[h]])    # Number of plots in stratum h 
#    Yv_i[h] <- w_h * sum(get(col.name, dat_vals)[dat_vals$stratum == strata[h]]) / n_h 
#    Zt_i[h] <- w_h * sum(get(col.name, dat_tot)[dat_tot$stratum == strata[h]]) / n_h     
#  }
#  Yv <- sum(Yv_i) 
#  Zt <- sum(Zt_i) 
#  R <- Yv / Zt

  # New data.table way
  dt_tot <- as.data.table(dat_tot)
  dt_vals <- as.data.table(dat_vals)
  
  tot_summary <- dt_tot[, .(w_h = w[1], n_h = n_h.plts[1], tot_sum = sum(get(col.name))), by = stratum]
  vals_summary <- dt_vals[, .(vals_sum = sum(get(col.name))), by = stratum]
  
  ests.out <- tot_summary[vals_summary, on = "stratum"]
  ests.out[, `:=`(Zt_i = w_h * tot_sum / n_h, Yv_i = w_h * vals_sum / n_h)]
  
  Zt <- ests.out[, sum(Zt_i)]
  Yv <- ests.out[, sum(Yv_i)]
  R <- Yv / Zt
  
  means <- list(Zt = Zt, Yv = Yv, R = R)
  return(means)
}


### This function will find, for any species, the estimate for quantile q1. Nested within bs.fcn and q_mort.grow.fcn.
quant.est.spp.fcn <- function(spp.sel, d.all_use, d.vals_use, samp, category.n, q1) {
  
  spp.id1 <- paste0("X", spp.sel)
  n_q <- get(spp.id1, category.n)[q1]   # Are there fewer than N.PLOT.LIM? If so, we can skip the calcs.
  
  # Means plus other metrics that are carried over into the calculation of the variance		
  if (n_q < N.PLOT.LIM) { # If too few points, just enter NA
    q_bs.mean <- NA 
  } else {
    q_bs.mean <- mean.q.fcn(dat_tot = d.all_use[samp$row.id, ], dat_vals = d.vals_use[samp$row.id, ], spp.sel1 = spp.sel)$R # Mean given the selected rows
  }
  return(q_bs.mean)
}


# Bootstrap function for obtaining an estimate of the mean for a single iteration.  Nested within q_mort.grow.fcn.
bs.fcn <- function(iter, d.all_use, d.vals_use, category.n, q1, spp.use) { # Iteration number, all of the relevant trees data, the tree values data (number dead, basal area growth)
  
  # First, sample the tables in use. This procedure obtains samples with replacement from each stratum and then combines the selected
  #  row numbers. 
  spp.use2 <- as.numeric(gsub("X", "", spp.use))

  strata.num.dt <- as.data.table(strata.num) # Using a data.table and then data.table call because much faster than dplyr.
  
  samp <- strata.num.dt[, .SD[sample(.N, .N, replace = TRUE)], by = stratum]
  

  bs.output <- unlist(spp.use2 %>% map(\(x) quant.est.spp.fcn(x, d.all_use, d.vals_use, samp, category.n, q1)))  # For a given iteration, finding the mean value for a given quantile (domain)
  return(bs.output)
}

## Base function for finding the estimated means for mortality/growth for a given quantile.
q_mort.grow.fcn <- function(q1, vals.dat, all.dat, array.name, category.n, selected.spp) {
  # Multiplying all of the species dead values (or species total values) by the respective species masks for "LiLc".
  #   With this step complete, the values for that quantile can be calculated for all species.
  
  vals.use <- vals.dat %>% select(starts_with("nd.")) * array.name[, q1, ] 
  dat.vals.use <- bind_cols(PlotDat, vals.use) 
  
  all.use <- all.dat %>% select(starts_with("nd.")) * array.name[, q1, ] 
  dat.all.use <- bind_cols(PlotDat, all.use) 
  
  # Running the parallel portion of the function
  furrr.out <- 1:BS.N %>% future_map(\(x) bs.fcn(iter = x, d.all_use = dat.all.use, d.vals_use = dat.vals.use, category.n = category.n, q1 = q1, spp.use = selected.spp), .options = furrr_options(seed = TRUE))
  
  # Processing the results and returning quantiles/mean
  furrr.table <- matrix(unlist(furrr.out), ncol = length(selected.spp), byrow = TRUE)
  quants <- data.frame(t(apply(furrr.table, 2, function(x) quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE))))
  names(quants) <- c("Median", "LCI.95", "UCI.95")
  quants$Means <- apply(furrr.table, 2, mean, na.rm = TRUE)
  quants$n.plts <- as.vector(category.n[q1, selected.spp ])
  quants$Species <- selected.spp
  quants$Quantile <- q1
  
  return(quants)
}












#### Graphical output functions -----------------------------------------------






### - Quantile distribution plot, used in multiple analysis pathways - ###
#       Obtain figure of quantile distribution relative to climate info
# plot.spp = species of interest (e.g., "X11"), quant.matrix = quantiles (1 thru 9)
# for each species, var1, var.delt = carried through from input for pairs.plts.fcn,
# quant.lims, quant.n = quantile boundaries and number of plots, established in parse.tree.clim.fcn,
# size.trees = one of "small diameter ", "large diameter ", and (for pathways 1 and 3) "".  
quant.dist.plt.fcn <- function(plot.spp, quant.matrix, var1, var.delt, quant.lims, n.plots.used, size.trees) {
  
  q_plot_spp <- quant.matrix %>% dplyr::select(all_of(plot.spp), all_of(var1), all_of(var.delt)) %>%
    filter(get(plot.spp) > 0) 
  names(q_plot_spp)[1] <- "Quantile"
  
  
  quant.lims.plt <- get(plot.spp, quant.lims)
  quant.lims.delt.plt <- get(plot.spp, quant.lims.delta)
  
  sppnum <- as.numeric(gsub("X", "", plot.spp))
  
  # Common and Genus/species name for plot title
  com.name <- spp.names$COMMON_NAME[spp.names$SPCD == sppnum]
  g.s.name <- paste(spp.names$GENUS[spp.names$SPCD == sppnum], spp.names$SPECIES[spp.names$SPCD == sppnum])
  
  # Need to isolate the correct colors if some quadrants are without plots:
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_quant)  # Get colors for plotting
  scatter.virid.use <- virid.use[n.plots.used$loc[n.plots.used$n > 0]]
  
  fig.lab <- if(ANALYSIS.PATHWAY == 2) {
    paste0("Plot conditions for ", size.trees, com.name, ",\n", g.s.name, ", ",  plot.spp)
  } else {
    paste0("Plot conditions for ", size.trees, com.name, ", ", g.s.name, ", ",  plot.spp)}
  
  plot.vals.plt <- ggplot(q_plot_spp, aes(get(var1), get(var.delt), color = factor(Quantile))) + 
    geom_point() + 
    #stat_ellipse(type = "norm", level = 0.95, col = "orange", lwd = 2) +
    geom_hline(yintercept = quant.lims.delt.plt[1], col = "blue") +
    geom_hline(yintercept = quant.lims.delt.plt[2], col = "blue") +
    geom_vline(xintercept = quant.lims.plt[1], col = "green") +
    geom_vline(xintercept = quant.lims.plt[2], col = "green") +
    geom_hline(yintercept = 0, col = "black") +
    theme_bw() + 
    scale_color_manual(values = scatter.virid.use, name = "Quantiles", labels = QUANT.LEVELS[which(n.plots.used$n > 0)]) + 
    labs(title = fig.lab,
         x = var.label, y = var.delt.label) +
    theme(text = element_text(size = 7),
          legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
  
  return(plot.vals.plt)
}





### -- Mapping elements -- ###

quant.map.fcn <- function(quant.matrix, spp.num, n.plots.used, virid.use) {
  
  map.dat.1 <- quant.matrix %>% 
    mutate(targ.spp = get(spp.num)) %>%
    filter(targ.spp > 0) %>%
    dplyr::select(LAT, LON, targ.spp)
  
  # Points defining the boundaries of the map
  maxlat <- max(map.dat.1$LAT); minlat <- min(map.dat.1$LAT)
  maxlong <- max(map.dat.1$LON); minlong <- min(map.dat.1$LON)
  
  # Refining the color palette based on the number of quadrants occupied.
  #  The scale is already truncated from 0.1 to 0.9, so if 3 is the minimum, 0.3 works well.
  #virid.min <- min(n.plots.used$loc[n.plots.used$n > 0]) / 10
  #virid.max <- max(n.plots.used$loc[n.plots.used$n > 0]) / 10
  map.virid.use <- virid.use[n.plots.used$loc[n.plots.used$n > 0]]
  
  ggplot(data = map.dat.1, aes(x = LON, y = LAT)) +
    coord_fixed(xlim = c(minlong - 1, maxlong + 1),  ylim = c(minlat - 1, maxlat + 1), ratio = 1.3) +
    geom_point(data = map.dat.1, aes(LON, LAT, color = factor(targ.spp)), size = 0.5) +
  #  geom_tile(data = map.dat.1, mapping = aes(x = LON, y = LAT, z = targ.spp), binwidth = 0.15, 
  #            stat = "summary_2d", fun = mean, na.rm = TRUE, show.legend = FALSE) + 
    scale_color_manual(values = map.virid.use) +  
    #                  labels = gsub(".", " ", QUANT.LEVELS[quants.used], fixed = TRUE)) +
    #scale_fill_viridis(option = "H", begin = virid.min, end = virid.max) +
    geom_polygon(data = west_df, mapping = aes(x = long, y = lat, group = group), color = "black", fill = "transparent") +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5),
          legend.position = "none")
}







# -- Plots for ANALYSIS.PATHWAY 1 -- #

## Function for plotting the mortality numbers by quantile and the distribution of plots in quantiles
pair.plts.fcn <- function(sppnum.to.plot, var.filename, quant.table, var.label, var.delt.label, 
                          var1, var.delt, quant.matrix, quant.lims, quant.n){
  
  plot.quant.dat <- quant.table %>% 
    filter(Species == sppnum.to.plot) %>%
    left_join(quant.level.table, by = c("Quantile" = "q.num"))
  
  
  # Details for plotting
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_quant)  # Get colors for plotting
  
  qt.max <- ceiling(max(plot.quant.dat$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
  
  
  # Create a grid of plots to match bivariate plot
  q.g.p.labs <- paste0(plot.quant.dat$var.deltvar, ", n = ", as.numeric(plot.quant.dat$n.plts)) # Labels for the quantile grid plot
  
  ylabs <- if (j == 1) "Mean Growth (inches^2)" else "Mean 10-yr Death Rate"
  
  quant.grid.plt.fcn <- function(quants, quant.index) {
    ggplot(data = plot.quant.dat %>% filter(Quantile %in% quant.index), aes(factor(Quantile), Means, fill = factor(Quantile))) + 
      geom_col() + 
      geom_errorbar(aes(ymax = UCI.95, ymin = LCI.95), width = 0.1) + 
      scale_fill_manual(values = virid.use[quant.index]) + 
      scale_y_continuous(limits = c(0, qt.max)) +
      scale_x_discrete(labels = q.g.p.labs[quant.index]) +
      theme_bw() +
      theme(legend.position = "none") + 
      labs(x = NULL, y = ylabs) +
      theme(text = element_text(size = 7))
  }
  # These will be arranged in a plot below
  p1 <- quant.grid.plt.fcn(c("LiHc", "MiHc", "HiHc"), 1:3) 
  p2 <- quant.grid.plt.fcn(c("LiMc", "MiMc", "HiMc"), 4:6) 
  p3 <- quant.grid.plt.fcn(c("LiLc", "MiLc", "HiLc"), 7:9) 
     

  quant.table1 <- quant.table %>% filter(Species == sppnum.to.plot)
  
  # Code to prep for next two figures - helps in reducing the color palette.
  n_plots <- get(sppnum.to.plot, quant.n)
  n_plots2 <- tibble(loc = 1:n_quant, n = n_plots) %>%
    left_join(tibble(Quantiles = plot.quant.dat$Quantiles, loc = 1:n_quant), by = "loc")
  
  # Plotting the scatterplot of points in quadrants:
  plot.vals.plt <- quant.dist.plt.fcn(plot.spp = sppnum.to.plot, quant.matrix = quant.matrix, 
                                      var1 = var1, var.delt = var.delt, quant.lims = quant.lims,
                                      n.plots.used = n_plots2, size.trees = "") 
  
  # Plotting the map of plot locations:
  quant.map <- quant.map.fcn(quant.matrix, spp.num = sppnum.to.plot, n.plots.used = n_plots2, virid.use = virid.use) 
  
  comb.plt <- ((p1/p2/p3 + plot_layout(axis_titles = "collect")) | plot.vals.plt | quant.map) #/ guide_area() + plot_layout(guides = 'collect', heights = c(10, 0.01)) 
  
  #ggsave(paste0(RESULTS.LOC, "Mort_figs_", var.filename, "/",sppnum.to.plot, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")
  ggsave(paste0(RESULTS.LOC, switch(ANALYSIS.PATHWAY, "1" = "Quantile_Only/", "2" = "Size_Class/", "3" = "Site_Class/"), 
                switch(ANALYSIS.TYPE[j], "mort" = "Mort", "grow" = "Growth"), 
                "_figs_", var.filename, "/", sppnum.to.plot, "_plots.png"), 
         comb.plt, device = "png", width = 7, height = 5, units = "in")
}



# -- Plots for ANALYSIS.PATHWAY 2 -- #

## Function for plotting the mortality numbers by quantile and the distribution of plots in quantiles

pair.plts2.fcn <- function(sppnum.to.plot, quant.table, var.label, var.delt.label, var1, var.delt,
                           quant.matrix.1, quant.matrix.2, quant.lims.1, quant.lims.2,
                           quant.n.1, quant.n.2){
  
  plot.quant.dat <- quant.table %>% 
    filter(Species == sppnum.to.plot) %>%
    left_join(quant.level.table, by = c("Quantile" = "q.num"))
  
  
  # Details for plotting
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_quant)  # Get colors for plotting
  
  qt.max <- ceiling(max(plot.quant.dat$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
  
  
  quant.grid.plt.fcn <- function(quants, quant.index) {
    
    num.labs <- plot.quant.dat %>% filter(Quantile %in% quant.index) %>%
      mutate(loc = 0.8 * LCI.95,
             lab = as.numeric(n.plts)) %>%
      arrange(Quantile, smaller_larger)
    
    
    ggplot(data = plot.quant.dat %>% filter(Quantile %in% quant.index), aes(Quantile, Means, fill = factor(smaller_larger), color = factor(Quantile))) + 
      geom_col(position = "dodge", width = 0.75, linewidth = 2) + 
      geom_errorbar(aes(ymax = UCI.95, ymin = LCI.95), width = 0.2, position = position_dodge(width = 0.75), color = 'black') + 
      geom_text(data = num.labs, aes(label = lab, x = Quantile, y = loc), position = position_dodge(width = 0.75), color = "black", size = 2) +
      scale_color_manual(values = virid.use[quant.index]) + 
      scale_fill_manual(values = c("#C1CDCD", "#838B8B")) +
      scale_y_continuous(limits = c(0, qt.max)) +
      scale_x_continuous(labels = quants, breaks = quant.index) +
      theme_bw() +
      theme(legend.position = "none") + 
      labs(x = NULL, y = "Mean") +
      theme(text = element_text(size = 7))
  }
  
  p1 <- quant.grid.plt.fcn(c("LiHc", "MiHc", "HiHc"), 1:3) 
  p2 <- quant.grid.plt.fcn(c("LiMc", "MiMc", "HiMc"), 4:6) 
  p3 <- quant.grid.plt.fcn(c("LiLc", "MiLc", "HiLc"), 7:9) 
  
  p_all <- plot_grid(p1, p2, p3, ncol = 1)
  
  plot.vals.plt.1 <- quant.dist.plt.fcn(plot.spp = sppnum.to.plot, quant.matrix = quant.matrix.1, 
                                        var1 = var1, var.delt = var.delt, quant.lims = quant.lims.1,
                                        quant.n = quant.n.1, size.trees = "smaller diameter ") 
  
  plot.vals.plt.2 <- quant.dist.plt.fcn(plot.spp = sppnum.to.plot, quant.matrix = quant.matrix.2, 
                                        var1 = var1, var.delt = var.delt, quant.lims = quant.lims.2,
                                        quant.n = quant.n.2, size.trees = "larger diameter ") 
  
  comb.plt <- p_all +(plot.vals.plt.1 + plot.vals.plt.2 + plot_layout(ncol = 1)) + plot_layout( widths = c(1, 1.2))
  
  #ggsave(paste0(RESULTS.LOC, "Mort_figs_", var.filename, "/",sppnum.to.plot, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")
  ggsave(paste0(RESULTS.LOC, switch(ANALYSIS.PATHWAY, "1" = "Quantile_Only/", "2" = "Size_Class/", "3" = "Site_Class/"), 
                switch(ANALYSIS.TYPE[j], "mort" = "Mort", "grow" = "Growth"), 
                "_figs_", var.filename, "/", sppnum.to.plot, "_plots.png"), 
         comb.plt, device = "png", width = 7, height = 5, units = "in")
  
}




# -- Plots for ANALYSIS.PATHWAY 3 -- #

## Function for plotting the mortality numbers by quantile and the distribution of plots in quantiles
pair.plts3.fcn <- function(sppnum.to.plot, quant.table, var.label, var.delt.label, 
                           var1, var.delt, quant.matrix, quant.n){
  
  plot.siteclass.dat <- quant.table %>% 
    filter(Species == sppnum.to.plot) 
  
  # Details for plotting
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_quant)  # Get colors for plotting
  
  qt.max <- ceiling(max(plot.siteclass.dat$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
  
  
  # Create a grid of plots to match bivariate plot
  #q.g.p.labs <- paste0("n = ", as.numeric(plot.siteclass.dat$n.plts)) # Text for plot number
  
  num.labs <- plot.siteclass.dat %>% filter(Quantile %in% quant.index) %>%
    mutate(loc = (qt.max / 30) + UCI.95,
           lab = paste0("n = ", as.numeric(n.plts)))
  
  ### --- Plotting Site Class
  
  siteclass.plt <- ggplot(plot.siteclass.dat, aes(factor(Quantile), Means, fill = factor(Quantile))) + 
    geom_col() + 
    geom_errorbar(aes(ymax = UCI.95, ymin = LCI.95), width = 0.1) + 
    geom_text(data = num.labs, aes(label = lab, x = Quantile, y = loc), color = "black", size = 2) +
    scale_fill_manual(values = virid.use[quant.index]) + 
    scale_y_continuous(limits = c(0, qt.max)) +
    #scale_x_discrete(labels = q.g.p.labs[quant.index]) +
    theme_bw() +
    theme(legend.position = "none") + 
    labs(x = "Site Class", y = "Mean") +
    theme(text = element_text(size = 7))
  
  
  ### --- Scatterplot of delta.var and var with site class = color ##
  classes.used <- plot.siteclass.dat$Quantile[is.na(plot.siteclass.dat$Means) == FALSE]
  
  q_plot_spp <- quant.matrix %>% dplyr::select(all_of(sppnum.to.plot), all_of(var1), all_of(var.delt)) %>%
    filter(get(all_of(sppnum.to.plot)) %in% classes.used) 
  names(q_plot_spp)[1] <- "SiteClass"
  
  sppnum <- as.numeric(gsub("X", "", sppnum.to.plot))
  
  # Common and Genus/species name for plot title
  com.name <- spp.names$COMMON_NAME[spp.names$SPCD == sppnum]
  g.s.name <- paste(spp.names$GENUS[spp.names$SPCD == sppnum], spp.names$SPECIES[spp.names$SPCD == sppnum])
  fig.lab <- paste0("Plot conditions for ", com.name, ", ", g.s.name, ", ",  sppnum.to.plot)
  
  
  plot.vals.plt <- ggplot(q_plot_spp, aes(get(var1), get(var.delt), color = factor(SiteClass))) + 
    geom_point() + 
    stat_ellipse(linewidth = 1) +  # Radius = level = 0.95 (default), type = "t" (default) = t-dist, 
    geom_hline(yintercept = 0) +
    theme_bw() + 
    scale_color_manual(values = virid.use[classes.used], name = "Site Classes", 
                       labels = gsub(".", " ", QUANT.LEVELS[classes.used], fixed = TRUE)) + 
    labs(title = fig.lab,
         x = var.label, y = var.delt.label) +
    theme(text = element_text(size = 7))
  
  
  # Now combining the two plots:
  comb.plt <- siteclass.plt + plot.vals.plt + plot_layout(widths = c(1, 1.3))
  
  #ggsave(paste0(RESULTS.LOC, "Mort_figs_", var.filename, "/",sppnum.to.plot, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")
  ggsave(paste0(RESULTS.LOC, switch(ANALYSIS.PATHWAY, "1" = "Quantile_Only/", "2" = "Size_Class/", "3" = "Site_Class/"), 
                switch(ANALYSIS.TYPE[j], "mort" = "Mort", "grow" = "Growth"), 
                "_figs_", var.filename, "/", sppnum.to.plot, "_plots.png"), 
         comb.plt, device = "png", width = 7, height = 5, units = "in")
}

