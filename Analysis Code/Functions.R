
##### ---- Functions ---- #####

### DATA PREP FUNCTIONS - MORTALITY ####


## Domain masks for each species. Used in clim.tree.resp.fcn (called in Data_Prep.R), operates on a single domain level and species 
domain.divide.fcn <- function(domain.level, spcd2, tree.plots2, PlotDat) {
  domain.name <- paste0(domain.level, ".",  spcd2)
  domain.died.out <- left_join(PlotDat, tree.plots2 %>% 
                                 dplyr::select(puid, var.deltvar) %>% # var.deltvar defined in clim.tree.resp.fcn, "LiHc" etc.
                                 filter(var.deltvar == domain.level), by = "puid") %>%
    mutate(var.deltvar = ifelse(is.na(var.deltvar), 0, 1)) %>%
    rename(!!domain.name := var.deltvar) 
  return(domain.died.out)
}



#### Data prep function for mortality --------------------------------------------------

## Function to determine the number of trees that died in a plot and their climate domains, used in Data_Prep.R.
clim.mort.resp.fcn <- function(spcd, clim.var, treedat.sel, clim.dat) {
  var1 <- "pre_mean" #paste0("pre.", clim.var)
  var.delt <- "difference"  #paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  # For alive/dead tree assessments
  Tree2 <- treedat.sel %>% 
    filter(SPCD == spcd) %>%
    dplyr::select(puid, REMPER, SUBP, TREE, STATUSCD, ntree.rep)
  
  tree.plots <- clim.dat %>% filter(puid %in% Tree2$puid) %>% 
    dplyr::select(puid, all_of(var1), all_of(var.delt))
  
  # Data quantiles
  # Setting the delta-var amounts to be fixed values.  See Global for setting VAR.DELTA.BOUNDARIES.
  quant.lims.delta.pos <- VAR.DELTA.BOUNDARIES$max.min[VAR.DELTA.BOUNDARIES$clim.var == clim.var] #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- -quant.lims.delta.pos #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
  quant.lims.delta <- if(n_domain ==9) c(quant.lims.delta.neg, quant.lims.delta.pos) else quant.lims.delta.pos
  
  # Actual quantiles are used if the Global settings for USE.QUANT.PROBS is TRUE.  Otherwise the Global-defined INIT.CLIM.BREAKS are used across species.
  quant.lims <- if(USE.QUANT.PROBS) {
    quantile(get(var1, tree.plots), probs = QUANT.PROBS)
  } else{
    INIT.CLIM.BREAKS
  }
  
  # Applying quantiles to a particular species
  tree.plots2 <- tree.plots %>% dplyr::select(puid, all_of(c(var1, var.delt))) %>% 
      {if (n_domain == 9) {
        mutate(., var.quants = ifelse(!!sym(var1) > quant.lims[2], 1, ifelse(!!sym(var1) < quant.lims[1], -1, 0)),
               deltvar.quants = ifelse(!!sym(var.delt) > quant.lims.delta[2], 1, ifelse(!!sym(var.delt) < quant.lims.delta[1], -1, 0)),
               vq = case_match(var.quants, -1 ~ "Li", 0 ~ "Mi", 1 ~ "Hi"),
               dvq = case_match(deltvar.quants, -1 ~ "Lc", 0 ~ "Mc", 1 ~ "Hc"),
               var.deltvar = factor(paste0(vq, dvq), levels = DOMAIN.LEVELS)) 
      } else if (n_domain == 6) {
        mutate(., var.quants = ifelse(!!sym(var1) > quant.lims[2], 1, ifelse(!!sym(var1) < quant.lims[1], -1, 0)),
               deltvar.quants = ifelse(!!sym(var.delt) > quant.lims.delta, 1, 0),
               vq = case_match(var.quants, -1 ~ "L", 0 ~ "M", 1 ~ "H"),
               dvq = case_match(deltvar.quants, 0 ~ "B", 1 ~ "A"),
               var.deltvar = factor(paste0(dvq, vq), levels = DOMAIN.LEVELS))
      }}
        
    # Not currently used, but available    
    # For centroid calcluation
    dat.cent <- tree.plots2 %>% group_by(var.deltvar) %>% # Quartile centers, good for finding centroids in each quadrant
      reframe(end.x = mean(get(var1)),        # Referred to as "end" because of distance measurements from the main centroid to each of the outer (end) centroids.
              end.y = mean(get(var.delt))) %>% 
      right_join(domain.level.table %>% select(-q.num), by = "var.deltvar") %>%
      arrange(var.deltvar)
    
   
  
  # Joining tree data with domain assignment, providing plot-level summary of the number of trees, number that died, and percent that died.
  
  # Mortality data
  Tree3 <- left_join(Tree2, tree.plots2, by = "puid") %>%
    group_by(puid) %>% #, var.deltvar, get(var1), get(var.delt)) %>%  # This will assign the names `get(var1)` and `get(var.delt)`, which are cleaned up below
    reframe(n.trees = sum(ntree.rep),
            #n.died.yr = sum(ntree.rep[STATUSCD == 2])/mean(REMPER),  # Annual mortality rate
            n.died.dcd = 10 * sum(ntree.rep[STATUSCD == 2])/mean(REMPER), # Decadal mortality rate
            #pct.died.yr = n.died.yr / n.trees)
            pct.died.dcd = n.died.dcd / n.trees)
  
  # Separate out a data set for those that died
  died.out <- left_join(PlotDat, Tree3 %>% dplyr::select(puid, n.died.dcd), by = "puid") %>%
    rename(!!return.col := n.died.dcd)
  died.out[is.na(died.out)] <- 0
  
  # Data set for all trees, alive and died
  all.trees <- left_join(PlotDat, Tree3 %>% dplyr::select(puid, n.trees), by = "puid") %>%
    rename(!!return.col := n.trees)
  all.trees[is.na(all.trees)] <- 0
  
  
  domain.out1 <- purrr::map(DOMAIN.LEVELS, domain.divide.fcn, spcd2 = spcd, tree.plots2 = tree.plots2, PlotDat = PlotDat)
  domain.out2 <- reduce(domain.out1, left_join, by = c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts")) %>%
    dplyr::select(-c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, died.out = died.out, all.trees = all.trees, 
                   domain.out = domain.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}


#### Data prep function for growth  --------------------------------------------------

## Function to determine the number of trees that died in a plot and their climate domains, used in Data_Prep.R.
clim.growth.resp.fcn <- function(spcd, clim.var, treedat.sel, clim.dat) {
  var1 <- "pre_mean" #paste0("pre.", clim.var)
  var.delt <- "difference"  #paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  Tree.spp <- treedat.sel %>% filter(SPCD == spcd, STATUSCD == 1)
  
  # For basal area growth/tree
  TreeG <- Tree.spp %>%
    # Amount of tree growth per year since the previous visit, multiplied by the number of trees per acre the tree represents. 
    #mutate(BA_Growth = ntree.rep * (pi * (DIA / 2)^2 - pi * (PREVDIA / 2)^2)/REMPER) %>%
    # Amount of tree growth per decade since the previous visit, multiplied by the number of trees per acre the tree represents. 
    mutate(BA_Growth = ntree.rep * 10 * (pi * (DIA / 2)^2 - pi * (PREVDIA / 2)^2)/REMPER) %>%
    group_by(puid) %>%
    #Goal: sum(change in BA/REMPER)/Total # trees
    reframe(n_g.trees = sum(ntree.rep),
            total.growth = sum(BA_Growth))
  
  tree.plots <- Tree.spp %>% 
    dplyr::select(puid, all_of(var1), all_of(var.delt)) %>% distinct()
  
  # Data quantiles
  quant.lims.delta.pos <- VAR.DELTA.BOUNDARIES$max.min[VAR.DELTA.BOUNDARIES$clim.var == clim.var] #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- -quant.lims.delta.pos #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
  quant.lims.delta <- if(n_domain ==9) c(quant.lims.delta.neg, quant.lims.delta.pos) else quant.lims.delta.pos
  # Actual quantiles are used if the Global settings for USE.QUANT.PROBS is TRUE.  Otherwise the Global-defined INIT.CLIM.BREAKS are used across species.
  quant.lims <- if(USE.QUANT.PROBS) {
    quantile(get(var1, tree.plots), probs = QUANT.PROBS)
  } else{
    INIT.CLIM.BREAKS
  }
  
  # Applying quantiles to a particular species
  if(ANALYSIS.PATHWAY != 3) {
    tree.plots2 <- tree.plots %>% dplyr::select(puid, all_of(c(var1, var.delt))) %>% 
      {if (n_domain == 9) {
        mutate(., var.quants = ifelse(!!sym(var1) > quant.lims[2], 1, ifelse(!!sym(var1) < quant.lims[1], -1, 0)),
               deltvar.quants = ifelse(!!sym(var.delt) > quant.lims.delta[2], 1, ifelse(!!sym(var.delt) < quant.lims.delta[1], -1, 0)),
               vq = case_match(var.quants, -1 ~ "Li", 0 ~ "Mi", 1 ~ "Hi"),
               dvq = case_match(deltvar.quants, -1 ~ "Lc", 0 ~ "Mc", 1 ~ "Hc"),
               var.deltvar = factor(paste0(vq, dvq), levels = DOMAIN.LEVELS)) 
      } else if (n_domain == 6) {
        mutate(., var.quants = ifelse(!!sym(var1) > quant.lims[2], 1, ifelse(!!sym(var1) < quant.lims[1], -1, 0)),
               deltvar.quants = ifelse(!!sym(var.delt) > quant.lims.delta, 1, 0),
               vq = case_match(var.quants, -1 ~ "L", 0 ~ "M", 1 ~ "H"),
               dvq = case_match(deltvar.quants, 0 ~ "B", 1 ~ "A"),
               var.deltvar = factor(paste0(dvq, vq), levels = DOMAIN.LEVELS))
      } else {
        values$errorMessage <- "Error: Expecting 6 or 9 estimation domains. Check Global.R:DOMAIN.LEVELS and Functions.R:clim.mort.resp.fcn."
        return(.)
      }}
    
    # For centroid calcluation
    dat.cent <- tree.plots2 %>% group_by(var.deltvar) %>% # Quartile centers, good for finding centroids in each quadrant
      reframe(end.x = mean(get(var1)),        # Referred to as "end" because of distance measurements from the main centroid to each of the outer (end) centroids.
              end.y = mean(get(var.delt))) %>% 
      right_join(domain.level.table %>% select(-q.num), by = "var.deltvar") %>%
      arrange(var.deltvar)
    
  } else {
    # For Analysis Pathway 3.  A different set-up - not using the 9 domains, but in effect 7 (for Site Classes). The plots will be much simpler.
    tree.plots2 <- tree.plots %>% left_join(PlotDat %>% select(puid, SITECLCD_plot), by = "puid") %>% # Adding on the Site Class variable
      filter(is.na(SITECLCD_plot) == FALSE) %>%   # Removing any NA's in case they snuck in there.
      mutate(var.deltvar = factor(paste0("Class.", SITECLCD_plot)))
    
    dat.cent <- NULL
    
    DOMAIN.LEVELS <- levels(tree.plots2$var.deltvar) # Renaming the DOMAIN.LEVELS away from "LiHc" and such
  }
  
  
  
  # Joining tree data with domain assignment, providing plot-level summary of the number of trees, number that died, and percent that died.
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
  
  
  domain.out1 <- purrr::map(DOMAIN.LEVELS, domain.divide.fcn, spcd2 = spcd, tree.plots2 = tree.plots2, PlotDat = PlotDat)
  domain.out2 <- reduce(domain.out1, left_join, by = c("STATECD", "puid", "n_h.plts", "ESTN_UNIT", "STRATUMCD", "w", "stratum")) %>%
    dplyr::select(-c("STATECD", "puid", "n_h.plts", "ESTN_UNIT", "STRATUMCD", "w", "stratum"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, growth.val = growth.val, growth.n.trees = growth.n.trees,
                   domain.out = domain.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}


#### Preparing mortality/growth output for analysis -------------------------------


## This is a large function that takes species-specific data and creates arrays to assist with analyses.
parse.tree.clim.fcn <- function(tree.dat, clim.var, analysis.type, resp.dat, tot.dat, selected.spp, clim.dat) {  
  # tree.dat = prepped list of data, clim.var = which climate variable name, 
  #  analysis.type = "mort" or "grow", resp.dat = response data set (i.e., 
  # "died.out" or "growth.val"), tot.dat = total trees data set (i.e., "all.trees" or "growth.n.trees")
  # Makes a list of all of the response tables and all-tree tables for the next step to operate on.
  
  extract.resp <- purrr::map(tree.dat , ~.[[resp.dat]]) 
  vals_dat <- reduce(extract.resp, left_join, by = c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts")) # Combining the list into a single tibble.
  #vals_dat <- reduce(extract.resp, left_join, by = c("STATECD", "PLOT_FIADB", "puid", "W_h", "STRATUM","SITECLCD_plot")) # Combining the list into a single tibble.
  
  extract.all <- purrr::map(tree.dat , ~.[[tot.dat]]) 
  all_dat <- reduce(extract.all, left_join, by = c("STATECD", "puid", "ESTN_UNIT", "STRATUMCD", "w", "stratum", "n_h.plts")) 
  
  
  # Now extracting the domain category information for each species 
  extract.domain <- purrr::map(tree.dat , ~.[["domain.out"]]) 
  # Creating a 3D array out of the domain info (plots x domains x spp)
  domain.array <- array(unlist(extract.domain), dim = c(nrow(extract.domain[[1]]), ncol(extract.domain[[1]]), length(extract.domain)), 
                        dimnames = list(c(1:nrow(extract.domain[[1]])),
                                        DOMAIN.LEVELS, 
                                        selected.spp))
  
  # Matrix of domain values for each species' plots joined with climate data
  # The values for plots with trees of a species are given the ID number for their respective quantile (e.g., i = 9 has plot values of 0 or 9).
  # Note that this is just for the construction of the quantile matrix (values 0 through 9), not quant.array.
  q2 <- domain.array
  for (i in 1:n_domain) {
    q2[, i, ] <- q2[, i, ] * i
  }
  domain.matrix <- apply(q2, c(1, 3), sum) %>%  # collapsing the domain values into a single matrix
    bind_cols(vals_dat %>% dplyr::select(puid)) %>% # It does not matter whether vals_dat or all_dat, just want puid
    left_join(clim.dat, by =  "puid")
  
  # Number of plots in each domain
  domain.n <- purrr::map(1:length(extract.domain), function(x) apply(domain.array[, , x], 2, sum))
  domain.n <- bind_cols(domain.n)
  names(domain.n) <- selected.spp
  
  extract.quant.lims <- purrr::map(tree.dat, ~.[["quant.lims"]]) 
  names(extract.quant.lims) <- selected.spp
  
  extract.delt.quant.lims <- purrr::map(tree.dat, ~.[["quant.lims.delta"]]) 
  names(extract.delt.quant.lims) <- selected.spp
  
  
  # Now extracting the centroid information for each domain/species)
  if(ANALYSIS.PATHWAY != 3) {
    extract.centroid <- purrr::map(tree.dat, ~.[["cent.loc"]]) 
    centroid.array <- array(unlist(extract.centroid), dim = c(nrow(extract.centroid[[1]]), ncol(extract.centroid[[1]]), length(extract.centroid)), 
                            dimnames = list(DOMAIN.LEVELS, 
                                            c("var.deltvar", "end.x", "end.y"),
                                            selected.spp))
  } else {
    extract.centroid <- centroid.array <-  NULL
  }
  
  ### State level ###
  state.matrix <- domain.matrix[, 1:length(selected.spp)]
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
  
  extract.state <- purrr::map(state.list, state.array.fcn, state.matrix)
  # Creating a 3D array out of the domain info (plots x domains x spp)
  state.array <- array(unlist(extract.state), dim = c(nrow(extract.state[[1]]), ncol(extract.state[[1]]), length(state.list)),
                       dimnames = list(c(1:nrow(extract.state[[1]])),
                                       selected.spp,
                                       state.list))
  state.array <- aperm(state.array, c(1, 3, 2))  # Changing the order of the dimensions to match those for domains 
  
  # Obtaining the number of plots per state.
  state.n.extract <- purrr::map(1:3, function(x) apply(extract.state[[x]], 2, sum)) #sum plots by state.
  state.n <- reduce(state.n.extract, bind_rows) %>%  # Combine state data into a tibble.
    mutate(state = state.list)
  
  
  
  ### All data level ###
  all.matrix <- domain.matrix[, 1:length(selected.spp)] 
  all.matrix <- apply(all.matrix, 2, function(x) ifelse(x > 0, 1, 0)) %>% 
    data.frame()
  all.array <- abind(all.matrix, all.matrix, along = 3)
  all.array <- aperm(all.array, c(1, 3, 2))  # Changing the order of the dimensions to match those for domains 
  
  
  all.n <- apply(all.array[, 1, ], 2, sum)
  all.n <- tibble(n = all.n, labs = selected.spp) %>% 
    pivot_wider(names_from = "labs", values_from = "n") %>%
    mutate(all = 1)
  
  ## We now have the data at the level needed to conduct the mortality analysis.  We need to select a species and a domain and find the 
  #  estimated mortality rate for that quantile.  The domain array (domain.array) offers a mask to remove all values other than the domain 
  #  plots of interest for the species. 
  
  
  mort.grow.dat <- list(
    spp.list = tree.dat[[1]]$spp.list,
    spp.id = tree.dat[[1]]$spp.id,
    vals_dat = vals_dat, 
    all_dat = all_dat, 
    domain.array = domain.array, 
    state.array = state.array,
    all.array = all.array,
    domain.matrix = domain.matrix,
    centroid.array = centroid.array,
    domain.n = domain.n,
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
  
  y <- parsed.dat$domain.matrix %>% dplyr::select(all_of(tspp), puid) %>%
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
  fire.frac.dat <- purrr::map(SEL.SPP, fire.frac.dat.fcn, treedat.sel  = treedat, parsed.dat = parseddat)
  fire.frac.dat2 <- rbindlist(fire.frac.dat) %>% data.frame()   # data.table function
  fwrite(fire.frac.dat2, file = paste0(tableloc, tablename))    # data.table csv write function
}




# Determine where to save results (data and plots)
save.loc.fcn <- function(i){
  paste0(RESULTS.LOC, "Quantile_Only/", 
         switch(i, "1" = "Growth", "2" = "Mort"), 
         "_figs_", CLIM.VAR.USE, "/")
}








### ------ Mortality and growth estimation -------  ####

#  This function is nested within calculate_domain_estimates.fcn, which is in turn nested within generate_bootstrap_array.fcn.  
#  This function finds the estimates for mu_hat (the ratio estimator), Y_hat (numerator, thing we're interested in), and 
#    X_Hat (denominator, density of things like total trees of a species). See Equation 6 and related.
mean.q.fcn <- function(dat_tot, dat_vals, spp.sel1) {
  col.name <- paste0("nd.", spp.sel1)
  
  dt_tot <- as.data.table(dat_tot)
  dt_vals <- as.data.table(dat_vals)
  
  tot_summary <- dt_tot[, .(w_h = w[1], n_h = n_h.plts[1], tot_sum = sum(get(col.name))), by = stratum]
  vals_summary <- dt_vals[, .(vals_sum = sum(get(col.name))), by = stratum]
  
  ests.out <- tot_summary[vals_summary, on = "stratum"]
  ests.out[, `:=`(Zt_i = w_h * tot_sum / n_h, Yv_i = w_h * vals_sum / n_h)]
  
  # Old slow code.  For understanding the data.table usage.
  #  Zt_i <- Yv_i <- rep(0, length(strata) )  # Weighted values for each strata
  #strat.num <- length(unique(strata))  
  #for (h in 1:strat.num) {                  # h indexes strata
  #    w_h <- unique(dat_tot$w[dat_tot$stratum == strata[h]])
  #    n_h <- unique(dat_tot$n_h.plts[dat_tot$stratum == strata[h]])    # Number of plots in stratum h 
  #    Yv_i[h] <- w_h * sum(get(col.name, dat_vals)[dat_vals$stratum == strata[h]]) / n_h 
  #    Zt_i[h] <- w_h * sum(get(col.name, dat_tot)[dat_tot$stratum == strata[h]]) / n_h     
  #  }
  
  Zt <- ests.out[, sum(Zt_i)]
  Yv <- ests.out[, sum(Yv_i)]
  R <- Yv / Zt
  
  means <- list(Zt = Zt, Yv = Yv, R = R)
  return(means)
}


# Calculate estimates for all species in a given domain using map
calculate_domain_estimates.fcn <- function(samp, domain_data, domain_idx, selected.spp, domain.n) {
  
  d.all_use <- domain_data[[domain_idx]]$all
  d.vals_use <- domain_data[[domain_idx]]$vals
  
  # Convert species names to numeric for domain.n lookup
  spp.use.numeric <- as.numeric(gsub("X", "", selected.spp))
  
  # Use map to calculate estimates for all species
  estimates <- spp.use.numeric %>% 
    map_dbl(\(spp.sel) {
      # Check if we have enough plots for this species in this domain. Otherwise returns NA.
      spp.id1 <- paste0("X", spp.sel)
      n_q <- get(spp.id1, domain.n)[domain_idx]
      
      if (n_q < N.PLOT.LIM) {
        return(NA_real_)
      } else {
       # return(mean.q.fcn(
        mean.q.fcn(
          dat_tot = d.all_use[samp$row.id, ],
          dat_vals = d.vals_use[samp$row.id, ],
          spp.sel1 = spp.sel
        )$R
      }
    })
  

 # df <- rbind(estimates, estimates.yv, estimates.zt) %>% data.frame()
#  colnames(df) <- SEL.SPP
 # write.csv(df, "Oregon.csv")

  return(estimates)
}


# Generate bootstrap sample.  Each stratum is separately sampled by the stratum size.   
generate_bootstrap_sample.fcn <- function() {
  strata.num.dt <- as.data.table(strata.num)
  
  #samp <- strata.num.dt[, .SD[sample(.N, .N, replace = TRUE)], by = stratum] # If we want to sample within strata
  samp <- strata.num.dt[, .SD[sample(.N, .N, replace = TRUE)],]  # Sampling all plots regardless of strata
  
  return(samp)
}


# Helper function to prepare data for all domains using map
prepare_domain_data.fcn <- function(vals.dat, all.dat, domain.array, n_domains) {
  1:n_domains %>% 
    purrr::map(\(d) {
      vals.use <- vals.dat %>% select(starts_with("nd.")) * domain.array[, d, ]
      dat.vals.use <- bind_cols(PlotDat, vals.use)
      
      all.use <- all.dat %>% select(starts_with("nd.")) * domain.array[, d, ]
      dat.all.use <- bind_cols(PlotDat, all.use)
      
      list(vals = dat.vals.use, all = dat.all.use)
    })
}


# Main function to generate the bootstrap array using parallel processing
generate_bootstrap_array.fcn <- function(vals.dat, all.dat, domain.array, domain.n, selected.spp, n_iter) {
  
  n_species <- length(selected.spp)
  n_domains <- if(ncol(domain.array) == 2) 1 else ncol(domain.array) # If == 2, then we are after
  # single estimates for species.  One of the two domains is there to keep the array three dimensional.
  
  # Pre-process data for all domains
  domain_data <- prepare_domain_data.fcn(vals.dat, all.dat, domain.array, n_domains)
  
  # Run all bootstrap iterations in parallel
  bootstrap_results <- 1:n_iter %>% 
    future_map(\(iter) {
      # Generate single bootstrap sample for this iteration
      samp <- generate_bootstrap_sample.fcn()
      
      # Calculate estimates for all domains for this iteration
      1:n_domains %>% 
        purrr::map(\(d) calculate_domain_estimates.fcn(
          samp = samp,
          domain_data = domain_data,
          domain_idx = d,
          selected.spp = selected.spp,
          domain.n = domain.n
        )) %>%
        do.call(rbind, .) # Combine domains into matrix
    }, .options = furrr_options(seed = TRUE))
  
  # Convert list of matrices to 3D array
  bootstrap_array <- array(
    data = unlist(bootstrap_results),
    dim = c(n_domains, n_species, n_iter),
    dimnames = list(
      domain = 1:n_domains,
      species = selected.spp,
      iteration = 1:n_iter
    )
  )
  
  bootstrap_array2 <- aperm(bootstrap_array, perm = c(3, 2, 1))  
  
  return(bootstrap_array2)
}





## Functions for preparing plotting data

# This function is used to find the quantiles of a bootstrap matrix. Used by both
#  domain.sum.fcn and domain.diff.fcn.
matrix.summ.fcn <- function(matrix.use, domain.id, domain_n){
  quants <- data.frame(t(apply(matrix.use, 2, function(x) quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE))))
  names(quants) <- c("Median", "LCI.95", "UCI.95")
  quants$Means <- apply(matrix.use, 2, mean, na.rm = TRUE)
  quants$n.plts <- unlist(as.vector(domain_n[domain.id, SEL.SPP ]))
  quants$Species <- SEL.SPP
  quants$Domain <- domain.id
  return(quants)
}

# Used in summarizing individual domains
domain.sum.fcn <- function(results.array, domain.num, domain_n) {
  matrix_use <- results.array[, , domain.num]
  matrix.summ.fcn(matrix_use, domain.num, domain_n)
}

# Used to summarize differences in two domains
domain.diff.fcn <- function(start.dom, sub.dom, results.array, domain_n) {
  start.matrix <- results.array[, , domain.level.table$q.num[domain.level.table$var.deltvar == start.dom]]
  sub.matrix <- results.array[, , domain.level.table$q.num[domain.level.table$var.deltvar == sub.dom]]
  
  result.matrix <- start.matrix - sub.matrix
  matrix.summ.fcn(result.matrix, paste(start.dom, "-", sub.dom), domain_n)
}


# This function binds differenced matrices for the creation of plots. Only used by 6-domain analysis.
diff.fig.prep.fcn <- function(diff.table, diff.levels){ 
  diff2.table <- diff.table %>%
    left_join(spp.names.fig, by = "Species") %>%
    select(-n.plts) %>%
    mutate(Domain = factor(Domain, levels = diff.levels),
           species_label = paste(GENUS, SPECIES), 
           significant = LCI.95 > 0 | UCI.95 < 0,
           y.offset = case_when(
             as.numeric(Domain) == 1 ~ 0.15,
             as.numeric(Domain) == 2 ~ 0, 
             as.numeric(Domain) == 3 ~ -0.15,
             .default = 0)
    ) %>%
    filter(is.na(significant) == FALSE) %>%
    # Order species by their position in the original species file
    arrange(match(Species, spp.names.fig$Species))
  
  diff2.table$species_label <- factor(diff2.table$species_label, 
                                      levels = rev(unique(diff2.table$species_label)))
  return(diff2.table)
}







#### Graphical output functions -----------------------------------------------

## Used in preparing the differenced plot, changing in^2 for Growth to cm^2
cm2.fcn <- function(k, results.table){
  if(k == 1){
    cols_sqrd <- which(names(results.table) %in% c("Median", "LCI.95", "UCI.95", "Means"))
    tx <- results.table %>% mutate(across(all_of(cols_sqrd), ~ .x * 6.4516))
    return(tx)
  } else {
    return(results.table)
  }
}



## This plot provides the 
diff.panel.fcn <- function(diff.dat, remove.y, fig.title, lab.right) {
  
  # Kludgy way to get a conifer/deciduous break.  X312 = bigleaf maple.
  decid.break <- length(SEL.SPP) - which(SEL.SPP == "X312") + 1.5
  
  ggplot(diff.dat, aes(y = as.numeric(species_label) + y.offset)) + 
    geom_vline(xintercept = 0, linewidth = 0.5) + 
    geom_segment(aes(x = LCI.95, xend = UCI.95, 
                     y = as.numeric(species_label) + y.offset, 
                     yend = as.numeric(species_label) + y.offset,
                     color = Domain), 
                 linewidth = 0.8) + 
    geom_point(aes(x = Means, y = as.numeric(species_label) + y.offset,
                   fill = ifelse(significant, "black", "white"),
                   color = Domain),
               size = 2, shape = 21, stroke = 0.8) +
    geom_hline(yintercept = decid.break, linetype = 2) +
    #scale_y_continuous(breaks = 1:length(levels(DvS2$species_label)),
    #                    labels = levels(DvS2$species_label)) +
    scale_y_continuous(breaks = 1:length(levels(diff.dat$species_label)),
                       labels = levels(diff.dat$species_label)) +
    scale_color_viridis(discrete = TRUE, name = "Domain", option = "B", begin = 0.2, end = 0.6) +
    scale_fill_identity() + 
    labs(x = NULL, y = NULL, title = fig.title) +
    theme_bw() + 
    theme(axis.text.y = element_text(size = 10, face = "italic"),
          legend.position = "inside",
          legend.position.inside = c(if(lab.right) 0.8 else 0.2, 0.95),
          legend.background = element_rect(fill = "white", color = "gray80"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          title = element_text(size = 10),
          legend.key.size = unit(0.4, "cm"))  + 
    {if(remove.y) {
      theme(axis.text.y = element_blank())
    }}
}

























### Plot of domain estimates for Growth or Mortality ###
## Used in pair.plts.fcn below.
domain.grid.plt.fcn <- function(domains, domain.index, use.dat2, qt.max, q.g.p.labs) {
  ggplot(data = use.dat2 %>% filter(Domain %in% domain.index), aes(factor(Domain), Means, fill = factor(Domain))) + 
    geom_col() + 
    geom_errorbar(aes(ymax = UCI.95, ymin = LCI.95), width = 0.1) + 
    scale_fill_manual(values = virid.use[domain.index]) + 
    scale_y_continuous(limits = c(0, qt.max)) +
    scale_x_discrete(labels = q.g.p.labs[domain.index]) +
    theme_bw() +
    theme(legend.position = "none") + 
    labs(x = NULL, y = ylabs) +
    theme(text = element_text(size = 7))
}




### - Domain distribution plot, used in multiple analysis pathways - ###
#       Obtain figure of climate domains
# plot.spp = species of interest (e.g., "X11"), domain.matrix = domains (1 thru 9)
# for each species, var1, var.delt = carried through from input for pairs.plts.fcn,
# quant.lims, domain.n = quantile boundaries and number of plots, established in parse.tree.clim.fcn,
# size.trees = one of "small diameter ", "large diameter ", and (for pathways 1 and 3) "".  
domain.dist.plt.fcn <- function(plot.spp, domain.matrix, var1, var.delt, quant.lims, n.plots.used, size.trees) {
  
  q_plot_spp <- domain.matrix %>% dplyr::select(all_of(plot.spp), all_of(var1), all_of(var.delt)) %>%
    filter(get(plot.spp) > 0) 
  names(q_plot_spp)[1] <- "Domain"
  
  
  quant.lims.plt <- get(plot.spp, quant.lims)
  quant.lims.delt.plt <- get(plot.spp, quant.lims.delta)
  
  sppnum <- as.numeric(gsub("X", "", plot.spp))
  
  # Common and Genus/species name for plot title
  com.name <- spp.names$COMMON_NAME[spp.names$SPCD == sppnum]
  g.s.name <- paste(spp.names$GENUS[spp.names$SPCD == sppnum], spp.names$SPECIES[spp.names$SPCD == sppnum])
  
  # Need to isolate the correct colors if some quadrants are without plots:
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_domain)  # Get colors for plotting
  scatter.virid.use <- virid.use[n.plots.used$loc[n.plots.used$n > 0]]
  
  fig.lab <- paste0("Plot conditions for ", size.trees, com.name, ", ", g.s.name, ", ",  plot.spp)
  
  plot.vals.plt <- ggplot(q_plot_spp, aes(get(var1), get(var.delt), color = factor(Domain))) + 
    geom_point() + 
    #stat_ellipse(type = "norm", level = 0.95, col = "orange", lwd = 2) +
    geom_hline(yintercept = quant.lims.delt.plt[1], col = "blue") +
    geom_hline(yintercept = quant.lims.delt.plt[2], col = "blue") +
    geom_vline(xintercept = quant.lims.plt[1], col = "green") +
    geom_vline(xintercept = quant.lims.plt[2], col = "green") +
    geom_hline(yintercept = 0, col = "black") +
    theme_bw() + 
    scale_color_manual(values = scatter.virid.use, name = "Domains", labels = DOMAIN.LEVELS[which(n.plots.used$n > 0)]) + 
    labs(title = fig.lab,
         x = var.label, y = var.delt.label) +
    theme(text = element_text(size = 7),
          legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
  
  return(plot.vals.plt)
}









### -- Mapping elements -- ###

domain.map.fcn <- function(domain.matrix, spp.num, n.plots.used, virid.use) {
  
  map.dat.1 <- domain.matrix %>% 
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
    #                  labels = gsub(".", " ", DOMAIN.LEVELS[domains.used], fixed = TRUE)) +
    #scale_fill_viridis(option = "H", begin = virid.min, end = virid.max) +
    geom_polygon(data = west_df, mapping = aes(x = long, y = lat, group = group), color = "black", fill = "transparent") +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5),
          legend.position = "none")
}







# -- Plots for ANALYSIS.PATHWAY 1 -- #





## Function for plotting the mortality numbers by quantile and the distribution of plots in quantiles
pair.plts.fcn <- function(sppnum.to.plot, use.dat, domain.matrix,
                          quant.lims, domain.n, k){
  
  use.dat2 <- use.dat %>% filter(Species == sppnum.to.plot)
  
  
  plot.domain.dat <- use.dat %>% 
    filter(Species == sppnum.to.plot) %>%
    left_join(domain.level.table, by = c("Domain" = "q.num"))
  
  # Details for plotting
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_domain)  # Get colors for plotting
  
  qt.max <- ceiling(max(use.dat2$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
  
  
  # Create a grid of plots to match bivariate plot
  q.g.p.labs <- paste0(DOMAIN.LEVELS, ", n = ", as.numeric(use.dat2$n.plts)) # Labels for the domain grid plot
  
  ylabs <- if (k == 1) "Mean Growth (inches^2)" else "Mean Annual Mortality Rate"
  
  
  # These will be arranged in a plot below
  if(n_domain == 9){
    p1 <- domain.grid.plt.fcn(c("LiHc", "MiHc", "HiHc"), 1:3, use.dat2, qt.max, q.g.p.labs) 
    p2 <- domain.grid.plt.fcn(c("LiMc", "MiMc", "HiMc"), 4:6, use.dat2, qt.max, q.g.p.labs) 
    p3 <- domain.grid.plt.fcn(c("LiLc", "MiLc", "HiLc"), 7:9, use.dat2, qt.max, q.g.p.labs) 
    p_all <- plot_grid(p1, p2, p3, ncol = 1)
  } else if(n_domain == 6) {
    p1 <- domain.grid.plt.fcn(c("DL", "DM", "DH"), 1:3, use.dat2, qt.max, q.g.p.labs) 
    p2 <- domain.grid.plt.fcn(c("SL", "SM", "SH"), 4:6, use.dat2, qt.max, q.g.p.labs) 
    p_all <- plot_grid(p1, p2, ncol = 1)
  }
  
  
  # Code to prep for next two figures - helps in reducing the color palette.
  n_plots <- get(sppnum.to.plot, domain.n)
  n_plots2 <- tibble(loc = 1:n_domain, n = n_plots) %>%
    left_join(tibble(Domain = plot.domain.dat$Domains, loc = 1:n_domain), by = "loc")
  
  # Plotting the scatterplot of points in quadrants:
  plot.vals.plt <- domain.dist.plt.fcn(plot.spp = sppnum.to.plot, domain.matrix = domain.matrix, 
                                       var1 = var1, var.delt = var.delt, quant.lims = quant.lims,
                                       n.plots.used = n_plots2, size.trees = "") 
  
  # Plotting the map of plot locations:
  domain.map <- domain.map.fcn(domain.matrix, spp.num = sppnum.to.plot, n.plots.used = n_plots2, virid.use = virid.use) 
  
  if(n_domain == 9) {
    comb.plt <- ((p1/p2/p3 + plot_layout(axis_titles = "collect")) | plot.vals.plt | domain.map) #/ guide_area() + plot_layout(guides = 'collect', heights = c(10, 0.01)) 
  } else if(n_domain == 6) {
    comb.plt <- ((p1/p2 + plot_layout(axis_titles = "collect")) | plot.vals.plt | domain.map) #/ guide_area() + plot_layout(guides = 'collect', heights = c(10, 0.01)) 
  }
  
  ggsave(paste0(save.loc.fcn(k), sppnum.to.plot, "_plots.png"), 
         comb.plt, device = "png", width = 7, height = 5, units = "in")
}





#### Manuscript RMarkdown functions ------------------------------

# Prepping data for the state mortality and growth plots
state.table.fcn <- function(state.dat, all.dat, sig.places) {
  s.1 <- dplyr::left_join(spp.names.use, state.dat, by = c("SpeciesCode" = "Species")) %>%
    dplyr::filter(is.na(n.plts) == FALSE) %>%
    dplyr::mutate(across(where(is.numeric), ~ sprintf(paste0("%.", sig.places, "f"), .x)),
                  n.plts = as.integer(n.plts),
                  results = case_when(
                    is.na(Median) ~ "", 
                    grepl("NA", paste0(Means, " (", LCI.95, ",", UCI.95, "; ", n.plts, ")")) ~ "",
                    !is.na(Median) ~ paste0(Means, " (", LCI.95, ",", UCI.95, "; ", n.plts, ")")
                  ))%>%
    dplyr::mutate(results = case_when(is.na(Median) ~ "", 
                                      Median == "NA" ~ "",
                                      !is.na(Median) ~ paste0(Means, " (", LCI.95, ",", UCI.95, "; ", n.plts, ")"))) %>%
    dplyr::select(sci_name, State, results) %>%
    tidyr::pivot_wider(names_from = State, values_from = results) %>%
    dplyr::rename("California" = "CA", "Washington" = "WA", "Oregon" = "OR") %>%
    relocate("Oregon", .after = "California")# %>%
  
  s.2 <- dplyr::left_join(spp.names.use, all.dat, by = c("SpeciesCode" = "Species")) %>%
    dplyr::filter(is.na(n.plts) == FALSE) %>%
    dplyr::mutate(across(where(is.numeric), ~ sprintf(paste0("%.", sig.places, "f"), .x)),
                  n.plts = as.integer(n.plts),
                  results = case_when(
                    is.na(Median) ~ "", 
                    grepl("NA", paste0(Means, " (", LCI.95, ",", UCI.95, "; ", n.plts, ")")) ~ "",
                    !is.na(Median) ~ paste0(Means, " (", LCI.95, ",", UCI.95, "; ", n.plts, ")")
                  )) %>%
    dplyr::select(sci_name, results) %>%
    dplyr::rename("All" = "results")
  
  
  s.3 <- dplyr::left_join(s.1, s.2, by = "sci_name") %>%
    rename("Species Name" = "sci_name")
  #tidyr::replace_na(list(California = "", Washington = "", Oregon = "")) %>%
  #dplyr::filter(All != "", !grep("NA", All))
  
  return(s.3)
}


# This function allows tables and species to be searched for maximum and minimum
#  mean values and report on the species, mean, and 95% CI.
maxmin.fcn <- function(table1, mean_col, lower_ci_col, upper_ci_col, 
                       species_col = "Species", extreme = "max", 
                       digits = 1,
                       cm2 = FALSE) {   # Should results be transformed to cm2 (growth)?
  
  table2 <- table1 %>% left_join(spp.names.use, by = c("Species" = "SpeciesCode"))
  
  required_cols <- c(mean_col, lower_ci_col, upper_ci_col, species_col)
  missing_cols <- required_cols[!required_cols %in% names(table2)]
  
  # Remove rows with missing values in key columns
  complete_data <- table2[complete.cases(table2[, required_cols]), ]
  
  if (nrow(complete_data) == 0) {
    stop("No complete cases found in the data")
  }
  
  # Find the species with extreme value
  mean_values <- complete_data[[mean_col]]
  
  if (extreme == "max") {
    target_idx <- which.max(mean_values)
  } else {
    target_idx <- which.min(mean_values)
  }
  
  # Extract the values for the target species
  target_species <- complete_data$sci_name[target_idx]
  target_mean <- complete_data[[mean_col]][target_idx]
  target_lower <- complete_data[[lower_ci_col]][target_idx]
  target_upper <- complete_data[[upper_ci_col]][target_idx]
  
  if(cm2 == FALSE) {
    cm2.use <- 1
  } else {
    cm2.use <- 6.4516
  }
  
  # Format the result
  maxmin_output <- list(spp = target_species,
                        mean = round(target_mean * cm2.use, digits),
                        mean.ci = sprintf("%.*f (%.*f, %.*f)", 
                                          digits, target_mean * cm2.use,
                                          digits, target_lower * cm2.use, 
                                          digits, target_upper * cm2.use),
                        # Add this new element:
                        mean.ci.spp = sprintf("%.*f (%.*f, %.*f; %s)", 
                                              digits, target_mean * cm2.use,
                                              digits, target_lower * cm2.use, 
                                              digits, target_upper * cm2.use,
                                              target_species),
                        mean.spp = sprintf("%.*f (%s)",
                                           digits, target_mean * cm2.use,
                                           target_species)
  )
  
  return(maxmin_output)
}



### Function for preparing supplementary tables of CWD max/min/quantiles.
###  Used in Manuscript_information.R
########################

summ.spp.fcn <- function(data.all, var.select, spp.names.use) {

summ.spp <- all_dat %>% left_join(climate.use, by = "puid") %>%
  dplyr::select(-STATECD, -ESTN_UNIT, -STRATUMCD, -w, -stratum, -n_h.plts, -LAT, -LON) %>%
  pivot_longer(cols = starts_with("nd."), names_to = "spp", values_to = "values") %>%
  mutate(spp = as.numeric(gsub("nd.", "", spp))) %>%
  filter(values > 0) %>%
  group_by(spp) %>%
  summarize(
    n = n(),
    Minimum = min(get(var.select)),
    q05 = quantile(get(var.select), 0.05),
    q25 = quantile(get(var.select), 0.25),
    Mean = mean(get(var.select)),
    Median = median(get(var.select)),
    q75 = quantile(get(var.select), 0.75),
    q95 = quantile(get(var.select), 0.95),
    Maximum = max(get(var.select))
  ) %>%
  left_join(spp.names.use %>% dplyr::select(SPCD, sci_name), by = c("spp" = "SPCD")) %>%
  relocate(sci_name) %>%
  rename("Species" = "sci_name") %>%
  dplyr::select(-spp) %>%
  mutate(across(c(-Species, -n), round, digits = 1))

return(summ.spp)
}


