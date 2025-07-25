### Code for finding the state level and overall growth/mortality for each species ###
## It doesn't matter what the climate variable is as the deepest level = state.###

for(k in 1:2){ # 1 = growth, 2 = mortality
  
  g.m.labels <- c("Growth", "Mortality")
  
  
  if(k == 1) {
    t.all <- map(as.numeric(gsub("X", "", SEL.SPP)), clim.growth.resp.fcn, clim.var = CLIM.VAR.USE, treedat.sel = treedat.use, clim.dat = climate.use)
    a.all <- parse.tree.clim.fcn(t.all, clim.var = CLIM.VAR.USE, analysis.type = "grow", resp.dat = "growth.val", tot.dat = "growth.n.trees", selected.spp = SEL.SPP, clim.dat = climate.use)
  } else { 
    t.all <- map(as.numeric(gsub("X", "", SEL.SPP)), clim.mort.resp.fcn, clim.var = CLIM.VAR.USE, treedat.sel = treedat.use, clim.dat = climate.use)
    a.all <- parse.tree.clim.fcn(t.all, clim.var = CLIM.VAR.USE, analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = SEL.SPP, clim.dat = climate.use)
  }
  
  
  
  # Now simplifying names
  vals_dat <- a.all$vals_dat
  all_dat <- a.all$all_dat
  domain.array <- a.all$domain.array
  state.array <- a.all$state.array
  state.n <- a.all$state.n
  all.array <- a.all$all.array
  all.n <- a.all$all.n
  centroid.array <- a.all$centroid.array
  domain.n <- a.all$domain.n 
  quant.lims <- a.all$quant.lims
  quant.lims.delta <- a.all$quant.lims.delt
  
  
  strata <- unique(vals_dat$stratum)
  strata.num <- vals_dat %>% dplyr::select(stratum, puid) %>%
    ungroup() %>% 
    mutate(row.id = row_number())
  
  
  
  # Setting up parallel computing
  plan(multisession, workers = n.cores) 
  
  # State-level mortality estimates.
  state.list <- unique(vals_dat$STATECD) # 6 = California, 41 = Oregon, 53 = Washington
  state.domain.index <- 1:length(state.list)
  
  ## Generating the bootstrap values: 
  state.bootstrap_results <- generate_bootstrap_array.fcn(
    vals.dat = vals_dat,
    all.dat = all_dat,
    domain.array = state.array,
    domain.n = state.n,
    selected.spp = SEL.SPP,
    n_iter = BS.N
  )
  
  # Finding and saving domain summaries
  state.domain.summaries <- state.domain.index %>% 
    map(\(d) domain.sum.fcn(state.bootstrap_results, d)) %>%
    do.call(rbind, .) %>%
    mutate(State = case_match(Domain, 1 ~ "CA", 2 ~ "WA", 3 ~ "OR")) %>%
    arrange(Species, Domain)
  
  write_csv(state.domain.summaries, file = paste0(save.loc.fcn(k), "State_", g.m.labels[k], "_", CLIM.VAR.USE, ".csv"))

  
  # Combined estimates of growth and mortality across OR/WA/CA
  all.bootstrap_results <- generate_bootstrap_array.fcn(
    vals.dat = vals_dat,
    all.dat = all_dat,
    domain.array = all.array,
    domain.n = all.n,
    selected.spp = SEL.SPP,
    n_iter = BS.N
  )
  
  all.species.summaries <- domain.sum.fcn(all.bootstrap_results, 1)
  
  # Writing values for species across states to another file
  write_csv(all.species.summaries, file = paste0(save.loc.fcn(k), "All_", g.m.labels[k], "_", CLIM.VAR.USE, ".csv"))
  
}        