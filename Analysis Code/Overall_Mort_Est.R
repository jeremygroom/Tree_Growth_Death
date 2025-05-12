### Code for finding the overall mortality for each species ###
## It doesn't matter what the climate variable is as the deepest level = state.###

tic()

# First creating the data sets, then joining the results in arrays
t.mort.all.vpdmin1 <- map(as.numeric(gsub("X", "", spp.id)), clim.mort.resp.fcn, clim.var = "vpdmin", treedat.sel = treedat.use)
a.mort.all.vpdmin1 <- parse.tree.clim.fcn(t.mort.all.vpdmin1, "vpdmin", analysis.type = "mort", resp.dat = "died.out", tot.dat = "all.trees", selected.spp = spp.id)


      # Now simplifying names
      vals_dat <- a.mort.all.vpdmin1$vals_dat
      all_dat <- a.mort.all.vpdmin1$all_dat
      quant.array <- a.mort.all.vpdmin1$quant.array
      state.array <- a.mort.all.vpdmin1$state.array
      state.n <- a.mort.all.vpdmin1$state.n
      all.array <- a.mort.all.vpdmin1$all.array
      all.n <- a.mort.all.vpdmin1$all.n
      quant.matrix <- a.mort.all.vpdmin1$quant.matrix
      centroid.array <- a.mort.all.vpdmin1$centroid.array
      quant.n <- a.mort.all.vpdmin1$quant.n 
      quant.lims <- a.mort.all.vpdmin1$quant.lims
      quant.lims.delta <- a.mort.all.vpdmin1$quant.lims.delt
      
      
      
      # Setting up parallel computing
      plan(multisession, workers = n.cores) 
      
      
      # State-level mortality estimates.
        state.list <- unique(vals_dat$STATECD) # 6 = California, 41 = Oregon, 53 = Washington
        state.ests <- map(1:length(state.list), q_mort.grow.fcn, vals.dat = vals_dat, 
                          all.dat = all_dat, array.name = state.array, category.n = state.n,
                          selected.spp = spp.id)
        state.ests2 <- rbindlist(state.ests) %>% data.frame() %>%  # data.table function
          mutate(State = case_match(Quantile, 1 ~ "CA", 2 ~ "OR", 3 ~ "WA")) %>%
          select(-Quantile)
        
          # Writing state values to file
        fwrite(state.ests2, file = paste0(RESULTS1.LOC, "Mort_figs_", CLIM.VAR[1], "/State_Mort_Ests_", CLIM.VAR[1], ".csv")) 
        
      # All states combined mortality
        all.mort.ests <- q_mort.grow.fcn(q1 = 1, vals.dat = vals_dat, all.dat = all_dat, 
                        array.name = all.array, category.n = all.n,
                        selected.spp = spp.id)
        # Writing values for species across states to another file
        fwrite(all.mort.ests, file = paste0(RESULTS1.LOC, "Mort_figs_", CLIM.VAR[1], "/All_Mort_Ests_", CLIM.VAR[1], ".csv")) 
toc()
        