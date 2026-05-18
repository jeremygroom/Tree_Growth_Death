### ==> Initial DBH: state-level estimates and multi-species panel plot ========
###  Manuscript revision context: characterizes initial tree size (PREVDIA) by
###  state to address reviewer concerns that BAI growth differences across CWD
###  domains may reflect differences in initial size. Same design-based ratio
###  estimator as the BAI analysis; same survivor cohort (STATUSCD == 1).


  # Build the initial-DBH arrays from the existing growth-tree data prep.
  # This relies on clim.growth.resp.fcn now returning init.dbh.val.
  arrays.init.dbh <- parse.tree.clim.fcn(
    tree.grow.dat,
    clim.var      = CLIM.VAR.USE,
    analysis.type = "init.dbh",
    resp.dat      = "init.dbh.val",
    tot.dat       = "growth.n.trees",
    selected.spp  = SEL.SPP,
    clim.dat      = climate.use
  )
  
  vals_dat.dbh    <- arrays.init.dbh$vals_dat
  all_dat.dbh     <- arrays.init.dbh$all_dat
  state.array.dbh <- arrays.init.dbh$state.array
  state.n.dbh     <- arrays.init.dbh$state.n
  all.array.dbh   <- arrays.init.dbh$all.array
  all.n.dbh       <- arrays.init.dbh$all.n
  
  strata.num.dbh <- vals_dat.dbh %>% dplyr::select(stratum, puid) %>%
    ungroup() %>% mutate(row.id = row_number())
  
  plan(multisession, workers = n.cores)
  
  state.list         <- unique(vals_dat.dbh$STATECD)  # 6=CA, 41=OR, 53=WA
  state.domain.index <- 1:length(state.list)
  
  # State-level bootstrap (analogous to state.bootstrap_results in Overall_Mort_Est.R)
  state.dbh.bs <- generate_bootstrap_array.fcn(
    vals.dat     = vals_dat.dbh,
    all.dat      = all_dat.dbh,
    domain.array = state.array.dbh,
    domain.n     = state.n.dbh,
    selected.spp = SEL.SPP,
    n_iter       = BS.N,
    strata.num   = strata.num.dbh,
    PlotDat      = PlotDat
  )
  
  state.init.dbh.summaries <- state.domain.index %>%
    purrr::map(\(d) domain.sum.fcn(state.dbh.bs, d, domain_n = state.n.dbh)) %>%
    do.call(rbind, .) %>%
    mutate(State = case_match(Domain, 1 ~ "CA", 2 ~ "WA", 3 ~ "OR")) %>%
    arrange(Species, Domain)
  
  # Three-state combined estimate (parallel to all.species.summaries)
  all.dbh.bs <- generate_bootstrap_array.fcn(
    vals.dat     = vals_dat.dbh,
    all.dat      = all_dat.dbh,
    domain.array = all.array.dbh,
    domain.n     = all.n.dbh,
    selected.spp = SEL.SPP,
    n_iter       = BS.N,
    strata.num   = strata.num.dbh,
    PlotDat      = PlotDat
  )
  all.init.dbh.summaries <- domain.sum.fcn(all.dbh.bs, 1, domain_n = all.n.dbh)
  
  # Save tables. Inches is the native FIA unit for PREVDIA; cm versions are
  # provided for direct comparison with the cm^2 BAI manuscript figures.
  init.dbh.outloc <- save.loc.fcn(1)   # alongside Growth outputs
  
  write_csv(state.init.dbh.summaries,
            file = paste0(init.dbh.outloc, "State_Initial_DBH_inches_", CLIM.VAR.USE, ".csv"))
  write_csv(state.init.dbh.summaries %>%
              mutate(across(c(Means, Median, LCI.95, UCI.95), ~ .x * 2.54)),
            file = paste0(init.dbh.outloc, "State_Initial_DBH_cm_",     CLIM.VAR.USE, ".csv"))
  write_csv(all.init.dbh.summaries,
            file = paste0(init.dbh.outloc, "All_Initial_DBH_inches_",   CLIM.VAR.USE, ".csv"))
  write_csv(all.init.dbh.summaries %>%
              mutate(across(c(Means, Median, LCI.95, UCI.95), ~ .x * 2.54)),
            file = paste0(init.dbh.outloc, "All_Initial_DBH_cm_",       CLIM.VAR.USE, ".csv"))
  
  # ----- Multi-species panel plot -----
  # spp.names.fig is otherwise built inside the n_domain==6 block further down,
  # so build a local copy here to keep this section self-contained.
  spp.names.fig.dbh <- spp.names %>% dplyr::select(SPCD, GENUS, SPECIES) %>%
    mutate(Species = paste0("X", SPCD)) %>%
    dplyr::select(-SPCD)
  
  init.dbh.panel <- state.var.panel.fcn(
    state.dat         = state.init.dbh.summaries,
    x.label           = "Mean initial DBH (cm)",
    fig.title         = NULL,
    scale.factor      = 2.54,
    spp.names.fig.use = spp.names.fig.dbh
  )
  
  ggsave(paste0(init.dbh.outloc, "Initial_DBH_State_Panel.jpeg"),
         plot = init.dbh.panel, device = ragg::agg_jpeg,
         width = 7, height = 8, units = "in", res = 1000)
  
