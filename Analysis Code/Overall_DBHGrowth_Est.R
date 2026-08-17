# Obtaining growth in DBH values, similar to Overall_Mort_est.R

### ==> Initial DBH: state-level estimates and multi-species panel plot ========
###  Manuscript revision context: characterizes initial tree size (PREVDIA) by
###  state to address reviewer concerns that BAI growth differences across CWD
###  domains may reflect differences in initial size. Same design-based ratio
###  estimator as the BAI analysis; same survivor cohort (STATUSCD == 1).


# Build the initial-DBH arrays from the existing growth-tree data prep.

vals_dat.dbh    <- arrays.dbh.grow$vals_dat
all_dat.dbh     <- arrays.dbh.grow$all_dat
state.array.dbh <- arrays.dbh.grow$state.array
state.n.dbh     <- arrays.dbh.grow$state.n
all.array.dbh   <- arrays.dbh.grow$all.array
all.n.dbh       <- arrays.dbh.grow$all.n

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

state.dbh.summaries <- state.domain.index %>%
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
all.dbh.summaries <- domain.sum.fcn(all.dbh.bs, 1, domain_n = all.n.dbh)

# Save tables. Inches is the native FIA unit for PREVDIA; cm versions are
# provided for direct comparison with the cm^2 BAI manuscript figures.
dbh.outloc <- save.loc.fcn(1)   # alongside Growth outputs

write_csv(state.dbh.summaries %>%
            mutate(across(c(Means, Median, LCI.95, UCI.95), ~ .x * 2.54)),
          file = paste0(dbh.outloc, "State_DBH_cm_",     CLIM.VAR.USE, ".csv"))
write_csv(all.dbh.summaries %>%
            mutate(across(c(Means, Median, LCI.95, UCI.95), ~ .x * 2.54)),
          file = paste0(dbh.outloc, "All_DBH_cm_",       CLIM.VAR.USE, ".csv"))

