
##### ---- Functions ---- #####

### DATA PREP FUNCTIONS - MORTALITY ####


## Quantile masks for each species. Used in clim.tree.resp.fcn (called in Data_Prep.R), operates on a single quantile level and species 
quant.divide.fcn <- function(quant.level, spcd2, tree.plots2, PlotDat) {
  quant.name <- paste0(quant.level, ".",  spcd2)
  quant.died.out <- left_join(PlotDat, tree.plots2 %>% 
                                dplyr::select(State_Plot, var.deltvar) %>% # var.deltvar defined in clim.tree.resp.fcn, "LiHc" etc.
                                filter(var.deltvar == quant.level), by = "State_Plot") %>%
    mutate(var.deltvar = ifelse(is.na(var.deltvar), 0, 1)) %>%
    rename(!!quant.name := var.deltvar) 
}


## Function to determine the number of trees that died in a plot and their climate quantiles, used in Data_Prep.R.
clim.tree.resp.fcn <- function(spcd, clim.var) {
  var1 <- paste0("pre.", clim.var)
  var.delt <- paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  Tree2 <- Tree1.2 %>% filter(SPCD == spcd) %>%
    dplyr::select(State_Plot, SUBP, TREE, SDN)
  
  tree.plots <- clim.dat %>% filter(State_Plot %in% Tree2$State_Plot) %>% 
    dplyr::select(State_Plot, all_of(var1), all_of(var.delt))
  
  # Data quantiles
  quant.lims.delta.pos <- quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
  quant.lims.delta <- c(quant.lims.delta.neg, quant.lims.delta.pos)
  quant.lims <- quantile(get(var1, tree.plots), probs = QUANT.PROBS)
  
  # Applying quantiles to a particular species
  tree.plots2 <- tree.plots %>% dplyr::select(State_Plot, all_of(c(var1, var.delt))) %>% 
    mutate(var.quants = ifelse(get(var1) > quant.lims[2], 1, ifelse(get(var1) < quant.lims[1], -1, 0)),
           deltvar.quants = ifelse(get(var.delt) > quant.lims.delta[2], 1, ifelse(get(var.delt) < quant.lims.delta[1], -1, 0)),
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
  
  
  quant.out1 <- map(QUANT.LEVELS, quant.divide.fcn, spcd2 = spcd, tree.plots2 = tree.plots2, PlotDat = PlotDat)
  quant.out2 <- reduce(quant.out1, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM")) %>%
    dplyr::select(-c("STATECD", "PLOT_FIADB", "State_Plot", "STRATUM", "W_h"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, died.out = died.out, all.trees = all.trees, 
                   quant.out = quant.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}



### ------ Mortality estimation -------  ####


# Find the estimates for mu_hat (the ratio estimator), Y_hat (numerator, thing we're interested in), and 
#    X_Hat (denominator, density of things like total trees of a species). See Equation 6 and related.
mean.q.fcn <- function(dat_tot, resp, spp.sel1) {     #dat_tot = data frame of species total for a quantile (denominator). resp = response data frame (numerator) for quantile.
  
  Zt <- Yv <- R <- 0 #  rep(0,length(sppname))    # Mean difference in response for a given species: Zt = Mean value for total trees, 
  #Yv = mean value for/of trees of interest (mean # that died), R = ratio of the two.
  #for (i in 1:length(sppname)) {                         # i indexes species
  col.spp <- grep(paste0("nd.", spp.sel1),  colnames(dat_tot))
  col.name <- paste0("nd.", spp.sel1)
  
  strat.vals <- table(dat_tot$STRATUM)/nrow(dat_tot)
  data.frame(strat.vals = strat.vals, w_h = unique(dat_tot$W_h)) %>% mutate(diff = strat.vals - w_h, pct.diff = diff/strat.vals)
  
  
  Zt_i <- Yv_i <- rep(0, strat.num)  # Weighted values for each strata
  for (h in 1:strat.num) {                  # h indexes strata
    n_h <- length(dat_tot$STRATUM[dat_tot$STRATUM == strat2$STRATUM[h]] )    # Number of plots in stratum h 
    Zt_i[h] <- strat2$W_h[h] * sum(get(col.name, dat_tot)[dat_tot$STRATUM == strat2$STRATUM[h]]) / n_h     
    Yv_i[h] <- strat2$W_h[h] * sum(get(col.name, resp)[resp$STRATUM == strat2$STRATUM[h]])/ n_h     
  }
  Zt <- sum(Zt_i) 
  Yv <- sum(Yv_i) 
  R <- Yv / Zt
  
  means <- list(Zt = Zt, Yv = Yv, R = R)
  return(means)
  
  #    dat.use <- left_join(dat_tot %>% select(STATECD, PLOT_FIADB, STRATUM, W_h, eval(col.name)),
  #                       resp %>% select(STATECD, PLOT_FIADB, STRATUM, W_h, eval(col.name)), 
  #                       by = c("STRATUM", "W_h", "PLOT_FIADB", "STATECD")) %>%
  #    rename("all" := (!!paste0(col.name, ".x")),
  #           "reduced" := (!!paste0(col.name, ".y"))) %>%
  #    group_by(STRATUM, W_h) %>%
  #    reframe(n_h = n(),
  #            Zt_i = W_h * sum(all)/n_h,
  #            Yv_i = W_h * sum(reduced)/n_h) %>% distinct()
  #  
  #  results.mean <- dat.use %>% group_by() %>%
  #    reframe(Zt = sum(Zt_i),
  #            Yv = sum(Yv_i),
  #            R = Yv/Zt)
  
}

# Helper functions for ratio.SE.fcn: 
# Function that finds variance or covariance for individual strata in the function ratio.SE.fcn
varcov.fcn <- function(xy.sum, xsum, ysum, nnh) {(1/(nnh * (nnh - 1))) * (xy.sum - (1/nnh) * xsum * ysum)} # Equations 10 and 11

# Function that finds the overall weighted variance and covariance in the function ratio.SE.fcn
wt.varcov.fcn <- function(vc, nnh, nn) {(1/nn) * (sum(strat2$W_h * nnh * vc) + sum((1 - strat2$W_h) * (nnh/nn) * vc))} # Equations 8 and 9


# Find the standard error of the variance for a species' quartile
ratio.SE.fcn <- function(d.all_z, d.ado_y, meandat)	{   # ii = column for response variable, meandat = list generated by actmean function
  
  col.spp <- which(colnames(d.all_z) == paste0("nd.", spp.sel)) # ID the column for the species in question
  
  Zt <- meandat$Zt   # Output from mean.q.fcn
  mu <- meandat$R    # Output from mean.q.fcn
  
  
  # for (i in 1:length(sppname))	{
  Zh <- Yh <- nn_h <- Zu2h <- Yu2h <- ZYh <- rep(0, length(strat2$STRATUM))  # Weighted values for each strata
  for(h in 1:length(strat2$STRATUM)){ 
    # Pieces for calculating stratum-level variances.
    # Number of plots in each stratum:
    nn_h[h] <- length(d.ado_y$STRATUM[d.ado_y$STRATUM == strat2$STRATUM[h]])
    # For each stratum h for each species, first calculate the Z and Y values.
    Zh[h] <- sum(d.all_z[d.all_z$STRATUM == strat2$STRATUM[h], col.spp])
    Yh[h] <- sum(d.ado_y[d.ado_y$STRATUM == strat2$STRATUM[h], col.spp])
    # Square values for Z and Y (for Step 6, 7, and 8)
    Zu2h[h] <- sum(d.all_z[d.all_z$STRATUM == strat2$STRATUM[h], col.spp] ^ 2)
    Yu2h[h] <- sum(d.ado_y[d.ado_y$STRATUM == strat2$STRATUM[h], col.spp] ^ 2)
    # ZY cross-products between first or second visits(for steps 9, 10)
    ZYh[h] <- sum(d.all_z[d.all_z$STRATUM == strat2$STRATUM[h], col.spp] * d.ado_y[d.ado_y$STRATUM == strat2$STRATUM[h], col.spp])
    
  }
  
  (7121 - 592*(0.00571^2))/(592 * (592 - 1))
  (Zu2h[18] - (1/nn_h[18])*Zh[18]^2)/(nn_h[18] * (nn_h[18] - 1))
  
  varcov.fcn <- function(xy.sum, xsum, ysum, nnh) {(1/(nnh * (nnh - 1))) * (xy.sum - (1/nnh) * xsum * ysum)} # Equations 10 and 11
  
  nn <- dim(d.all_z)[1]   # nn = N = total plots including ones not sampled
  
  # Weighed variances (Step 8 of algorithm sheet)
  varZ <- varcov.fcn(Zu2h, Zh, Zh, nn_h)
  varY <- varcov.fcn(Yu2h, Yh, Yh, nn_h)
  covZY <- varcov.fcn(ZYh, Zh, Yh, nn_h)
  
  w.varZ <- wt.varcov.fcn(varZ, nn_h, nn)
  w.varY <- wt.varcov.fcn(varY, nn_h, nn)
  
  # Covariances for visit-pair estimator (Step 11)
  w.covZY <- wt.varcov.fcn(covZY, nn_h, nn)
  
  # Variance of ratios (Step 12)
  
  varR <- (1/Zt^2) * (w.varY + mu^2 * w.varZ - 2 * mu * w.covZY)  ## Clean up for production.  Can verify using monte carlo
  
  varR
}



