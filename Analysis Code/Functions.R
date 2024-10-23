
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



#### Data prep function for mortality --------------------------------------------------

## Function to determine the number of trees that died in a plot and their climate quantiles, used in Data_Prep.R.
clim.mort.resp.fcn <- function(spcd, clim.var) {
  var1 <- paste0("pre.", clim.var)
  var.delt <- paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  # For alive/dead tree assessments
  Tree2 <- Tree1.2 %>% filter(SPCD == spcd) %>%
    dplyr::select(State_Plot, SUBP, TREE, SDN)
  
  tree.plots <- clim.dat %>% filter(State_Plot %in% Tree2$State_Plot) %>% 
    dplyr::select(State_Plot, all_of(var1), all_of(var.delt))
  
  # Data quantiles
    # Setting the delta-var amounts to be fixed values.  See Global for setting VAR.DELTA.BOUNDARIES.
  quant.lims.delta.pos <- VAR.DELTA.BOUNDARIES$max.min[VAR.DELTA.BOUNDARIES$clim.var == clim.var] #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- -quant.lims.delta.pos #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
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

  # Mortality data
  Tree3 <- left_join(Tree2, tree.plots2, by = "State_Plot") %>%
    group_by(State_Plot) %>% #, var.deltvar, get(var1), get(var.delt)) %>%  # This will assign the names `get(var1)` and `get(var.delt)`, which are cleaned up below
    reframe(n.trees = n(),
            n.died = length(SDN[SDN == 2]),
            pct.died = n.died / n.trees)# %>%
    #rename(!!var1 := `get(var1)`,
     #      !!var.delt := `get(var.delt)`) 

       # Separate out a data set for those that died
  died.out <- left_join(PlotDat, Tree3 %>% dplyr::select(State_Plot, n.died), by = "State_Plot") %>%
    rename(!!return.col := n.died)
  died.out[is.na(died.out)] <- 0
  
       # Data set for all trees, alive and died
  all.trees <- left_join(PlotDat, Tree3 %>% dplyr::select(State_Plot, n.trees), by = "State_Plot") %>%
    rename(!!return.col := n.trees)
  all.trees[is.na(all.trees)] <- 0
  
  
  quant.out1 <- map(QUANT.LEVELS, quant.divide.fcn, spcd2 = spcd, tree.plots2 = tree.plots2, PlotDat = PlotDat)
  quant.out2 <- reduce(quant.out1, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM", "SITECLCD")) %>%
    dplyr::select(-c("STATECD", "PLOT_FIADB", "State_Plot", "STRATUM", "W_h", "SITECLCD"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, died.out = died.out, all.trees = all.trees, 
                   quant.out = quant.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}


#### Data prep function for growth  --------------------------------------------------

## Function to determine the number of trees that died in a plot and their climate quantiles, used in Data_Prep.R.
clim.growth.resp.fcn <- function(spcd, clim.var) {
  var1 <- paste0("pre.", clim.var)
  var.delt <- paste0("delt.", clim.var)
  return.col <- paste0("nd.", spcd)
  
  # For basal area growth/tree
  Treediam.1 <- Tree1.2 %>% filter(SPCD == spcd, STATUSCD == 1, INVYR < 2011) %>%
    select(State_Plot, SUBP, TREE, SPCD, DIA)
  
  Treediam.2 <-  Tree1.2 %>% filter(SPCD == spcd, STATUSCD == 1, INVYR > 2010) %>%
    select(State_Plot, SUBP, TREE, SPCD, DIA)
  
  TreeG <- inner_join(Treediam.1, Treediam.2, by = c("State_Plot",  "SUBP", "TREE", "SPCD")) %>% #The inner join reduces the two data frames to trees that occur alive at times 1 and 2.
    mutate(BA_Growth = pi * (DIA.y / 2)^2 - pi * (DIA.x / 2)^2) %>%
    group_by(State_Plot) %>%
    reframe(n_g.trees = n(),
            total.growth = sum(BA_Growth))  # Will separate the number of trees and the sum of the growth to obtain a weighted mean of growth per plot and number of trees per plot 
  
  tree.plots <- clim.dat %>% filter(State_Plot %in% TreeG$State_Plot) %>% 
    dplyr::select(State_Plot, all_of(var1), all_of(var.delt))
  
  # Data quantiles
  quant.lims.delta.pos <- VAR.DELTA.BOUNDARIES$max.min[VAR.DELTA.BOUNDARIES$clim.var == clim.var] #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) > 0], probs = QUANT.DELTA.PROBS.POS)
  quant.lims.delta.neg <- -quant.lims.delta.pos #quantile(get(var.delt, tree.plots)[get(var.delt, tree.plots) < 0], probs = QUANT.DELTA.PROBS.NEG)
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
  # Growth data  
  # Overall data:
  Growth.dat <- left_join(PlotDat, TreeG, by = "State_Plot") 
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
  quant.out2 <- reduce(quant.out1, left_join, by = c("STATECD", "PLOT_FIADB", "State_Plot", "W_h", "STRATUM", "SITECLCD")) %>%
    dplyr::select(-c("STATECD", "PLOT_FIADB", "State_Plot", "STRATUM", "W_h", "SITECLCD"))
  
  out.list <- list(spp.list = spp.list, spp.id = spp.id, died.out = died.out, all.trees = all.trees, growth.val = growth.val, growth.n.trees = growth.n.trees,
                   quant.out = quant.out2, quant.lims = quant.lims, quant.lims.delta = quant.lims.delta, cent.loc = dat.cent)
  
  return(out.list)
}




### ------ Mortality and growth estimation -------  ####

#  This function is nested within quant.est.spp.fcn, which is in turn nested within bs.fcn and finally q_mort.grow.fcn.  It is called in DeadTree_Analysis
# Find the estimates for mu_hat (the ratio estimator), Y_hat (numerator, thing we're interested in), and 
#    X_Hat (denominator, density of things like total trees of a species). See Equation 6 and related.
mean.q.fcn <- function(dat_tot, dat_vals, spp.sel1) {     #dat_tot = data frame of species total for a quantile (denominator). dat_vals = response data frame (numerator) for quantile.
  
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
    Yv_i[h] <- strat2$W_h[h] * sum(get(col.name, dat_vals)[dat_vals$STRATUM == strat2$STRATUM[h]]) / n_h     
  }
  Zt <- sum(Zt_i) 
  Yv <- sum(Yv_i) 
  R <- Yv / Zt
  
  means <- list(Zt = Zt, Yv = Yv, R = R)
  return(means)
  
  #    dat.use <- left_join(dat_tot %>% select(STATECD, PLOT_FIADB, STRATUM, W_h, eval(col.name)),
  #                       dat_vals %>% select(STATECD, PLOT_FIADB, STRATUM, W_h, eval(col.name)), 
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


### This function will find, for any species, the estimate for quantile q1. Nested within bs.fcn and q_mort.grow.fcn.
quant.est.spp.fcn <- function(spp.sel, d.all_use, d.vals_use, samp) {
  
  spp.id1 <- paste0("X", spp.sel)
  n_q <- get(spp.id1, category.n)[q1]   # Are there fewer than N.PLOT.LIM? If so, we can skip the calcs.
  
  # Means plus other metrics that are carried over into the calculation of the variance		
  if (n_q < N.PLOT.LIM) { # If too few points, just enter NA
    q_bs.mean <- NA 
  } else {
    q_bs.mean <- mean.q.fcn(dat_tot = d.all_use[samp, ], dat_vals = d.vals_use[samp, ], spp.sel1 = spp.sel)$R # Mean given the selected rows
  }
  return(q_bs.mean)
}


# Bootstrap function for obtaining an estimate of the mean for a single iteration.  Nested within q_mort.grow.fcn.
bs.fcn <- function(iter, d.all_use, d.vals_use, e=parent.frame()) { # Iteration number, all of the relevant trees data, the tree values data (number dead, basal area growth)
  
  # First, sample the tables in use. This procedure obtains samples with replacement from each stratum and then combines the selected
  #  row numbers.  
  samp <- unlist(map(strata, function(x) sample(strata.num$val[strata.num$stratum == x], replace = TRUE)))

  temp.spp <- rep(NA, length(SEL.SPP)) # Empty storage vector
  
  bs.out[iter, ] <- unlist(spp.list %>% map(\(x) quant.est.spp.fcn(x, d.all_use, d.vals_use, samp)))  # For a given iteration, finding the mean value for a given quantile (domain)
}

## Base function for finding the estimated means for mortality/growth for a given quantile.
q_mort.grow.fcn <- function(q1, vals.dat, all.dat, array.name, category.n) {
  # Multiplying all of the species dead values (or species total values) by the respective species masks for "LiLc".
  #   With this step complete, the values for that quantile can be calculated for all species.
  
  vals.use <- vals.dat %>% select(starts_with("nd.")) * array.name[, q1, ] 
  dat.vals.use <- bind_cols(orig[, 1:6], vals.use) 
  
  all.use <- all.dat %>% select(starts_with("nd.")) * array.name[, q1, ] 
  dat.all.use <- bind_cols(orig[, 1:6], all.use) 
  
  # Bootstrap estimate holding matrix
  bs.out <- data.frame(matrix(NA, nrow = BS.N, ncol = length(spp.list)))
  
  # Running the parallel portion of the function
  furrr.out <- 1:BS.N %>% future_map(\(x) bs.fcn(iter = x, d.all_use = dat.all.use, d.vals_use = dat.vals.use), .options = furrr_options(seed = TRUE))
  
  # Processing the results and returning quantiles/mean
  furrr.table <- matrix(unlist(furrr.out), ncol = length(spp.list), byrow = TRUE)
  quants <- data.frame(t(apply(furrr.table, 2, function(x) quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE))))
  names(quants) <- c("Median", "LCI.95", "UCI.95")
  quants$Means <- apply(furrr.table, 2, mean, na.rm = TRUE)
  quants$n.plts <- as.vector(category.n[q1, SEL.SPP ])
  quants$Species <- SEL.SPP
  quants$Quantile <- q1
  
  return(quants)
}
















# MOTHBALLED VARIANCE ESTIMATION FUNCTIONS
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














#### Graphical output functions ---------------------------------


## Function for plotting the mortality numbers by quantile and the distribution of plots in quantiles
pair.plts.fcn <- function(sppnum.to.plot){
  
  plot.spp <- sppnum.to.plot
  plot.quant.dat <- quant.table %>% 
    filter(Species == plot.spp) %>%
    #bind_cols(quant.n %>% dplyr::select(plot.spp)) %>%
    #rename("n" = plot.spp) %>%
    mutate(Quantiles = case_match(Quantile, 1:length(QUANT.LEVELS) ~ QUANT.LEVELS),
           Quantiles = factor(Quantiles, levels = QUANT.LEVELS))
  
  
  # Details for plotting
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_quant)  # Get colors for plotting
  
  qt.max <- ceiling(max(plot.quant.dat$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
  
  
  # Create a grid of plots to match bivariate plot
  q.g.p.labs <- paste0(plot.quant.dat$Quantiles, ", n = ", as.numeric(plot.quant.dat$n.plts)) # Labels for the quantile grid plot
  
  quant.grid.plt.fcn <- function(quants, quant.index) {
    ggplot(data = plot.quant.dat %>% filter(Quantiles %in% quants), aes(Quantiles, Means, fill = Quantiles)) + 
      geom_col() + 
      geom_errorbar(aes(ymax = UCI.95, ymin = LCI.95), width = 0.1) + 
      scale_fill_manual(values = virid.use[quant.index]) + 
      scale_y_continuous(limits = c(0, qt.max)) +
      scale_x_discrete(labels = q.g.p.labs[quant.index]) +
      theme_bw() +
      theme(legend.position = "none") + 
      labs(x = NULL, y = "Mean") +
      theme(text = element_text(size = 7))
  }
  
  p1 <- quant.grid.plt.fcn(c("LiHc", "MiHc", "HiHc"), 1:3) 
  p2 <- quant.grid.plt.fcn(c("LiMc", "MiMc", "HiMc"), 4:6) 
  p3 <- quant.grid.plt.fcn(c("LiLc", "MiLc", "HiLc"), 7:9) 
  
  p_all <- plot_grid(p1, p2, p3, ncol = 1)
  
  quant.table1 <- quant.table %>% filter(Species == plot.spp)
  
  
  ## - Obtain figure of quantile distribution relative to climate info - ##
  q_plot_spp <- quant.matrix %>% dplyr::select(all_of(plot.spp), all_of(VAR1), all_of(VAR.DELT)) %>%
    filter(get(all_of(plot.spp)) > 0) 
  names(q_plot_spp)[1] <- "Quantile"
  
    n_plots <- get(all_of(plot.spp), quant.n)#  table(get(all_of(spp.id1), q_plot_spp))
    n_plots2 <- tibble(loc = 1:n_quant, n = n_plots) %>%
    left_join(tibble(Quantiles = quant.table1$Quantiles, loc = 1:n_quant), by = "loc")
  
  quant.lims.plt <- get(plot.spp, quant.lims)
  quant.lims.delt.plt <- get(plot.spp, quant.lims.delta)
  
  sppnum <- as.numeric(gsub("X", "", sppnum.to.plot))
  
  # Common and Genus/species name for plot title
  com.name <- spp.names$COMMON_NAME[spp.names$SPCD == sppnum]
  g.s.name <- paste(spp.names$GENUS[spp.names$SPCD == sppnum], spp.names$SPECIES[spp.names$SPCD == sppnum])
  
  # Need to isolate the correct colors if some quadrants are without plots:
  scatter.virid.use <- virid.use[n_plots2$loc[n_plots2$n > 0]]
  
  
  plot.vals.plt <- ggplot(q_plot_spp, aes(get(VAR1), get(VAR.DELT), color = factor(Quantile))) + 
    geom_point() + 
    #stat_ellipse(type = "norm", level = 0.95, col = "orange", lwd = 2) +
    geom_hline(yintercept = quant.lims.delt.plt[1], col = "blue") +
    geom_hline(yintercept = quant.lims.delt.plt[2], col = "blue") +
    geom_vline(xintercept = quant.lims.plt[1], col = "green") +
    geom_vline(xintercept = quant.lims.plt[2], col = "green") +
    theme_bw() + 
    scale_color_manual(values = scatter.virid.use, name = "Quantiles", labels = QUANT.LEVELS) + 
    labs(title = paste0("Plot conditions for ", com.name, ", ", g.s.name, ", ",  sppnum.to.plot),
         x = var.label, y = var.delt.label) +
    theme(text = element_text(size = 7))
  #  scale_color_viridis_d(name = "Quartiles", option = "H", begin = 0.1, end = 0.9) #+
  #labs(x = var1.lab, y = var.delt.lab, title = spp.name)
  
  
  comb.plt <- p_all + plot.vals.plt + plot_layout(widths = c(1, 1.3))
  
  #ggsave(paste0(RESULTS.LOC, "Mort_figs_", var.filename, "/",sppnum.to.plot, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")
  ggsave(paste0(RESULTS.LOC, switch(ANALYSIS.TYPE, "mort" = "Mort", "grow" = "Growth"), "_figs_", var.filename, "/", sppnum.to.plot, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")
  
}




