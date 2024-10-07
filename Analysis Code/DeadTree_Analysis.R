source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))

# Obtaining the data to work with: 
# Opening the zipped RDS file, output from Data_Prep.R
dieout.dat <- readr::read_rds(unzip(paste0(DATA.LOC, "dieout_dat.zip"), paste0(DATA.LOC,"dieout_dat.rds")))

# Deleting written RDS file (~150 MB)
if (file.exists(paste0(DATA.LOC, "dieout_dat.rds")) == TRUE) {
  file.remove(paste0(DATA.LOC, "dieout_dat.rds"))
}

spp.id <- dieout.dat$spp.id
spp.list <- dieout.dat$spp.list
ado_vals <- dieout.dat$ado_vals
all_vals <- dieout.dat$all_vals
quant.array <- dieout.dat$quant.array
state.array <- dieout.dat$state.array
state.n <- dieout.dat$state.n
quant.matrix <- dieout.dat$quant.matrix
centroid_array <- dieout.dat$centroid.array
centroid.array <- dieout.dat$centroid.array
quant.n <- dieout.dat$quant.n
quant.lims <- dieout.dat$quant.lims

## First, adjusting the species list
spp.list <- sort(as.numeric(gsub("X", "", SEL.SPP)))
spp.id <- paste0("X", spp.list)

# Questions about the plot data
# 1) Are there intensification plots
#sum(orig$intensification)  # Yes, 6327
# 2) How many forested plots?
#length(orig$propfor[orig$propfor > 0])  # 18672 
# 3) Do our trees exclusively occur in forested plots?
#orig2 <- orig %>% mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
#  left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
#orig2$trees <- apply(orig2[, 10:59], 1, sum)
#summary(orig2$propfor[orig2$trees > 0])
#orig2[orig2$propfor == 0 & orig2$trees > 0,] # 79 plots
# It is evident that the proportion of forest/non-forest/missed areas of plots is not necessarily constant between the two visits. 
# 4) We have climate data for 37,612 plots and plot data for 38656 plots.  Are we missing climate data for any plots with our species of interest?
#dim(clim.dat) ; dim(orig)   # 37612 13 , 38656 59
#orig2 <- orig %>% mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
#left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
#length(orig2$pre.temp[is.na(orig2$pre.temp) == FALSE])  # 37612; all climate variable plot values had matches
#x1 <- orig2[is.na(orig2$pre.temp) == TRUE,  ]
#sum(apply(x1[, 10:59], 1, sum))   # result = zero.  None of these plots contain our trees of interest.



VAR1 <- "pre.vpdmin"
VAR.DELT <- "delt.vpdmin"

# temporary: will change as this is automated # 
#quant.sel <- "LiLc"                            # Select the quantile
#spp.sel <- as.numeric(gsub("X", "", SEL.SPP))  # Select the species (set above)
#spp.sel.index <- which(colnames(orig)[grep("X", colnames(orig))] == SEL.SPP)
n.spp <- length(spp.id)



#quant.vals <- which(QUANT.LEVELS == quant.sel) 

n_quant <- length(QUANT.LEVELS)

q_means <- q_SE <- q_bs.UCI <- q_bs.LCI <- rep(NA, n_quant)

# Strata for bootstrap resampling procedure
strata <- unique(ado_vals$STRATUM)
strata.num <- data.frame(val = 1:nrow(ado_vals), stratum = ado_vals$STRATUM) 


q_dieout.fcn <- function(q1, ado.dat, all.dat, array.name, category.n) {
  # Multiplying all of the species dead values (or species total values) by the respective species masks for "LiLc".
  #   With this step complete, the values for that quantile can be calculated for all species.
  ado.vals.use <- ado.dat[6:ncol(ado.dat)] * array.name[, q1, ] 
  
  dat.ado.use <- bind_cols(orig[, 1:6], ado.vals.use) #%>% 
  #    mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
  #    left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
  
  all.vals.use <- all.dat[6:ncol(all.dat)] * array.name[, q1, ] 
  
  dat.all.use <- bind_cols(orig[, 1:6], all.vals.use) # %>% 
  #    mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
  #    left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))

  bs.out <- data.frame(matrix(NA, nrow = BS.N, ncol = length(spp.list)))
  # Bootstrap function for obtaining an estimate of the mean

  bs.fcn <- function(iter, d.all_use, d.sub_use) { # Iteration number, all of the relevant data, the subset of the relevant data

    # First, sample the tables in use. This procedure obtains samples with replacement from each stratum and then combines the selected
    #  row numbers.  
    samp <- unlist(map(strata, function(x) sample(strata.num$val[strata.num$stratum == x], replace = TRUE)))
    

    temp.spp <- rep(NA, length(spp.id))
    
    # This function will find, for any species, the estimate for quantile q1.
    quant.est.spp <- function(spp.sel) {
      
      spp.id1 <- paste0("X", spp.sel)
      n_q <- get(spp.id1, category.n)[q1]   # Are there fewer than N.PLOT.LIM? If so, we can skip the calcs.
      
      # Means plus other metrics that are carried over into the calculation of the variance		
      if (n_q < N.PLOT.LIM) { # If too few points, just enter NA
        q_bs.mean <- NA 
      } else {
        q_bs.mean <- mean.q.fcn(dat_tot = d.all_use[samp, ], resp = d.sub_use[samp, ], spp.sel1 = spp.sel)$R # Mean given the selected rows
      }
      return(q_bs.mean)
    }
    bs.out[iter, ] <- unlist(map(spp.list, quant.est.spp))
   }
      # Running the parallel portion of the function
  furrr.out <- future_map(1:BS.N, bs.fcn, d.all_use = dat.all.use, d.sub_use = dat.ado.use, .options = furrr_options(seed = TRUE))  
  
  # Processing the results and returning quantiles/mean
  furrr.table <- matrix(unlist(furrr.out), ncol = length(spp.list), byrow = TRUE)
  quants <- data.frame(t(apply(furrr.table, 2, function(x) quantile(x, probs = c(0.5, 0.025, 0.975), na.rm = TRUE))))
  names(quants) <- c("Median", "LCI.95", "UCI.95")
  quants$Means <- apply(furrr.table, 2, mean, na.rm = TRUE)
  quants$n.plts <- as.vector(category.n[q1, 1:length(spp.id) ])
  quants$Species <- spp.id
  quants$Quantile <- q1
  
  return(quants)
}


plan(multisession, workers = n.cores) # Setting up parallel computing


# Running for all species, all quantiles: 30 seconds
y <- Sys.time()  # 16 minutes on new computer at 1000 iterations (9.4 min for 9 spp)
all.quants <- map(1:n_quant, q_dieout.fcn, ado.dat = ado_vals, all.dat = all_vals, array.name = quant.array, category.n = quant.n)
quant.table <- bind_rows(all.quants) %>% arrange(Species, Quantile)
Sys.time() - y


# Running for states.  11.5 min on new machine, 16 min on old machine  (2.4 minutes using new machine, 9 species)
state.list <- unique(ado_vals$STATECD) # 6 = California, 41 = Oregon, 53 = Washington
y <- Sys.time()  # 16 minutes on new computer at 1000 iterations
state.ests <- map(1:length(state.list), q_dieout.fcn, ado.dat = ado_vals, all.dat = all_vals, array.name = state.array, category.n = state.n)
Sys.time() - y
#state.table <- bind_rows(state.ests)
#write_rds(state.ests, file = paste0(RESULTS.LOC, "Estimates/State_Ests_2024.rds"))


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
q.g.p.labs <- paste0(plot.quant.dat$Quantiles, ", n = ", plot.quant.dat$quant.n) # Labels for the quantile grid plot

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


# Obtain figure of quantile distribution relative to climate info


#spp.id1 <- spp.id[which(spp.list == plot.spp)]
q_plot_spp <- quant.matrix %>% dplyr::select(all_of(plot.spp), pre.vpdmin, delt.vpdmin) %>%
  filter(get(all_of(plot.spp)) > 0) 
names(q_plot_spp)[1] <- "Quantile"
#q_plot_spp %>% mutate(Quantile = case_match(Quantile,  (1:n_quant) ~ QUANT.LEVELS))



n_plots <- get(all_of(plot.spp), quant.n)#  table(get(all_of(spp.id1), q_plot_spp))
n_plots2 <- tibble(loc = 1:n_quant, n = n_plots) %>%
  left_join(tibble(Quantiles = quant.table1$Quantiles, loc = 1:n_quant), by = "loc")

#dim(quant.array)

quant.lims <- dieout.dat$quant.lims

quant.lims.plt <- get(plot.spp, quant.lims)

sppnum <- as.numeric(gsub("X", "", sppnum.to.plot))

# Common and Genus/species name for plot title
com.name <- spp.names$COMMON_NAME[spp.names$SPCD == sppnum]
g.s.name <- paste(spp.names$GENUS[spp.names$SPCD == sppnum], spp.names$SPECIES[spp.names$SPCD == sppnum])

plot.vals.plt <- ggplot(q_plot_spp, aes(pre.vpdmin, delt.vpdmin, color = factor(Quantile))) + 
  geom_point() + 
  #stat_ellipse(type = "norm", level = 0.95, col = "orange", lwd = 2) +
  geom_hline(yintercept = quant.lims.plt[1, 2], col = "blue") +
  geom_hline(yintercept = quant.lims.plt[2, 2], col = "blue") +
  geom_vline(xintercept = quant.lims.plt[1, 1], col = "green") +
  geom_vline(xintercept = quant.lims.plt[2, 1], col = "green") +
  theme_bw() + 
  scale_color_manual(values = virid.use, name = "Quantiles", labels = QUANT.LEVELS) + 
  labs(title = paste0("Plot conditions for ", com.name, ", ", g.s.name, ", ",  spp.id1)) +
  theme(text = element_text(size = 7))
#  scale_color_viridis_d(name = "Quartiles", option = "H", begin = 0.1, end = 0.9) #+
#labs(x = var1.lab, y = var.delt.lab, title = spp.name)


comb.plt <- p_all + plot.vals.plt + plot_layout(widths = c(1, 1.3))

ggsave(paste0(RESULTS.LOC, "Test_figs_vpdmin/",sppnum.to.plot, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")

}

map(spp.id, pair.plts.fcn)

