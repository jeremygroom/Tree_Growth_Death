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
quant.matrix <- dieout.dat$quant.matrix
centroid_array <- dieout.dat$centroid.array
centroid.array <- dieout.dat$centroid.array
quant.n <- dieout.dat$quant.n
quant.lims <- dieout.dat$quant.lims




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

q_means <- q_SE <- rep(NA, n_quant)


#for (q in 1:n_quant) {

q_dieout.fcn <- function(q1, ado.dat, all.dat) {
  # Multiplying all of the species dead values (or species total values) by the respective species masks for "LiLc".
  #   With this step complete, the values for that quantile can be calculated for all species.
  ado.vals.use <- ado.dat[6:ncol(ado.dat)] * quant.array[, q1, ] 
  
  dat.ado.use <- bind_cols(orig[, 1:6], ado.vals.use) #%>% 
  #    mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
  #    left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
  
  all.vals.use <- all.dat[6:ncol(all.dat)] * quant.array[, q1, ] 
  
  dat.all.use <- bind_cols(orig[, 1:6], all.vals.use) # %>% 
  #    mutate(State_Plot2 = as.numeric(paste0(PLOT_FIADB, STATECD))) %>%
  #    left_join(clim.dat, by = c("State_Plot2" = "State_Plot"))
  
  sel.spp <- 1
  # This function will find, for any species, the estimate for quantile q1.
  quant.est.spp <- function(spp.sel) {
    
    spp.id1 <- paste0("X", spp.sel)
    n_q <- get(spp.id1, quant.n)[q1]   # Are there fewer than N.PLOT.LIM? If so, we can skip the calcs.
    
    
    # Means plus other metrics that are carried over into the calculation of the variance		
    if (n_q < N.PLOT.LIM) { # If too few points, just enter NA
      q_means$R <- NA; q_means$Zt <- NA; q_means$Yv <- NA; 
    } else {
      q_means <- mean.q.fcn(dat_tot = dat.all.use, resp = dat.ado.use, spp.sel1 = spp.sel)
    }
    # Standard errors
    if (n_q < N.PLOT.LIM) {
      q_SE <- NA
    } else {
      q_SE <- sqrt(ratio.SE.fcn(d.all_z = dat.all.use, d.ado_y = dat.ado.use, meandat = q_means))
    }
    
    list(q_means.R = q_means$R, q_means.Zt = q_means$Zt, q_means.Yv = q_means$Yv,
         q_SE = q_SE, quant.n = n_q)
    
  }
  
  qX_ests <- map(spp.list, quant.est.spp)
  qX_table <- bind_cols(Quantiles = q1, spp_id = spp.list, bind_rows(qX_ests))
  return(qX_table)
}

# Running for all species, all quantiles: 30 seconds
all.quants <- map(1:n_quant, q_dieout.fcn, ado.dat = ado_vals, all.dat = all_vals)
quant.table <- bind_rows(all.quants) %>% arrange(spp_id, Quantiles)

## Function for plotting the mortality numbers by quantile and the distribution of plots in quantiles
pair.plts.fcn <- function(sppnum.to.plot){

plot.spp <- sppnum.to.plot
plot.quant.dat <- quant.table %>% 
  filter(spp_id == plot.spp) %>%
  #bind_cols(quant.n %>% dplyr::select(plot.spp)) %>%
  #rename("n" = plot.spp) %>%
  mutate(q_uci = q_means.R + 1.96 * q_SE,
         q_lci = q_means.R - 1.96 * q_SE,
         Quantiles = case_match(Quantiles, 1:9 ~ QUANT.LEVELS),
         Quantiles = factor(Quantiles, levels = QUANT.LEVELS)) %>%
  dplyr::select(-q_means.Yv, -q_means.Zt, -q_SE) 


# Details for plotting
virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_quant)  # Get colors for plotting

qt.max <- ceiling(max(plot.quant.dat$q_uci, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights


#ggplot(data = plot.quant.dat, aes(Quantiles, q_means.R, fill = Quantiles)) + 
#  geom_col() + 
#  geom_errorbar(aes(ymax = q_uci, ymin = q_lci), width = 0.25) + 
#  scale_fill_viridis_d(name = "Quartiles", option = "H", begin = 0.1, end = 0.9)


# Create a grid of plots to match bivariate plot

q.g.p.labs <- paste0(plot.quant.dat$Quantiles, ", n = ", plot.quant.dat$quant.n) # Labels for the quantile grid plot

quant.grid.plt.fcn <- function(quants, quant.index) {
  ggplot(data = plot.quant.dat %>% filter(Quantiles %in% quants), aes(Quantiles, q_means.R, fill = Quantiles)) + 
    geom_col() + 
    geom_errorbar(aes(ymax = q_uci, ymin = q_lci), width = 0.1) + 
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

quant.table1 <- quant.table %>% filter(spp_id == spp.sel)


# Obtain figure of quantile distribution relative to climate info


spp.id1 <- spp.id[which(spp.list == plot.spp)]
q_plot_spp <- quant.matrix %>% dplyr::select(spp.id1, pre.vpdmin, delt.vpdmin) %>%
  filter(get(all_of(spp.id1)) > 0) 
names(q_plot_spp)[1] <- "Quantile"
#q_plot_spp %>% mutate(Quantile = case_match(Quantile,  (1:n_quant) ~ QUANT.LEVELS))



n_plots <- get(all_of(spp.id1), quant.n)#  table(get(all_of(spp.id1), q_plot_spp))
n_plots2 <- tibble(loc = 1:n_quant, n = n_plots) %>%
  left_join(tibble(Quantiles = quant.table1$Quantiles, loc = 1:n_quant), by = "loc")

dim(quant.array)

quant.lims <- dieout.dat$quant.lims

quant.lims.plt <- get(spp.id1, quant.lims)

# Common and Genus/species name for plot title
com.name <- spp.names$COMMON_NAME[spp.names$SPCD == sppnum.to.plot]
g.s.name <- paste(spp.names$GENUS[spp.names$SPCD == sppnum.to.plot], spp.names$SPECIES[spp.names$SPCD == sppnum.to.plot])

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

ggsave(paste0(RESULTS.LOC, "Test_figs_vpdmin/",spp.id1, "_plots.png"), comb.plt, device = "png", width = 7, height = 4, units = "in")

}

map(spp.list, pair.plts.fcn)

