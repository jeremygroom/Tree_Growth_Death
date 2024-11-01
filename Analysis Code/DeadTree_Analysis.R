source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))

# Set the climate varible to examine:

for (i in 1:length(CLIM.VAR)){
  for(j in 1:length(ANALYSIS.TYPE)) {
    
    var1 <- paste0("pre.", CLIM.VAR[i])
    var.delt <- paste0("delt.", CLIM.VAR[i])
    
    # Need climate names for files and axes.
    clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE)
    
    var.label <- clim.names$label[clim.names$values == var1]
    var.filename <- clim.names$filename[clim.names$values == var1]
    var.delt.label <- clim.names$label[clim.names$values == var.delt]
    
    # Obtaining the data to work with: 
    # Opening the zipped RDS file, output from Data_Prep.R
    mort.grow.dat <- readr::read_rds(unzip(paste0(DATA.LOC, ANALYSIS.TYPE[j], "_dat_", var.filename, ".zip"), paste0(DATA.LOC, ANALYSIS.TYPE[j], "_dat_", var.filename, ".rds")))
    
    # Deleting written RDS file (~150 MB)
    if (file.exists(paste0(DATA.LOC, ANALYSIS.TYPE[j], "_dat_" , var.filename, ".rds")) == TRUE) {
      file.remove(paste0(DATA.LOC, ANALYSIS.TYPE[j], "_dat_" , var.filename, ".rds"))
    }
    
    #spp.id <- mort.grow.dat$spp.id
    spp.list <- mort.grow.dat$spp.list
    vals_dat <- mort.grow.dat$vals_dat
    all_dat <- mort.grow.dat$all_dat
    quant.array <- mort.grow.dat$quant.array
    state.array <- mort.grow.dat$state.array
    state.n <- mort.grow.dat$state.n
    quant.matrix <- mort.grow.dat$quant.matrix
    centroid.array <- mort.grow.dat$centroid.array
    quant.n <- mort.grow.dat$quant.n 
    quant.lims <- mort.grow.dat$quant.lims
    quant.lims.delta <- mort.grow.dat$quant.lims.delt
    
    ## First, adjusting the species list
    #spp.id <- paste0("X", spp.list)     # Can use SEL.SPP
    spp.list <- as.numeric(gsub("X", "", SEL.SPP))
    
    # Need climate names for files and axes.
    clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE)
    
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
    
    
    
    n.spp <- length(SEL.SPP)
    
    n_quant <- length(QUANT.LEVELS)
    
    q_means <- q_SE <- q_bs.UCI <- q_bs.LCI <- rep(NA, n_quant)
    
    # Strata for bootstrap resampling procedure
    strata <- unique(vals_dat$STRATUM)
    strata.num <- data.frame(val = 1:nrow(vals_dat), stratum = vals_dat$STRATUM) 
    
    
    plan(multisession, workers = n.cores) # Setting up parallel computing
    
    
    ## Quantile estimates for mortality or growth, across select species.
    y <- Sys.time()  # 16 minutes on new computer at 1000 iterations (9.4 min for 9 spp)
    quant.index <- 1:n_quant
    all.quants <- quant.index %>% map(\(x) q_mort.grow.fcn(q1 = x, vals.dat = vals_dat, all.dat = all_dat, array.name = quant.array, category.n = quant.n) )  # Alternate mapping formulation, only advantage is clarity.
    quant.table <- bind_rows(all.quants) %>% arrange(Species, Quantile)
    Sys.time() - y
    
    
    # Running for states.  11.5 min on new machine, 16 min on old machine  (2.4 minutes using new machine, 9 species)
    state.list <- unique(vals_dat$STATECD) # 6 = California, 41 = Oregon, 53 = Washington
    y <- Sys.time()  # 16 minutes on new computer at 1000 iterations
    state.ests <- map(1:length(state.list), q_mort.grow.fcn, vals.dat = vals_dat, all.dat = all_dat, array.name = state.array, category.n = state.n)
    Sys.time() - y
    #state.table <- bind_rows(state.ests)
    #write_rds(state.ests, file = paste0(RESULTS.LOC, "Estimates/State_Ests_2024.rds"))
    
    
    # Plotting paired plots of mortality/growth by quantile and a scatterplot of plot distribution by quantiles.
    map(SEL.SPP, pair.plts.fcn, quant.table = quant.table, var.label, var.delt.label)
  }
}
