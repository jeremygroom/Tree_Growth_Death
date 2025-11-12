## This code provides information summaries for the manuscript and generates 
#   non-estimate figures and (maybe) tables.
#   The file Data_Prep_Analysis.R sources this code.  
#   The file Manuscript.Rmd will open and use its outputs.


### Q: How many trees did we track? How many plots?
tree.num <- dim(tree.plt.data %>% filter(SPCD %in% sel.spp))[1]
plot.num <- dim(PlotDat)[1]
strat.num <- length(unique(PlotDat$stratum))

## Q: What percentage of plots with 9, 10, and 11 yr remeasure?
remeas <- tree.plt.data %>% dplyr::select(puid, REMPER) %>% distinct() %>% mutate(REMPER2 = round(REMPER, digits = 0))
pct.remeas <- round(100 * table(remeas$REMPER2)/nrow(remeas), digits = 1)[4:6]


## Number of Gymnosperms?  Which ones?
gym.angio.dat <- spp.names %>% filter(SPCD %in% sel.spp) %>%
  mutate(sci_name = paste(GENUS, SPECIES),
         SpeciesCode = paste0("X", SPCD)) %>%
  dplyr::select(SpeciesCode, SPCD, COMMON_NAME, SPECIES_SYMBOL, SFTWD_HRDWD, sci_name)
n.gym.spp <- length(gym.angio.dat$SPCD[gym.angio.dat$SFTWD_HRDWD == "S"])


### Q: What percentage of change in CWD is negative & positive?
# Quantiles of data below different thresholds
n.cwd <- nrow(climate.use)
n.zero <- nrow(climate.use %>% filter(difference < 0))
n.min.3 <- nrow(climate.use %>% filter(difference < -3))

cwd.lt.0.pct <- round(n.zero/n.cwd, 3) * 100 #Percent CWD diff < 0  =  8%
cwd.gt.0.pct <- round(1 - n.zero/n.cwd, 3) * 100 # Percent CWD diff > 0 (91.6%)
cwd.lt.min3.pct <- n.min.3/n.cwd #Percent CWD diff < -3  = 0.6% (This is useful for the mapping)



mort.out <- readRDS(paste0(save.loc.fcn(2), "Domain_Analysis_Output.RDS"))

init.cwd.stats <- summ.spp.fcn(data.all = mort.out$all_dat, "pre_mean", gym.angio.dat) 
diff.cwd.stats <- summ.spp.fcn(data.all = mort.out$all_dat, "difference", gym.angio.dat)

init.cwd.min.mean <- which(init.cwd.stats$Mean == min(init.cwd.stats$Mean))
init.cwd.max.mean <- which(init.cwd.stats$Mean == max(init.cwd.stats$Mean))
init.range = init.cwd.stats$q95 - init.cwd.stats$q05
init.cwd.90range.min <- which(init.range == min(init.range))
init.cwd.90range.max <- which(init.range == max(init.range))


diff.cwd.min.mean <- which(diff.cwd.stats$Mean == min(diff.cwd.stats$Mean))
diff.cwd.max.mean <- which(diff.cwd.stats$Mean == max(diff.cwd.stats$Mean))
diff.range = as.vector(diff.cwd.stats$q95 - diff.cwd.stats$q05)
diff.cwd.90range.min <- diff.range[which(diff.range == min(diff.range))]
diff.cwd.90range.max <- diff.range[which(diff.range == max(diff.range))]
diff.cwd.95max <- which(diff.cwd.stats$q95 == max(diff.cwd.stats$q95))

# For supplemental, MC permutation p-values
psig.dat <- read_csv(paste0(RESULTS.OTHER, "Permutation_results.csv"), show_col_types = FALSE)
psig.dat2 <- psig.dat %>% dplyr::select(-order.c) %>%
  left_join(gym.angio.dat %>% dplyr::select(SpeciesCode, sci_name, SFTWD_HRDWD), by = c("Species" = "SpeciesCode")) %>%
  arrange(desc(SFTWD_HRDWD), sci_name) %>%
  dplyr::select(-SFTWD_HRDWD, -Species) %>%
  relocate(sci_name)

num.cols <- grep("sig.dist", names(psig.dat2) )

psig.dat2[, num.cols] <- apply(psig.dat2[, num.cols], 2, as.character)
psig.dat2[is.na(psig.dat2)] <- ""

analysis.stats <- list(tree.num = tree.num,
                       plot.num = plot.num,
                       strat.num = strat.num,
                       pct.remeas = pct.remeas,
                       cwd.lt0.pct = cwd.lt.0.pct,
                       cwd.gt0.pct = cwd.gt.0.pct,
                       gym.angio.dat = gym.angio.dat,
                       n.gym.spp = n.gym.spp,
                       init.cwd.stats = init.cwd.stats,
                       diff.cwd.stats = diff.cwd.stats,
                       psig.dat2 = psig.dat2,
                       init.cwd.min.mean = init.cwd.min.mean, 
                       init.cwd.max.mean = init.cwd.max.mean, 
                       init.cwd.90range.min = init.cwd.90range.min,
                       init.cwd.90range.max = init.cwd.90range.max,
                       diff.cwd.min.mean = diff.cwd.min.mean, 
                       diff.cwd.max.mean = diff.cwd.max.mean, 
                       diff.cwd.90range.min = diff.cwd.90range.min,
                       diff.cwd.90range.max = diff.cwd.90range.max,
                       diff.cwd.95max = diff.cwd.95max
) 

write_rds(analysis.stats, paste0(RESULTS.OTHER, "analysis.stats.RDS"))




##### Figure: Example growth, mortality and CWD scatterplot -----------------------
domain.figs <- list()
for(k in 1:2){
  plt.dat <- readRDS(paste0(save.loc.fcn(k), "Domain_Analysis_Output.RDS"))
  plt.dat2 <- plt.dat$domain.summaries
  
  # Need climate names for files and axes.
  clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE) %>%
    filter(filename == CLIM.VAR.USE)
  
  var1 <- clim.names$values[grep("pre", clim.names$values)]
  var.delt <- clim.names$values[grep("d", clim.names$values)]
  var.label <- clim.names$label[clim.names$values == var1]
  #var.filename <- clim.names$filename[clim.names$values == var1]
  var.delt.label <- clim.names$label[clim.names$values == var.delt]
  
  
  sppnum.to.plot <- "X81"
  mort.grow.dat <- get(paste0("arrays.", ANALYSIS.TYPE[k])) 
  
  domain.matrix <- mort.grow.dat$domain.matrix
  domain.n <- mort.grow.dat$domain.n 
  quant.lims <- mort.grow.dat$quant.lims
  quant.lims.delta <- mort.grow.dat$quant.lims.delt
  
  
  use.dat2 <- plt.dat2 %>% filter(Species == sppnum.to.plot)
  use.dat2 <- cm2.fcn(k, use.dat2)  # Transforming growth from inches2 to cm2
  
  plot.domain.dat <- plt.dat2 %>% 
    filter(Species == sppnum.to.plot) %>%
    left_join(domain.level.table, by = c("Domain" = "q.num"))
  
  # Details for plotting
  virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_domain)  # Get colors for plotting
  
  qt.max <- ceiling(max(use.dat2$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
  
  # Create a grid of plots to match bivariate plot
  q.g.p.labs <- paste0(DOMAIN.LEVELS, ", n = ", as.numeric(use.dat2$n.plts)) # Labels for the domain grid plot
  
  ylabs <- if (k == 1) "Mean Growth (cm\u00B2/decade)" else "Mean Decadal Mortality Rate"
  
  text.size <- 9
  
  p1 <- domain.grid.plt.fcn(c("DL", "DM", "DH"), 1:3, use.dat2, qt.max, q.g.p.labs, ylabs, text.size) 
  p2 <- domain.grid.plt.fcn(c("SL", "SM", "SH"), 4:6, use.dat2, qt.max, q.g.p.labs, ylabs, text.size) 
  #p_all <- plot_grid(p1, p2, ncol = 1)
  list1 <- list(p1 = p1, p2 = p2)
  fig.letter <- switch(k, "A", "B")
  domain.figs[[paste0("p1_", k)]] <- p1 + annotate("text", x = 0.65 , y = qt.max * 0.95, label = fig.letter, parse = TRUE, size = 6) + theme(text = element_text(size = text.size))
  domain.figs[[paste0("p2_", k)]] <- p2 + theme(text = element_text(size = text.size))
  
  
  # Plotting the scatterplot of points in quadrants; only doing so for Mortality.
  if(k == 2) {
    # Code to prep for scatterplot figure - helps in reducing the color palette.
    n_plots <- get(sppnum.to.plot, domain.n)
    n_plots2 <- tibble(loc = 1:n_domain, n = n_plots) %>%
      left_join(tibble(Domain = plot.domain.dat$Domains, loc = 1:n_domain), by = "loc")
    
    
    plot.vals.plt <- domain.dist.plt.fcn(plot.spp = sppnum.to.plot, domain.matrix = domain.matrix, 
                                         var1 = var1, var.delt = var.delt, quant.lims = quant.lims,
                                         n.plots.used = n_plots2, size.trees = "", text.size) 
    
    q_plot_spp <- domain.matrix %>% dplyr::select(all_of(sppnum.to.plot), all_of(var1), all_of(var.delt)) %>%
      filter(get(sppnum.to.plot) > 0) 
    pts.max.y <- ceiling(layer_scales(plot.vals.plt)$y$range$range[2] * 10) / 10 # Extracts max Y extent of ggplot
    pts.min.y <- ceiling(layer_scales(plot.vals.plt)$y$range$range[1] * 10) / 10 # Extracts max Y extent of ggplot
    
    #ceiling(max(use.dat2$UCI.95, na.rm = TRUE) * 10) / 10 # For consistent y-axis heights
    
    
    range.x <- range(q_plot_spp$pre_mean, na.rm = TRUE)
    pts.min.x <- min(q_plot_spp$pre_mean) + 0.065 * (range.x[2] - range.x[1]) 
    pts.max.x <- max(q_plot_spp$pre_mean)
    
    # Locations for domain text in figure:
    q.lim.n <- which(names(quant.lims) == sppnum.to.plot)
    q.lim.plot <- quant.lims[[q.lim.n]]
    # Low/med/high positions
    l.txt.loc <- q.lim.plot[1] - (q.lim.plot[1] - pts.min.x)/2
    m.txt.loc <- q.lim.plot[1] + (q.lim.plot[2] - q.lim.plot[1])/2
    h.txt.loc <- q.lim.plot[2] + (pts.max.x - q.lim.plot[2])/2
    # Above/below positions
    a.txt.loc <- pts.max.y - (pts.max.y - 5)/2
    b.txt.loc <- pts.min.y + (5 - pts.min.y)/2
    
    txt.loc.x.vec <- as.vector(rep(c(l.txt.loc, m.txt.loc, h.txt.loc), 2))
    txt.loc.y.vec <- as.vector(c(rep(a.txt.loc, 3), rep(b.txt.loc, 3)))
    
    plot.vals.plt <- plot.vals.plt + 
      labs(title = NULL) + 
      annotate("text", x = pts.min.x , y = pts.max.y, label = "C", parse = TRUE, size = 6) +
      annotate("text", x = txt.loc.x.vec , y = txt.loc.y.vec, label = DOMAIN.LEVELS, parse = TRUE, size = 4) +
      theme(text = element_text(size = text.size)) 
    
    plot.vals.plt.no_legend <- plot.vals.plt + theme(legend.position = 'none')
  }
  
}

p1_1 <- domain.figs$p1_1
p1_2 <- domain.figs$p1_2
p2_1 <- domain.figs$p2_1
p2_2 <- domain.figs$p2_2 

comb.plt <- ((p1_1/p2_1 + plot_layout(axis_titles = "collect")) | (p1_2/p2_2 + plot_layout(axis_titles = "collect")) | plot.vals.plt.no_legend ) #/ guide_area() + plot_layout(guides = 'collect', heights = c(10, 0.01)) 

comb.plt.legend <- get_legend(plot.vals.plt)

comb.plt <- plot_grid(comb.plt, comb.plt.legend, ncol = 1, rel_heights = c(1, 0.03))


ggsave(paste0(RESULTS.OTHER, "Example_spp_figure.png"), comb.plt, device = 'png',
       width = 8, height = 6, dpi = 300, units = "in", bg = "white")






##### Figure: Map of plots displaying initial and delta CWD -----------------------
us_states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
west_states <- us_states[us_states$ID %in% c("california", "oregon", "washington"), ]

# Transform to a suitable projection for the west coast (e.g., NAD83 / California Albers)
west_states_proj <- st_transform(west_states, crs = 3310)

# Prepare climate data
map.dat.1 <- climate.use %>%
  filter(is.na(difference) == FALSE,
         difference >= -3) %>% 
  ungroup()

# Convert climate data to sf object and transform to same projection
climate_sf <- st_as_sf(map.dat.1, coords = c("LON", "LAT"), crs = 4326)
climate_sf_proj <- st_transform(climate_sf, crs = 3310)

# Extract coordinates in projected system
coords_proj <- st_coordinates(climate_sf_proj)
climate_proj_df <- climate_sf_proj %>%
  st_drop_geometry() %>%
  bind_cols(data.frame(X = coords_proj[,1], Y = coords_proj[,2]))

# Get bounding box for the map extent
bbox <- st_bbox(west_states_proj)
x_range <- c(bbox["xmin"] - 50000, bbox["xmax"] + 50000)
y_range <- c(bbox["ymin"] - 50000, bbox["ymax"] + 50000)

# Create inset map of US with west states highlighted
us_states_full <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
# Transform to Albers Equal Area for US
us_states_albers <- st_transform(us_states_full, crs = 5070)
west_states_albers <- us_states_albers[us_states_albers$ID %in% c("california", "oregon", "washington"), ]

inset_map <- ggplot() +
  geom_sf(data = us_states_albers, fill = "white", color = "gray60", size = 0.3) +
  geom_sf(data = west_states_albers, fill = "gray60", color = "black", size = 0.5) +
  theme_void() +
  theme(panel.background = element_rect(fill = "lightblue", color = "black", linewidth = 0.5),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5))


# Create a data frame with mountain range names and approximate coordinates
mountain_ranges <- data.frame(
  name = c("Cascade Range", "Sierra Nevada", "Coast Ranges", 
           "Blue Mountains", "Klamath\nMountains"),
  #lon = c(-121.5, -119.5, -123.5, -123.5, -118.5, -123.0),
  lon = c(-122, -119.5, -123.8,  -118.5, -123.0),
  lat = c(45.5, 37.5, 44.0, 45.2, 42)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 3310)  # Match your projection

# Extract coordinates for plotting
coords <- st_coordinates(mountain_ranges)
mountain_ranges_df <- cbind(mountain_ranges, coords) %>%
  st_drop_geometry() %>%
  mutate(ang = c(82, -55, 82, 65, 0),
         #colr = c("white", "black", "white", "white", "black", "black"))
         colr = rep("black", 5))



# Enhanced mapping function
cwd.map.fcn <- function(var.use, title.txt, lab.use, add_extras = FALSE) {
  p <- ggplot() +
    # Add state boundaries
    # Add climate data points
    geom_point(data = climate_proj_df, 
               aes(x = X, y = Y, color = get(var.use)), 
               size = 0.8, alpha = 0.8) +
    scale_color_viridis_c(option = "H", begin = 0.1, end = 0.9, 
                          name = lab.use,
                          guide = guide_colorbar(
                            title.position = "top",
                            title.hjust = 0.5,
                            barwidth = 8,
                            barheight = 0.8,
                            frame.colour = "black",
                            ticks.colour = "black"
                          )) +
    geom_sf(data = west_states_proj, fill = "transparent", color = "black", size = 0.8) +
    coord_sf(xlim = x_range, ylim = y_range, expand = FALSE) +
    geom_text(data = mountain_ranges_df,
              aes(x = X, y = Y, label = name, angle = ang),
              size = 3, fontface = "bold.italic", color = mountain_ranges_df$colr,
              check_overlap = TRUE) +
    theme_void() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.box.margin = margin(t = -5, b = 5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(title = title.txt)
  
  # Add extras to the first subplot
  if (add_extras) {
    # Add north arrow
    p <- p + annotation_north_arrow(
      location = "bl", 
      which_north = "true",
      pad_x = unit(0.4, "npc"),
      pad_y = unit(0.02, "npc"),
      style = north_arrow_fancy_orienteering(
        fill = c("white", "black"),
        line_col = "black",
        text_size = 8
      ),
      height = unit(1.2, "cm"),
      width = unit(1.2, "cm")
    )
    
    # Add scale bar
    p <- p + annotation_scale(
      location = "bl",
      width_hint = 0.25,
      pad_x = unit(0.02, "npc"),
      pad_y = unit(0.02, "npc"),
      style = "ticks",
      line_col = "black",
      text_col = "black",
      height = unit(0.3, "cm")
    )
    
    # Add inset map
    p <- p + annotation_custom(
      ggplotGrob(inset_map),
      xmin = x_range[1] + 0.5 * diff(x_range),
      xmax = x_range[1] + 0.97 * diff(x_range),
      ymin = y_range[1] + 0.37 * diff(y_range),
      ymax = y_range[1] + 0.55 * diff(y_range)
    )
  }
  
  return(p)
}

# Create the two subplots
cwd.plt.startcwd <- cwd.map.fcn(
  'pre_mean', 
  title.txt = NULL, 
  lab.use = "Initial Climatic Water Deficit (mm)",
  add_extras = TRUE
) + 
  annotate("text", x = 0.05 * (x_range[2] - x_range[1]) + x_range[1], 
           y = 0.97 * (y_range[2] - y_range[1]) + y_range[1],
           label = "bold(A)", parse = TRUE, size = 6) +
  annotate("text", x = 0.54 * (x_range[2] - x_range[1]) + x_range[1], 
           y = 0.415 * (y_range[2] - y_range[1]) + y_range[1],
           label = "bold(B)", parse = TRUE, size = 6)


cwd.plt.deltcwd <- cwd.map.fcn(
  'difference', 
  title.txt = NULL, 
  lab.use = "Climatic Water Deficit Change (mm)"
) + 
  annotate("text", x = 0.05 * (x_range[2] - x_range[1]) + x_range[1], 
           y = 0.97 * (y_range[2] - y_range[1]) + y_range[1],
           label = "bold(C)", parse = TRUE, size = 6)

# Combine the plots
cwd.maps <- plot_grid(
  cwd.plt.startcwd, 
  cwd.plt.deltcwd, 
  nrow = 1, 
  align = "hv"
) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "null"))


cwd.maps <- ggdraw() +
  draw_plot(cwd.plt.startcwd, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(cwd.plt.deltcwd, x = 0.5, y = 0, width = 0.5, height = 1) 


# Print the final figure
print(cwd.maps)

# Saving a png of the figure and a RDS version as well. The MS relies on the RDS.
ggsave(paste0(RESULTS.OTHER, "cwd_maps_figure.png"), cwd.maps, device = 'png',
       width = 6, height = 8.81, dpi = 300, units = "in",
       bg = "white")

ms.figures.other <- list(cwd.maps = cwd.maps)

write_rds(ms.figures.other, paste0(RESULTS.OTHER, "ms.figures.other.RDS"))




