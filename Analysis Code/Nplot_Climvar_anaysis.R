## Species distributions, plot numbers, and climate variable distributions ##

# --> Run Data_Prep_Analysis.R up through data.use (around line 91)

### -- Number of plots per species -- ###

n.spp.plts <- data.use %>% select(SPCD, puid) %>%
  distinct() %>%
  group_by(SPCD) %>%
  reframe(n = n())

hist(n.spp.plts$n, breaks = 15)
data.frame(n.spp.plts %>% arrange(n))
ggplot(n.spp.plts %>% arrange(n), aes(x = 1:nrow(n.spp.plts), y = n)) +
  geom_point() + 
  geom_line() + 
  labs(x = "Negative species rank by plot number", y = "Number of plots", 
       title = "Number of plots per species. Horizontal line at 400 plots") + 
  theme_bw() + 
  geom_hline(yintercept = 400)

spp.400 <- n.spp.plts %>% filter(n > 400) %>% 
  left_join(spp.names %>% select(-VARIETY, -SUBSPECIES, -WEST, -EXISTS_IN_PNWRS, -STOCKING_SPGRPCD, -FOREST_TYPE_SPGRPCD), by = "SPCD")

write_csv(spp.400, "spp.400.csv")


# Elevation data
elev.dat <- tree.plt.data %>% select(puid, ELEV) %>% distinct() %>% tibble()


# Exploring CWD
ggplot(climate.use, aes(difference)) + 
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0)

## Map of CWD

# Histogram of the 'diff' in CWD:
map.dat.1 <- climate.use %>%
  filter(is.na(difference) == FALSE,
         difference >= -3) %>% 
  ungroup() %>%
  left_join(elev.dat, by = 'puid')

# Quantiles of data below different thresholds
n.cwd <- nrow(climate.use)
n.zero <- nrow(climate.use %>% filter(difference < 0))
n.min.3 <- nrow(climate.use %>% filter(difference < -3))

cwd.lt.0.pct <- n.zero/n.cwd #Percent CWD diff < 0  =  8%
cwd.lt.min3.pct <- n.min.3/n.cwd #Percent CWD diff < -3  = 0.6%



# Points defining the boundaries of the map
maxlat <- max(map.dat.1$LAT); minlat <- min(map.dat.1$LAT)
maxlong <- max(map.dat.1$LON); minlong <- min(map.dat.1$LON)

# Refining the color palette based on the number of quadrants occupied.
#  The scale is already truncated from 0.1 to 0.9, so if 3 is the minimum, 0.3 works well.
#virid.min <- min(n.plots.used$loc[n.plots.used$n > 0]) / 10
#virid.max <- max(n.plots.used$loc[n.plots.used$n > 0]) / 10
#map.virid.use <- virid.use[n.plots.used$loc[n.plots.used$n > 0]]


cwd.map.fcn <- function(var.use, title.txt, lab.use){
ggplot(data = map.dat.1, aes(x = LON, y = LAT)) +
  coord_fixed(xlim = c(minlong - 1, maxlong + 1),  ylim = c(minlat - 1, maxlat + 1), ratio = 1.3) +
  geom_point(data = map.dat.1, aes(LON, LAT, color = get(var.use)), size = 0.5) +
  #  geom_tile(data = map.dat.1, mapping = aes(x = LON, y = LAT, z = targ.spp), binwidth = 0.15, 
  #            stat = "summary_2d", fun = mean, na.rm = TRUE, show.legend = FALSE) + 
  scale_color_viridis(option = "H", begin = 0.1, end = 0.9, name = lab.use) +  
  #                  labels = gsub(".", " ", QUANT.LEVELS[quants.used], fixed = TRUE)) +
  #scale_fill_viridis(option = "H", begin = virid.min, end = virid.max) +
  geom_polygon(data = west_df, mapping = aes(x = long, y = lat, group = group), color = "black", fill = "transparent") +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5)) + 
  labs(title = title.txt)
}

cwd.plt.startcwd <- cwd.map.fcn('pre_mean', title.txt = "FIA plot initial CWD values (10-yr mean prior to first visit)", lab.use = "Initial Climatic\nWater Deficit")
cwd.plt.deltcwd <- cwd.map.fcn('difference', title.txt = "FIA plot 10-year difference in CWD values", lab.use = "Climatic Water\nDeficit Change")

cwd.maps <- plot_grid(cwd.plt.startcwd, cwd.plt.deltcwd, nrow = 1)
ggsave("CWD_map.png", cwd.maps, device = "png")


#Elevation and change in Climatic Water Deficit
ggplot(map.dat.1, aes(ELEV, difference)) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()




### --- Number of plots per INVYR per spp
# Are plot-years for species balanced? That is, are the same number of plots with a species on them visited
#  every year?  
spp.consid <- as.numeric(c("11", "15", "17", "19", "64", "73", "81", "93", "108", "116", "117", "119", "122", "202", "242", "263", "264", "312", "351", "361", "631", "805", "818"))

plts.yr.fcn <- function(spp.x){
n.use <- paste0("n", spp.x) 
data.use %>% filter(SPCD == spp.x) %>% 
  select(INVYR, puid) %>% 
  distinct() %>%
  group_by(INVYR) %>%
  reframe(n = n()) %>%
  rename(!!n.use := n)
}

plts.yr <- map(spp.consid, plts.yr.fcn)
plts.yr.all <- reduce(plts.yr, left_join, by = "INVYR")
# Yes, aside from 2011 (which is a spill-over year) there is balance.  Roughly the same number of plots are sampled per year.  





# Investigating CWD over time
# Starting with summer_means from Climate_Dat_Preparation


summer_means %>% filter(is.na(mean) == FALSE)

ggplot(summer_means %>% filter(puid == "10094_41_57"), aes(year, mean)) + geom_line()

#ggplot(summer_means, aes(year, mean)) + geom_hex() + theme_bw() + geom_smooth(method = "lm")

# Trying all 47k plots with geom_path and alpha transparency
tic()
g_alpha <- ggplot(summer_means, aes(year, mean, group = puid)) + geom_path(alpha = 0.01)
ggsave("g_alpha.png", g_alpha, dev = "png")
toc()


# Let's sample 1k plots...
samp_puid <- sample(unique(summer_means$puid), 1000, replace = FALSE) 
s2 <- summer_means %>% filter(puid %in% samp_puid)

g_alpha2 <- ggplot(s2, aes(year, mean, group = puid)) + geom_path(alpha = 0.03) + theme_bw()
ggsave("g_alpha_1k.png", g_alpha2, dev = "png")


# I'm seeing that during some years CWD goes down sharply for a lot of plots, but goes up for others.  
# I could.... plot by year, plot by starting CWD (0 to 100, 100 to 175, 175 and above), etc.  
# I am tempted to see which plots geographically are in each category.  


#I think we may want to rerun this material for all plots with elevation and lat/lon

puid_use <- unique(tree.plt.data$puid)

s3 <- summer_means %>% filter(puid %in% puid_use) %>%
  left_join(latlon_puid, by = "puid") %>%
  group_by(puid) %>% 
  mutate(delta = mean - lag(mean)) %>%
  ungroup()

samp_puid2 <-  sample(unique(s3$puid), 1000, replace = FALSE) 
s4 <- s3 %>% filter(puid %in% samp_puid2)

g_alpha_s3 <- ggplot(s3, aes(year, mean, group = puid)) + geom_path(alpha = 0.01)
ggsave("g_alpha_s3.png", g_alpha_s3, dev = "png")

ggplot(s4, aes(year, mean, group = puid)) + geom_path(alpha = 0.03) + theme_bw()



# How does the change in CWD vary by year?

ggplot(s3 %>% filter(year > 1990), aes(year, delta, group = year)) + geom_violin()



g_line_violin <- ggplot(s3, aes(year, mean, group = puid)) + 
  geom_path(alpha = 0.01) + 
  geom_violin(data = s3 %>% filter(year > 1990), aes(year, delta, group = year), fill = "blue")

ggsave("g_s3_line_violin.png", g_line_violin, dev = "png")


g_line_delta <- ggplot(s4, aes(year, delta, group = puid)) + geom_path(alpha = 0.01) + theme_bw()
ggsave("g_delta_s4.png", g_line_delta, dev = "png")




