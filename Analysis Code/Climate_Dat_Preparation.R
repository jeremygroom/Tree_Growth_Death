### Actual Evapotranspiration data set manipulation for the analysis ###

clim_dat <- read_rds(file.path(CLIMATE.LOC, "extraction_070825.RDS")) # This is external of the GitHub repo

# The selected climate variable is CLIM.VAR.USE, set in Global.R
clim_var_name <- CLIM.VAR.USE
clim_var_listed <- if(clim_var_name == "cwd") "def" else clim_var_name

clim_dat <- clim_dat %>% filter(year > 1980, var == clim_var_listed) 

# Use the tree dataset from Data_Prep_Analysis.R to get the years needed for each plot
data_clim_yrs <- tree.plt.data %>% 
  select(puid, MEASYEAR, REMPER) %>%
  distinct() %>%
  mutate(YEAR1 = MEASYEAR - round(REMPER, digits = 0) + 1,
         PREYR = YEAR1 - 19) %>%
  select(puid, REMPER, PREYR, YEAR1, MEASYEAR)



#### Check on plot REMPER distribution------------------
#ggplot(data_clim_yrs, aes(REMPER)) + geom_histogram()
# plots:
#n.plt <- nrow(data_clim_yrs)
# Plots between 9 & 11 yrs
#n.9.11 <- nrow(data_clim_yrs %>% filter(REMPER >= 9 & REMPER <= 11))
#n.9.11/n.plt
# 93.1% of plots have REMPER between 9 and 11 yrs

#### Check on missing clim data -----------------------
#sub_clim <- clim_dat %>% dplyr::select(puid, var) %>% distinct()
#data_clim_check <- left_join(data_clim_yrs, sub_clim, by = "puid")  %>%
#  mutate(state = substr(puid, 7, 8),
#         state = as.numeric(gsub("_", "", state)))
#nrow(data_clim_check)      #12298 plots
#missing.clim <- data_clim_check %>% filter(is.na(var)) # 0 rows

#data_clim_check2 <- data_clim_check %>% tibble() %>% filter(is.na(var))



#### Analyzing climate variable  data  ------------------------------------------------

var1 <- clim_dat
data_yrs <- data_clim_yrs


# Convert columns to appropriate types if needed
var1$value <- as.numeric(var1$value)
var1$year <- as.numeric(var1$year)
var1$month <- as.numeric(var1$month)

data_yrs$YEAR1 <- as.numeric(data_yrs$YEAR1)
data_yrs$MEASYEAR <- as.numeric(data_yrs$MEASYEAR)

# Calculate annual means
annual_means <- var1 %>%
  group_by(puid, year) %>%
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")

# Calculate summer means (months 7-9)
summer_means <- var1 %>%
  filter(month >= 7, month <= 9) %>%
  group_by(puid, year) %>%
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")



#### Function to find the mean climate var value and extract results for a single puid
analyze_single_puid <- function(single_puid, data_type = c("annual", "summer")) {
  # Get year range for this puid
  puid_years <- data_yrs %>%
    filter(puid == single_puid)
  
  if(nrow(puid_years) == 0) {
    return(NULL)  # Skip if puid not found in years data
  }
  
  pre_year <- puid_years$PREYR[1]
  start_year <- puid_years$YEAR1[1]
  end_year <- puid_years$MEASYEAR[1]
  
  # Select appropriate data based on data_type
  if(data_type == "annual") {
    data.var <- annual_means %>%
      filter(puid == single_puid, year >= pre_year, year <= end_year)
    value_col <- "annual_mean"
  } else {
    data.var <- summer_means %>%
      filter(puid == single_puid, year >= pre_year, year <= end_year)
    value_col <- "summer_mean"
  }
  
  if(nrow(data.var) < 2) {
    return(data.frame(
      puid = single_puid,
      data_type = data_type,
      pre_year = pre_year,
      start_year = start_year,
      end_year = end_year,
      pre_mean = NA,
      end_mean = NA,
      difference = NA,
      yrs_btwn_visits = yrs_btwn_visits
    ))
  }
  
  
  # Calculate means
  pre_mean <- mean(data.var$mean[data.var$year >= pre_year & data.var$year <= start_year])
  end_mean <- mean(data.var$mean[data.var$year >= start_year & data.var$year <= end_year])
  yrs_btwn_visits <- length(start_year:end_year)
  difference <- end_mean - pre_mean
  
  return(data.frame(
    puid = single_puid,
    data_type = data_type,
    pre_year = pre_year,
    start_year = start_year,
    end_year = end_year,
    pre_mean = pre_mean,
    end_mean = end_mean,
    difference = difference,
    yrs_btwn_visits = yrs_btwn_visits
  ))
}

# Get all unique puids from years data
all_puids <- unique(data_yrs$puid)
indices <- 1:length(all_puids)
# batch_size is for the parallel processing. Each 
#  operation does not take much computing power, so breaking into 
#   batches.  
batch_size <- 300
batches <- split(indices, ceiling(seq_along(indices) / batch_size))

# Analyze each puid for both annual and summer data
plan(multisession, workers = n.cores)

parallel.var.fcn <- function(batch.def, data.type) {
  annual_result <- future_map(batch.def, function(idx) {
    map_df(all_puids[idx], analyze_single_puid, data_type = data.type)
  }) %>% 
    bind_rows()
}

# Takes about three minutes on my computer....
tic()
annual_results <- parallel.var.fcn(batches, "annual")
summer_results <- parallel.var.fcn(batches, "summer")
toc()


# Combine all results
all_results <- list(annual_results = annual_results, summer_results = summer_results)

# Save results
write_rds(all_results, file.path(DATA.LOC, "Climate_plot_results.rds"))

#all_results <- read_rds(file.path(DATA.LOC, "Climate_plot_results.rds"))



