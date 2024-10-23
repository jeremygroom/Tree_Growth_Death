#### Obtaining and processing PLOT, TREE, and COND tables from FIA databases for Oregon, Washington, and California #####


statelist <- c("CA", "WA", "OR")

# This function will extract the zipped databases, reduce their tables, and provide the tables as a list.
db.extract.fcn <- function(state, db.loc) {
  db.use <- dbConnect(RSQLite::SQLite(), unzip(paste0(db.loc, "SQLite_FIADB_", state, ".zip"), paste0("SQLite_FIADB_", state, ".db")))
  plot <- dbReadTable(db.use, "PLOT", row.names = FALSE) %>% select(CN, STATECD, PLOT, INVYR, MEASYEAR, REMPER)
  tree <- dbReadTable(db.use, "TREE", row.names = FALSE) %>% select(CN, STATECD, PLOT, TREE, SUBP, INVYR, AGENTCD)
  cond <- dbReadTable(db.use, "COND", row.names = FALSE) %>% select(CN, STATECD, PLOT, PLT_CN, INVYR, STATECD, PLOT, MICRPROP_UNADJ, MACRPROP_UNADJ, SITECLCD)
  
  list.out <- list(plot = plot, tree = tree, cond = cond)
  dbDisconnect(db.use)  # Need to disconnect the databases so that it can be reconnected to a different state.
  file.remove(paste0("SQLite_FIADB_", state, ".db"))
  return(list.out)
  
}

all.states <- map(statelist, db.extract.fcn, db.loc = SQL.LOC)  # Appling the function to all three states

extract.plots <- map(all.states , ~.[["plot"]])   # Binding by rows all of the plot data across states
plot_vals <- bind_rows(extract.plots)

extract.trees <- map(all.states , ~.[["tree"]]) 
tree_vals <- bind_rows(extract.trees)

extract.cond <- map(all.states , ~.[["cond"]]) 
cond_vals <- bind_rows(extract.cond)

# This RDS file is written, zipped, and then removed.
write_rds(list(plot_vals = plot_vals, tree_vals = tree_vals, cond_vals = cond_vals),
          paste0(DATA.LOC, "Addl_PlotTreeInfo.rds"))

# This is the zip file the rest of the code will operate with.
zip(zipfile = paste0(DATA.LOC, "Addl_PlotTreeInfo.zip"), files = paste0(DATA.LOC, "Addl_PlotTreeInfo.rds"),
    flags = '-r9Xj' )  # This weird bit keeps the parent directory from being included in the zip folder

# Deleting written RDS file (~80 MB)
if (file.exists(paste0(DATA.LOC, "Addl_PlotTreeInfo.rds")) == TRUE) {
  file.remove(paste0(DATA.LOC, "Addl_PlotTreeInfo.rds"))
}
