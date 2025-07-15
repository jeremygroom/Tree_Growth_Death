#### Obtaining and processing PLOT, TREE, and COND tables from FIA databases for Oregon, Washington, and California #####


statelist <- c("CA", "WA", "OR")

# This function will extract the zipped databases, reduce their tables, and provide the tables as a list.
db.extract.fcn <- function(state, db.loc) {
  db.use <- dbConnect(RSQLite::SQLite(), unzip(paste0(db.loc, "SQLite_FIADB_", state, ".zip"), paste0("SQLite_FIADB_", state, ".db")))
  plot <- dbReadTable(db.use, "PLOT", row.names = FALSE) %>% select(CN, STATECD, COUNTYCD, PLOT, PLOT_STATUS_CD,MACRO_BREAKPOINT_DIA,MEASYEAR, REMPER, LAT, LON, ELEV, INTENSITY)
  tree <- dbReadTable(db.use, "TREE", row.names = FALSE) %>% select(PLT_CN, PLOT, STATECD, SUBP, CONDID, PREVCOND, TREE, PREV_TRE_CN, STATUSCD, DRYBIO_AG, CARBON_AG, TPA_UNADJ, 
                                                                    DIA, PREVDIA, SPCD, STANDING_DEAD_CD, PREV_STATUS_CD,AGENTCD, RECONCILECD, INVYR) %>%
    mutate(AGENTCD2 = ifelse(is.na(AGENTCD), 0, AGENTCD),
           DIA2 = ifelse(is.na(DIA) & AGENTCD2 > 0, 0, DIA)) %>% 
    filter(!RECONCILECD %in% 5:9,    # Removing 1120 trees w/ RECONCILECD values
           !is.na(TREE),             # Remove entries where no TREE ID number
           !is.na(PREV_STATUS_CD),   # Remove entries where no previous STATUSCD
           PREV_STATUS_CD != 2,     # Remove entries where tree previously found dead (dead --> dead)
           STATUSCD != 0,           # Remove entries where tree is not presently in sample
           AGENTCD2 != 80,
           PREVDIA >= 5,             # Remove microplot trees that became ingrowth.
           DIA2 == 0 | DIA2 >= 5     # Keep DIA = NA in case the tree burned and died, no DIA.  Remove the small trees.
    )
  cond <- dbReadTable(db.use, "COND", row.names = FALSE) %>% select(PLT_CN, PLOT, STATECD, CONDID, CONDPROP_UNADJ, PROP_BASIS, COND_STATUS_CD, OWNGRPCD) #%>%
    #filter(COND_STATUS_CD == 1)     # Only forested plots
  sccm <- dbReadTable(db.use, "SUBP_COND_CHNG_MTRX", row.names=FALSE) %>% select(PLT_CN, SUBP, STATECD, CONDID, PREVCOND, SUBPTYP_PROP_CHNG, SUBPTYP) %>% filter(SUBPTYP != 2)
  peu <- dbReadTable(db.use,"POP_ESTN_UNIT", row.names=FALSE) %>% select(CN, STATECD, ESTN_UNIT, EVAL_CN, AREA_USED, P1PNTCNT_EU)
  pev <- dbReadTable(db.use,"POP_EVAL",row.names=FALSE) %>% select(EVALID, STATECD, EVAL_GRP_CN, ESTN_METHOD, CN, END_INVYR, REPORT_YEAR_NM)
  peg <- dbReadTable(db.use,"POP_EVAL_GRP", row.names = FALSE)  # Not sure this one is used
  pet <- dbReadTable(db.use,"POP_EVAL_TYP", row.names=FALSE) %>% select(EVAL_TYP, EVAL_CN)
  ppsa <- dbReadTable(db.use, "POP_PLOT_STRATUM_ASSGN",row.names=FALSE) %>% 
    select(STRATUM_CN, PLT_CN, EVALID) %>%
    mutate(id_eval  = substr(EVALID, nchar(EVALID) - 1, nchar(EVALID))) %>%
    filter(id_eval == "03") %>% # Plots sampled for growth and mortality
    dplyr::select(-id_eval) %>%
    group_by(PLT_CN) %>%
    summarize(STRATUM_CN = first(STRATUM_CN, .multi = FAIL))
  popstrm <- dbReadTable(db.use, "POP_STRATUM",row.names=FALSE) %>% select(STRATUMCD, STATECD, ESTN_UNIT_CN, EXPNS, P2POINTCNT, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR, CN, P1POINTCNT)
  
  
  list.out <- list(plot = plot, tree = tree, cond = cond, sccm = sccm, peu = peu, pev = pev, peg = peg, pet = pet, ppsa = ppsa, popstrm = popstrm)
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

extract.sccm <- map(all.states , ~.[["sccm"]]) 
sccm_vals <- bind_rows(extract.sccm)

extract.peu <- map(all.states , ~.[["peu"]]) 
peu_vals <- bind_rows(extract.peu)

extract.pev <- map(all.states , ~.[["pev"]]) 
pev_vals <- bind_rows(extract.pev)

extract.peg <- map(all.states , ~.[["peg"]]) 
peg_vals <- bind_rows(extract.peg)

extract.pet <- map(all.states , ~.[["pet"]]) 
pet_vals <- bind_rows(extract.pet)

extract.ppsa <- map(all.states , ~.[["ppsa"]]) 
#extract.ppsa <- lapply(extract.ppsa, function(df) {
#  df$MODIFIED_DATE <- as.character(df$MODIFIED_DATE)  # Needed consistency in column type. Differed across states.
#  return(df) 
#})
ppsa_vals <- bind_rows(extract.ppsa)

extract.popstrm <- map(all.states , ~.[["popstrm"]]) 
popstrm_vals <- bind_rows(extract.popstrm)



# This RDS file is written, zipped, and then removed.
write_rds(list(plot_vals = plot_vals, tree_vals = tree_vals, cond_vals = cond_vals, sccm_vals = sccm_vals,
               peu_vals = peu_vals, pev_vals = pev_vals, peg_vals = peg_vals,
               pet_vals = pet_vals, ppsa_vals = ppsa_vals, popstrm_vals = popstrm_vals),
          paste0(DATA.LOC, "Addl_PlotTreeInfo.rds"))

rm(extract.popstrm, extract.ppsa, extract.pet, extract.peg, extract.peu, extract.cond, extract.trees, extract.plots,
   popstrm_vals, ppsa_vals, pet_vals, peg_vals, peu_vals, cond_vals, tree_vals, plot_vals, all.states)


# This is the zip file the rest of the code will operate with.
zip(zipfile = paste0(DATA.LOC, "Addl_PlotTreeInfo.zip"), files = paste0(DATA.LOC, "Addl_PlotTreeInfo.rds"),
    flags = '-r9Xj' )  # This weird bit keeps the parent directory from being included in the zip folder

# Deleting written RDS file (~200 MB)
if (file.exists(paste0(DATA.LOC, "Addl_PlotTreeInfo.rds")) == TRUE) {
  file.remove(paste0(DATA.LOC, "Addl_PlotTreeInfo.rds"))
}
