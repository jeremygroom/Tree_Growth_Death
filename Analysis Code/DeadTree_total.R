## This code runs A. Yost's mortality code on our tree data.  The purpose is to find where our 
#   code and his diverge.  


Tree.use <- Tree1 %>% filter(STATECD == 41) %>%
  select("PLOT_FIADB", "TREE", "SUBP") %>% 
  mutate(usetree = 1) %>% 
  distinct()




library(plyr)
library(stringr)
library(dplyr)

library(RSQLite)
library(furrr)   # parallel processing for purrr::map()-like functions
library(parallel)  # for finding out how many cores we have
library(tictoc)   # evaluate how long something takes

con <- dbConnect(RSQLite::SQLite(), "G:/My Drive/Consulting Practice/Contracts/ODF_SAS_to_R_2023/OR_FIA/SQLite_FIADB_OR.db")

plot<-dbReadTable(con,"PLOT",row.names=FALSE)
tree<-dbReadTable(con,"TREE",row.names=FALSE)
cond<-dbReadTable(con,"COND",row.names=FALSE)
peu<-dbReadTable(con,"POP_ESTN_UNIT", row.names=FALSE)
pev<-dbReadTable(con,"POP_EVAL",row.names=FALSE)
pet<-dbReadTable(con, "POP_EVAL_TYP",row.names=FALSE)
popstrm<-dbReadTable(con, "POP_STRATUM",row.names=FALSE)
ppsa<-dbReadTable(con, "POP_PLOT_STRATUM_ASSGN",row.names=FALSE)


PLOT <- select(plot, CN, PLOT, PLOT_STATUS_CD,MACRO_BREAKPOINT_DIA,ECOSUBCD, MEASYEAR, LAT, LON, REMPER)
   # CN = plot unique ID, Plot = phase 2 plot number, PLOT_STATUS_CD = status code (1 = some forest accessible), MACRO_BREAKPOINT_DIA = 24 or 30 inches measured for macroplots,
  #  ECOSUBCD = ecological subsection code, not sure why useful, MEASYEAR = yr measurement completed (not necessarily INVYR), REMPER = remeasure period in yrs and decimal yrs
COND <- select(cond, PLT_CN, CONDID, CONDPROP_UNADJ, PROP_BASIS, COND_STATUS_CD, OWNGRPCD)
  # PLT_CN = plot seq. num, CONDID = Condition class number - conditions across plot, CONDPROPR_UNADJ = Unadjusted proportion of a plot in a condition - for area calcs, PROP_BASIS =
#    type of subplots installed when sampled, COND_STATUS_CD = basic land classifcation, OWNGRPCD = broad ownership classification.
TREE <- select(tree, CN, PLOT, INVYR, SUBP, PREV_TRE_CN, PLT_CN, CONDID, TREE, STATUSCD, DRYBIO_AG, 
               CARBON_AG, TPA_UNADJ, DIA, SPCD, STANDING_DEAD_CD, PREV_STATUS_CD,AGENTCD,RECONCILECD) %>%
  # CN = unique tree ID, PLOT = Phase2 plot number, INVYR, SUBP = subplot, PREV_TRE_CN = previous CN, PLT_CN, CONDID = condition class #, TREE, STATUSCD = live, dead, removed, 
  #  DRYBIO_AG = aboveground dry biomass, CARBON_AG, TPA_UNADJ = TPA without accounting for things, DIA, SPCD, STANDING_DEAD_CD = yes or no, PREV_STATUS_CD, AGENTCD, RECONCILED = 
  #    ingrowth, throughgrowth, missed, moved, etc.
  left_join(Tree.use, by = c("PLOT" = "PLOT_FIADB", "TREE", "SUBP")) %>% filter(usetree == 1)
POP_ESTN_UNIT <- select(peu, CN, ESTN_UNIT, EVAL_CN, AREA_USED, P1PNTCNT_EU)
  # CN = Population estimation unit #, ESTN_UNIT = Geographic area being stratified, EVAL_CN = eval # (next), AREA_USED = area used to calc expansion factors, P1PNTCNT_EU = # pixels
#   in the estimation unit (phase 1)
POP_EVAL <- select(pev, EVALID, EVAL_GRP_CN, ESTN_METHOD, CN, END_INVYR, REPORT_YEAR_NM)
   #    EVALID = This plus RSCD code used to identify field plots and associated Phase 1 summary data.  EVAL_GRP_CN = Key links pop eval record to pop eval group record,
  #     ESTN_METHOD = Estimation method (SRS, Stratified, double sampling, post-strat, subsampling), CN, END_INVYR, REPORT_YEAR_NM = data collection yrs.
POP_EVAL_TYP <- select(pet, EVAL_TYP, EVAL_CN)
     # EVAL_TYP = includes EXPMORT which is plots used for tree mortality estimates, EVAL_CN = link to POP_EVAL
POP_PLOT_STRATUM_ASSGN <- select(ppsa, STRATUM_CN, PLT_CN)
     #  STRATUM_CN = stratum sequence #
POP_STRATUM <- select(popstrm, STRATUMCD, ESTN_UNIT_CN, EXPNS, P2POINTCNT, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR, CN, P1POINTCNT)
     # STRATUMCD = 3214 stratum numbers, ESTN_UNIT_CN = link to pop estimation unit number, EXPNS = expansion factor, P2POINTCNT = # field plots within stratum, 
    #    ADJ_FACTOR_MICR = adjust factor for microplot in case of partial sampling, ADJ_FACTOR_SUBP = same for subplots, ADJ_FACTOR_MACR = macroplots


data <- PLOT %>%
  ## Add a PLT_CN column for easy joining
  mutate(PLT_CN = CN) %>%
  ## Join COND & TREE
  left_join(COND, by = 'PLT_CN') %>%
  left_join(TREE, by = c('PLT_CN', 'CONDID')) %>%
  ## Population tables
  left_join(POP_PLOT_STRATUM_ASSGN, by = 'PLT_CN') %>%
  left_join(POP_STRATUM, by = c('STRATUM_CN' = 'CN')) %>%
  left_join(POP_ESTN_UNIT, by = c('ESTN_UNIT_CN' = 'CN')) %>%
  left_join(POP_EVAL, by = c('EVAL_CN' = 'CN')) %>%
  left_join(POP_EVAL_TYP, by = 'EVAL_CN')


df2003<-data %>% 
  filter(EVALID==412003,
         EVAL_TYP=="EXPMORT")

###########################################
# Estimating Mortality Rate

#Live Trees PREV_STATUS_CD ==1 & STATUSCD==1

#Move these two columns (for better organization) that aren't used in the mortality calculation
df2003<-df2003 %>% 
  relocate(c(DRYBIO_AG, CARBON_AG),.after = RECONCILECD)


df2003$TPA_UNADJt1<-0
df2003$DIA1<-0
df2003<-df2003 %>% 
  relocate(DIA1,.after = DIA)
library(svMisc)
#Move tpa_unadj values for trees with DIA>=5 from inventory 1 to df2003thpl.
#This code takes about .5 hrs to run


# Super long time....
# for(i in 1:nrow(df2003)){
#   
#   if(is.na(df2003$PREV_TRE_CN[i])){next} else if(df2003$SPCD[i]!=242) {next} else{a<-df2003$PREV_TRE_CN[i]}
#   #if_else(!is.na(df2003b$PREV_TRE_CN[i]) & df2003b$SPCD[i]==242, df2003b$TPA_UNADJt1[i]<-tree$TPA_UNADJ[which(tree$CN==df2003b$PREV_TRE_CN[i])],0)
#   c<-0
#   for(j in 1:nrow(tree)){
#     
#     if(!is.na(tree$CN[j]) & tree$CN[j]==a & tree$DIA[j]>=5)
#     {df2003$TPA_UNADJt1[i]<-tree$TPA_UNADJ[j]; df2003$DIA1[i]<-tree$DIA[j];c<-1}
#     if(c==1){break}
#   }
#   progress(i,nrow(df2003))
# }


# number of rows of df2003
nrow(df2003)  # 347421

# number of rows of df2003 without NAs or SPCD = 242
nrow(df2003 %>% filter(is.na(df2003$PREV_TRE_CN) == FALSE, df2003$SPCD == 242)) # 3743

df2003.2 <- df2003 %>% filter(is.na(df2003$PREV_TRE_CN) == FALSE, df2003$SPCD == 242) # getting rid of those unwanteds

# How many unique rows of df2003's PREV_TRE_CN
length(unique(df2003.2$PREV_TRE_CN)) # 3743 lines are unique, so no NAs or duplicates

# Which rows are duplicated and why? (Turns out there are no duplicates)
df2 <- df2003.2 %>% dplyr::select(PREV_TRE_CN) %>% 
  group_by(PREV_TRE_CN) %>% 
  summarize(n = n())  
max(df2$n) # = 1

# The tree data are only wanted for adjusting the df2003 data.  So, let's create a smaller tree dataset.
ptc <- unique(df2003.2$PREV_TRE_CN)  # Grabbing the codes we want ('a')
tree2 <- tree %>% filter(CN %in% ptc & DIA >= 5) # This has 255k rows instead of 708k
tree3 <- tree2 %>% dplyr::select(CN, TPA_UNADJ, DIA) %>%
  rename("TPA_UNADJ_tree" = "TPA_UNADJ", "DIA_tree" = "DIA")

# Bringing it back to the whole dataset... notice that df2003 is used as the starting point, NAs and all
df2003.3 <- df2003 %>% left_join(tree3, by = c("PREV_TRE_CN" = "CN")) %>%
  mutate(TPA_UNADJt1 = ifelse(is.na(TPA_UNADJ_tree) == FALSE, TPA_UNADJ_tree, TPA_UNADJt1),
         DIA1 = ifelse(is.na(DIA_tree) == FALSE, DIA_tree, DIA1))

df2003 <- df2003.3

df2003 <- df2003 %>%
  #  df1001 <- df1001 %>%
  mutate(
    ## AREA
    aAdj = case_when(
      ## When NA, stay NA
      is.na(PROP_BASIS) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
      ## Otherwise, use the subpplot value
      PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
    ## TREE
    tAdj = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA1) ~ NA, #ADJ_FACTOR_SUBP,
      ## When DIA is less than 5", use microplot value
      DIA1 < 5 ~ ADJ_FACTOR_MICR,
      ## When DIA is greater than 5", use subplot value
      MACRO_BREAKPOINT_DIA==30 & DIA1 >= 5 & DIA1 <30 ~ ADJ_FACTOR_SUBP,
      MACRO_BREAKPOINT_DIA==24 & DIA1 >= 5 & DIA1 <24 ~ ADJ_FACTOR_SUBP,
      ## DIA is greater than macro breakpoint diameter
      MACRO_BREAKPOINT_DIA==30 & DIA1 >=30 ~ ADJ_FACTOR_MACR,
      MACRO_BREAKPOINT_DIA==24 & DIA1 >=24 ~ ADJ_FACTOR_MACR
    ))

#Create the live and dead tree indicators. including or excluding fire and harvest trees
df2003$aDI<-0
df2003$LtDI<-0 #Live tree indicator
df2003$DtDI<-0 #Dead tree indicator
df2003$Nlive<-0
df2003$Ndead<-0

#plot_status_cd=1 is sampled with forest present, 2-sampled with no forest, 3 nonsampled
#cond_status_cd=1 is accessible forest land. 2-other wooded land, 3-noncensus water, 4-census water, 5-not sampled.

df2003<-df2003 %>% 
  mutate(aDI=if_else(COND_STATUS_CD == 1, 1,0),
         
         LtDI=ifelse(DIA1>=5 & SPCD==242 &  STATUSCD==1 & !is.na(PREV_STATUS_CD) & PREV_STATUS_CD==1,1*aDI,0), #DIA1 = first inventory diameter
         
         DtDI=ifelse(TPA_UNADJt1>0 &DIA1>=5 & SPCD==242 & STATUSCD==2 & !is.na(PREV_STATUS_CD) & PREV_STATUS_CD==1 & AGENTCD !=30 & AGENTCD!=80,1*aDI,0),
         
         Nlive= (LtDI*TPA_UNADJt1*tAdj*EXPNS),
         Ndead=(DtDI*TPA_UNADJt1*tAdj*EXPNS)) %>% 
  
  #organize columns  
  relocate(c(aAdj,tAdj,aDI,LtDI,DtDI,Nlive,Ndead),.after = RECONCILECD) %>% 
  relocate(STATUSCD,.before =PREV_STATUS_CD) %>% 
  relocate(TPA_UNADJt1,.after = TPA_UNADJ) %>% 
  relocate(TPA_UNADJ_tree,.after=TPA_UNADJt1) %>% 
  relocate(DIA_tree,.after=DIA1)
#build a new dataframe grouped by 

#df2003$DtDI<-if_else(df2003$AGENTCD==80,0,df2003$DtDI)

Mrate <- df2003 %>%
  rename("PLOT" = "PLOT.x") %>%
  group_by(PLOT,REMPER) %>%
  summarize(NdeadPlot=sum(Ndead),
            NdeadPlotRP=NdeadPlot/first(REMPER),
            NlivePlot=sum(Nlive))

Mrate$NdeadPlot<-if_else(is.na(Mrate$NdeadPlot),0,Mrate$NdeadPlot)
Mrate$NdeadPlotRP<-if_else(is.na(Mrate$NdeadPlotRP),0,Mrate$NdeadPlotRP)
Mrate$NlivePlot<-if_else(is.na(Mrate$NlivePlot),0,Mrate$NlivePlot)

Mrate$Total<-Mrate$NdeadPlotRP+Mrate$NlivePlot
Mrate242<-mean(Mrate$NdeadPlotRP)/(mean(Mrate$Total))*100

Mrate242




rm(tree2,tree3,df2003.3,df2003.2)

