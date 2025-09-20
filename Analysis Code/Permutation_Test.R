## Monte Carlo permutation analysis ##

source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))

library(furrr)
library(parallel)
library(tictoc)  # For development, timing routines.
library(data.table)
library(abind) # for combining matrices into arrays
library(simctest) # Gandy 2009 permutation reduction approach

## Constants

# Number of maximum MC permutation iterations:
PERM.ITER.N <- 10




## Loading data to be used.  Assuming analysis.arrays.RDS exist. If not, Data_Prep_Analysis will need to be run.
if(file.exists(paste0(DATA.LOC, "analysis.arrays.zip")) == TRUE){
  arrays.use <- read_rds(unzip(paste0(DATA.LOC, "analysis.arrays.zip"), "analysis.arrays.RDS"))
  if(file.exists("analysis.arrays.RDS")) {
    file.remove("analysis.arrays.RDS") # Too big for GitHub
  }  
} else {
  print("File analysis.arrays.zip not found")
}


# Number of MC permutation iterations:
arrays.grow <- arrays.use$arrays.grow
arrays.mort <- arrays.use$arrays.mort
PlotDat <- arrays.use$PlotDat
# BS.N <- 100


# Finding the analysis results to pull out the minimum number of significant findings
#  per contrast. This is used in the Gandy simctest() function at the end.
contrast.min <- rep(NA, 6)
for(k in 1:2) {
plt.dat <- readRDS(paste0(save.loc.fcn(k), "Processed_6Domain_Data.RDS"))

sig.AvB2 <- length(plt.dat$AvB2$Species[plt.dat$AvB2$significant == TRUE])
sig.A_LMH2 <- length(plt.dat$A_LMH2$Species[plt.dat$A_LMH2$significant == TRUE])
sig.B_LMH2 <- length(plt.dat$B_LMH2$Species[plt.dat$B_LMH2$significant == TRUE])

contrast.min[(k-1) * 3 + 1:3] <- c(sig.AvB2, sig.A_LMH2, sig.B_LMH2)
}
contrast.min.val <- min(contrast.min)                # Minimum number of significant contrasts
contrast.min.position <- which(contrast.min == min(contrast.min)) # Position of the minimum number (id which contrast)


# Need climate names for files and axes.
clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE) %>%
  filter(filename == CLIM.VAR.USE)

var1 <- clim.names$values[grep("pre", clim.names$values)]
var.delt <- clim.names$values[grep("d", clim.names$values)]


var.label <- clim.names$label[clim.names$values == var1]
#var.filename <- clim.names$filename[clim.names$values == var1]
var.delt.label <- clim.names$label[clim.names$values == var.delt]


## Setting up some matrices to collect values while the Gandy procedure is underway.
#  We can see what the levels of false positives are for other contrasts.
perm.grow <- perm.mort <- matrix(data = NA, nrow = 1, ncol = 3, dimnames = list(NULL, c("AvB", "Above", "Below")))
perm.output <- list(perm.grow = perm.grow, perm.mort = perm.mort)



## Here's the main permutation function.  
perm.test.fcn <- function() {
  
  gandy.out <- c(rep(NA, 6))
  for(k in 1:length(ANALYSIS.TYPE)) {  # 1 = grow, 2 = mortality
    
    # Obtaining the data to work with: Grabbing the list objects from above
    mort.grow.dat <- get(paste0("arrays.", ANALYSIS.TYPE[k])) 
    
    vals_dat <- mort.grow.dat$vals_dat
    all_dat <- mort.grow.dat$all_dat
    #domain.array <- mort.grow.dat$domain.array
    state.array <- mort.grow.dat$state.array
    state.n <- mort.grow.dat$state.n
    domain.matrix <- mort.grow.dat$domain.matrix
    centroid.array <- mort.grow.dat$centroid.array
    domain.n <- mort.grow.dat$domain.n 
    quant.lims <- mort.grow.dat$quant.lims
    quant.lims.delta <- mort.grow.dat$quant.lims.delt
    
    
    ## First, adjusting the species list
    #spp.id <- paste0("X", spp.list)     # Can use SEL.SPP
    spp.list <- as.numeric(gsub("X", "", SEL.SPP))
    
    strat.num <- length(unique(vals_dat$stratum))
    
    n.spp <- length(SEL.SPP)
    
    
    # Strata for bootstrap resampling procedure
    strata <- unique(vals_dat$stratum)
    strata.num <- vals_dat %>% dplyr::select(stratum, puid) %>%
      ungroup() %>% 
      mutate(row.id = row_number())
    
    ## THIS IS WHERE THE DOMAIN ASSIGNMENTS ARE SAMPLED FOR THE PERMUTATION:
    domain.samp.array <- samp.domain.array.fcn(domain.matrix)
      
      
    ## Quantile estimates for mortality or growth, across select species.
    domain.index <- 1:n_domain

    
    
    ## Generating the bootstrap values: 
    # 200 iterations = 104.99 s, or 1.75 min, or 29.2 hours
    # 1000 iterations = 426.21 s, or 7.1 min, or 118.3 hrs or 5 days....
    bootstrap_results <- generate_bootstrap_array.fcn(
      vals.dat = vals_dat,
      all.dat = all_dat,
      domain.array = domain.samp.array, # This now uses the sampled array
      domain.n = domain.n,
      selected.spp = SEL.SPP,
      n_iter = BS.N, 
      strata.num = strata.num,
      PlotDat = PlotDat
    )
    
    # Finding and saving domain summaries
    domain.summaries <- domain.index %>% 
      purrr::map(\(d) domain.sum.fcn(bootstrap_results, d, domain_n = domain.n)) %>%
      do.call(rbind, .) %>%
      arrange(Species, Domain)
    
    
    ## The code in this section is dedicated to preparing data files for plotting the differences in six-domain results.
    ## This process has two steps: first, the domain matrix pairs are subtracted from one another 
    #  using the domain.diff.fcn function.  Then, using the differenced matrices, the matrices
    #  are processed and prepared for figure creation.  Along the way the species' Latin names are
    #  prepped for inclusion.  
      
      # Domains compared: Above vs. Below threshold 
      A.vec <- c("AL", "AM", "AH")
      B.vec <- c("BL", "BM", "BH")
      
      AvB <- bind_rows(map2(A.vec, B.vec, domain.diff.fcn, results.array = bootstrap_results, domain_n = domain.n)) %>%
        arrange(Species, Domain)
      
      
            # Domains compared: Above threshold, Low/Med/High comparisons
      A_LMH.vec1 <- c("AH", "AM", "AH")
      A_LMH.vec2 <- c("AL", "AL", "AM")
      
      A_LMH <- bind_rows(map2(A_LMH.vec1, A_LMH.vec2, domain.diff.fcn, results.array = bootstrap_results, domain_n = domain.n)) %>%
        arrange(Species, Domain)
      
      # Domains compared: Below threshold, Low/Med/High comparisons
      B_LMH.vec1 <- c("BH", "BM", "BH")
      B_LMH.vec2 <- c("BL", "BL", "BM")
      
      B_LMH <- bind_rows(map2(B_LMH.vec1, B_LMH.vec2, domain.diff.fcn, results.array = bootstrap_results, domain_n = domain.n)) %>%
        arrange(Species, Domain)
      

    n.sig.AvB <- sig.n.fcn(AvB)
    n.sig.A_LMH <- sig.n.fcn(A_LMH)
    n.sig.B_LMH <- sig.n.fcn(B_LMH)
    
    perm.output[[k]] <<- perm.output[[k]] %>% rbind(matrix(c(n.sig.AvB, n.sig.A_LMH, n.sig.B_LMH), ncol = 3))
    gandy.out[(k-1) * 3 + 1:3] <- c(n.sig.AvB, n.sig.A_LMH, n.sig.B_LMH)

  }   # end k 
  return(gandy.out)
}

# The Gandy simctest needs a function that makes a comparison.  In this case the 
# comparison is to see if each iteration produces a TRUE value.
gen <- function(){perm.test.fcn()[contrast.min.position] >= contrast.min.val}

# Setting up parallel computing. See Global.R for n.cores.
plan(multisession, workers = n.cores) 


tic()
simctest(gen, maxsteps = PERM.ITER.N)  # Here's the permutation test.
toc()

# Closing parallel workers
plan(sequential)



#1781.76 sec elapsed, 30 min for 10 iterations of the permutation, 100 bootstraps
# per permutation.  So, 3 min/iteration.  It does not appear that simctest is 
#  making the process unduly long. The process does that all by itself : ) 


#> perm.output
#$perm.grow
#AvB Above Below
#[1,]  NA    NA    NA
#[2,]   6     6     6
#[3,]   9     9     1
#[4,]   3     5     4
#[5,]   6     8     5
#[6,]   4     6     6
#[7,]   8     7     5
#[8,]   3     9     6
#[9,]   4     2     7
#[10,]   3     3     2
#[11,]  14     8    14

#$perm.mort
#     AvB Above Below
#[1,]  NA    NA    NA
#[2,]   4     7     7
#[3,]   7     6     3
#[4,]   4     5     9
#[5,]   2     5     6
#[6,]   6     4     7
#[7,]   1     8     9
#[8,]   5     5     8
#[9,]   4     6     4
#[10,]   6     4     2
#[11,]   6     6     7
















