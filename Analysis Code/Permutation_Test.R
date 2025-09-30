## Monte Carlo permutation analysis ##

library(furrr)
library(parallel)
library(future) # Involved in parallel computing
library(tictoc)  # For development, timing routines.
library(data.table)
library(matrixStats) # For fast iteration processing
library(abind) # for combining matrices into arrays
library(simctest) # Gandy 2009 permutation reduction approach

source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))


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
#  per contrast and the distance of the CI from zero.
sig.dist <- NULL
contrast.min <- rep(NA, 6)
for(k in 1:2) {
plt.dat <- readRDS(paste0(save.loc.fcn(k), "Processed_6Domain_Data.RDS"))

# Counting the number of sig findings by column for k (growth, mortality)
sig.AvB2 <- length(plt.dat$AvB2$Species[plt.dat$AvB2$significant == TRUE])
sig.A_LMH2 <- length(plt.dat$A_LMH2$Species[plt.dat$A_LMH2$significant == TRUE])
sig.B_LMH2 <- length(plt.dat$B_LMH2$Species[plt.dat$B_LMH2$significant == TRUE])

contrast.min[(k-1) * 3 + 1:3] <- c(sig.AvB2, sig.A_LMH2, sig.B_LMH2)

# Finding the CI distance from zero by column for k (growth, mortality)
sig.dist.AvB <- sig.dist.fcn(plt.dat$AvB)
sig.dist.A_LMH <- sig.dist.fcn(plt.dat$A_LMH)
sig.dist.B_LMH <- sig.dist.fcn(plt.dat$B_LMH)

sig.dist.k <- data.frame(sig.dist.AvB, sig.dist.A_LMH, sig.dist.B_LMH) %>%
  mutate(order.c = seq_along(1:n())) %>%
  relocate(order.c)
if(is.null(sig.dist)){
sig.dist <- sig.dist.k %>% cbind(Species = plt.dat$AvB$Species) %>% relocate(Species)
} else {
  if(all(row.names(sig.dist.k) == row.names(sig.dist)) & all(sig.dist.k$order.c == sig.dist$order.c)) {
    sig.dist <- left_join(sig.dist, sig.dist.k, by = "order.c") 
    colnames(sig.dist)[3:ncol(sig.dist)] <- paste0(c("Domain", "sig.dist"), rep(1:6, each = 2))  
    
  } else {
    print("Error: misaligned species or comparison numbers.")
  }
}

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
perm.test.fcn <- function(arrays.grow, arrays.mort, PlotDat) {
  gandy.out <- numeric(6)
  sig.tests.out <- NULL#matrix(NA, nrow = ncol(arrays.grow$domain.n) * 3, ncol = 14 )
  #colnames(sig.tests.out) <- c("Species", "G.AvB.domain", "G.AvB.sig", "G.A_LMH.domain", "G.A_LMH.sig", "G.B_LMH.domain", "G.B_LMH.sig",
   #                         "M.AvB.domain", "M.AvB.sig", "M.A_LMH.domain", "M.A_LMH.sig", "M.B_LMH.domain", "M.B_LMH.sig")
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
    strata.num <- vals_dat %>% dplyr::select(stratum, puid, w, n_h.plts) %>%
      ungroup() %>% 
      mutate(row.id = row_number())
    
    ## THIS IS WHERE THE DOMAIN ASSIGNMENTS ARE SAMPLED FOR THE PERMUTATION:
    domain.samp.array <- samp.domain.array.fcn(domain.matrix)

    ## Quantile estimates for mortality or growth, across select species.
    domain.index <- 1:n_domain

    
    
    ## Generating the bootstrap values: 
    # 200 iterations = 104.99 s, or 1.75 min, or 29.2 hours
    # 1000 iterations = 426.21 s, or 7.1 min, or 118.3 hrs or 5 days....
    
    
    # 100 iterations, 217 sec / 3.5 min using future_map
    # 100 iterations, 227 s/3.8 min using map (no parallel processing): 
    # 200 iterations, 447.29 s/ 7.45 min using map 
    
    #plan(multisession, workers = 5)
    
 
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
        arrange(Species) %>%
        mutate(Significant = ifelse(LCI.95 < 0 & UCI.95 < 0, 1, 
                                    ifelse(LCI.95 > 0 & UCI.95 > 0, 1, 0)))
      
      
            # Domains compared: Above threshold, Low/Med/High comparisons
      A_LMH.vec1 <- c("AH", "AM", "AH")
      A_LMH.vec2 <- c("AL", "AL", "AM")
      
      A_LMH <- bind_rows(map2(A_LMH.vec1, A_LMH.vec2, domain.diff.fcn, results.array = bootstrap_results, domain_n = domain.n)) %>%
        arrange(Species) %>%
        mutate(Significant = ifelse(LCI.95 < 0 & UCI.95 < 0, 1, 
                                    ifelse(LCI.95 > 0 & UCI.95 > 0, 1, 0)))
      
      # Domains compared: Below threshold, Low/Med/High comparisons
      B_LMH.vec1 <- c("BH", "BM", "BH")
      B_LMH.vec2 <- c("BL", "BL", "BM")
      
      B_LMH <- bind_rows(map2(B_LMH.vec1, B_LMH.vec2, domain.diff.fcn, results.array = bootstrap_results, domain_n = domain.n)) %>%
        arrange(Species) %>%
        mutate(Significant = ifelse(LCI.95 < 0 & UCI.95 < 0, 1, 
                                    ifelse(LCI.95 > 0 & UCI.95 > 0, 1, 0)))
      

    n.sig.AvB <- sig.n.fcn(AvB)
    n.sig.A_LMH <- sig.n.fcn(A_LMH)
    n.sig.B_LMH <- sig.n.fcn(B_LMH)
    

    # This code populates sig.tests.out.  The function sig.dist.fcn is structured to
    #   enable left joins and verify that the results will be correctly interpreted
    #   relative to the observed results, sig.dif.  
    sig.dist.use <- list(
      sig.dist.AvB = sig.dist.fcn2(AvB, 1 + (k-1)*3),
      sig.dist.A_LMH = sig.dist.fcn2(A_LMH, 2 + (k-1)*3),
      sig.dist.B_LMH = sig.dist.fcn2(B_LMH, 3 + (k-1)*3)
    )
    
    if(k == 1) {
      sig.tests.out <- purrr::reduce(sig.dist.use, dplyr::left_join, by = c("Species", "order.c"))
    } else {
      k2.tests.out <- purrr::reduce(sig.dist.use, dplyr::left_join, by = c("Species", "order.c"))
      sig.tests.out <- left_join(sig.tests.out, k2.tests.out, by = c("Species", "order.c"))
    }
    
#    perm.output[[k]] <<- perm.output[[k]] %>% rbind(matrix(c(n.sig.AvB, n.sig.A_LMH, n.sig.B_LMH), ncol = 3))
    gandy.out[(k-1) * 3 + 1:3] <- c(n.sig.AvB, n.sig.A_LMH, n.sig.B_LMH)
  }   # end k 
  return(list(gandy.out = gandy.out, sig.tests.out = sig.tests.out))
}


## Some run times...
# plan(4, 5), BATCH.SIZE = 200 = 246 s, 4 outputs, 4.1 min, 61.5s/iteration
# plan(5, 4), BATCH.SIZE = 250 = 290 s, 5 outputs, 4.83 min, 58s/iteration
#plan(10, 2), BATCH.SIZE = 500 =  533, 10 outputs.  8.88 min, 53.3 sec/iteration
#plan(20, 1), BATCH.SIZE = 1000 = 1110 , 20 outputs.  18.5 min, 55.5 sec/iteration

#plan(2, 10), BATCH.SIZE = 500 =  613s, 10 outputs.  10.2 min, 61.3 sec/iteration
#plan(5, 4) BATCH.SIZE = 200 = 305 s, 5 outputs, 5.08 min, 61 s/iteration

# plan(10, 2), batch size 500 = 1048.78 s, 20 outputs, 17.48 min, 52.45 sec/iteration. 14.6 hrs on my machine for 1000 outputs




## With 76 cores...
#  I recommend trying BATCH.SIZE = 500, n.output = 38, n.process = 2:
n.output <- 38   # Number of output per cycle
n.process <- 76/n.output
BATCH.SIZE <- 1000/n.process
# I estimate your machine will complete the task in about 4 hours.


total.output <- 1000  # Maybe try with 10 first?  There will be some warnings about unused cores, I think.

#n.output <- 4
#n.process <- 20/n.output
#BATCH.SIZE <- 1000/n.process

# Here we set up worker cores for multi-level parallel processing.  
plan(list(
  tweak(multisession, workers =  n.output),  # Outer level (outputs)
  tweak(multisession, workers = I(n.process)) # Inner level (iterations)
))

#print(plan('list'))

tic()
x <- furrr::future_map(1:total.output, function(i) {
  J <-i 
  y <- perm.test.fcn(arrays.grow = arrays.grow,
                     arrays.mort = arrays.mort, 
                     PlotDat = PlotDat)
  return(y)
  
}, .options = furrr_options(seed = TRUE))
toc()


plan(sequential)  # This closes parallel workers



## Code for opening previously-run permutation tests
x <- read_rds(paste0(DATA.LOC, "perm_test_out_092625.RDS"))

gandy_table <- data.frame(do.call(rbind, lapply(x, function(xi) xi$gandy.out)))

n.sto <- length(x)
sto.rows <- nrow(x[[1]]$sig.tests.out)
sto.cols <- ncol(x[[1]]$sig.tests.out)

sig_tests_table <- map_dfr(seq_along(x), ~ {
  as.data.frame(x[[.x]]$sig.tests.out) %>%
    mutate(id = .x)
})


sig_array <- array(dim = c(sto.rows, sto.cols, n.sto))

# Fill each slice with the corresponding matrix
for(i in 1:n.sto) {
  sig_array[, , i] <- x[[i]]$sig.tests.out
}


odd_cols <- seq(3, 13, by = 2)  # Creates c(3, 5, 7, 9, 11, 13)

sig_array2 <- sig_array[, odd_cols, ]
  
sig_array2 <- apply(sig_array2, c(1, 2, 3), as.numeric)

sig_array_prop <- apply(sig_array2, c(1, 2), function(x) mean(x == 1, na.rm = TRUE))

# Convert to data frame
proportions_df <- as.data.frame(proportions)
colnames(proportions_df) <- paste0("col_", odd_cols)


#tic()
#perm_results <- 1:PERM.ITER.N %>% 
#  future_map(\(perm_n) perm.test.fcn(arrays.grow = arrays.grow,
#                                     arrays.mort = arrays.mort, 
#                                     PlotDat = PlotDat
#                                     ), .options = furrr_options(seed = TRUE))
#toc()


# The Gandy simctest needs a function that makes a comparison.  In this case the 
# comparison is to see if each iteration produces a TRUE value.
#gen <- function(){perm.test.fcn()[contrast.min.position] >= contrast.min.val}

# Setting up parallel computing. See Global.R for n.cores.
#plan(multisession, workers = n.cores) 


#tic()
#simctest(gen, maxsteps = PERM.ITER.N)  # Here's the permutation test.
#toc()



