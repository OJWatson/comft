library(magenta)
library(tidyverse)

## ----------------------------------------------------o
## 1. Setting up a cluster configuration --------------
## ----------------------------------------------------o

# Setting Up Cluster From New

# Log in to didehpc
credentials = "C:/Users/ow813/.smbcredentials"
options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "ow813")

# not if T is not mapped then map network drive
didehpc::didehpc_config_global(temp=didehpc::path_mapping("tmp",
                                                          "T:",
                                                          "//fi--didef3.dide.ic.ac.uk/tmp",
                                                          "T:"),
                               home=didehpc::path_mapping("OJ",
                                                          "Z:",
                                                          "//wpia-hn/Malaria/nas5_malaria",
                                                          "Z:"),
                               credentials=credentials,
                               cluster = "wpia-hn")
#cluster = "fi--didemrchnb")

# Creating a Context
context_name <- "analysis/context-dide"

ctx <- context::context_save(
  path = context_name,
  package_sources = conan::conan_sources(
    packages= c(
      "local::scripts/binaries/magenta_1.3.1.zip",
      "local::scripts/binaries/dde_1.0.2.zip",
      "rrq"),
    repos = "https://mrc-ide.github.io/drat/"
  )
)

# set up a specific config for here as we need to specify the large RAM nodes
config <- didehpc::didehpc_config(use_workers = TRUE, parallel = FALSE, cores = 3)

# Configure the Queue
obj <- didehpc::queue_didehpc(ctx, config = config)


## ----------------------------------------------------o
## 2. Create cluster param grid --------------
## ----------------------------------------------------o

#### Parameter Grid ####

# parameter set up
N <- 100000
tl <- 20
hd <- 20
nl <- 7
m <- 0.00084

# Default starting resistance
crt_76T <- 0.1
k13_valid <- 0.1
mdr1_184F <- 0.2
mdr1_86Y <- 0.02
mdr1_CNV <- 0.05
pfpm23_CNV <- 0.02

# create our drugs
aspy <- magenta:::drug_create_dhappq()
dhappq <- magenta:::drug_create_dhappq()
asaq <- magenta:::drug_create_asaq()
al <- magenta:::drug_create_al()

# give it the same properties as AL but genetically spread like dhappq (i.e. single locus pd)
aspy$dur_P <- al$dur_P
aspy$dur_SPC <- al$dur_SPC
aspy$hill_n <- al$hill_n
aspy$hill_kA <- al$hill_kA
aspy$hill_res_n <- al$hill_res_n
aspy$hill_res_kA <- al$hill_res_kA

# but give it lpc properties that are similar to asaq but genetically spread like dhappq (i.e. single locus pd
# that is not yet known.
aspy$lpf[1:16] <- max(al$lpf)
aspy$lpf[17:32] <- 0.895737
aspy$lpf[33:48] <- max(al$lpf)
aspy$lpf[49:64] <- 0.895737
aspy$lpf[65:80] <- 0.891124
aspy$lpf[81:96] <- 0.735268
aspy$lpf[97:112] <- 0.891124
aspy$lpf[113:128] <- 0.735268
aspy$prophylactic_positions <- 6

# update positions for 7 loci
aspy$barcode_positions <- 0:6
dhappq$barcode_positions <- 0:6
asaq$barcode_positions <- 0:6
al$barcode_positions <- 0:6

# and extend lpf
dhappq$lpf <- c(dhappq$lpf, dhappq$lpf)
asaq$lpf <- c(asaq$lpf, asaq$lpf)
al$lpf <- c(al$lpf, al$lpf)

# all drugs
dl <- list(al, asaq, dhappq, aspy)

## Now create what we scan over
EIRs <- c(0.5,1,2,4,8,16,32,64,128)
fts <- c(.2,.4,.6,.8)
num_drugs <- c(2, 3, 4)

# then specific scenario scans
cycle_delays <- c(1, 2, 4)
mft_prop <- c("even", "al_heavy")
p_tests <- c(.2,.4,.6,.8,.1)
al_heavy_ratios <- list(
  c(0.8, 0.2),
  c(4/6, 1/6, 1/6),
  c(4/7, 1/7, 1/7, 1/7)
)

# create our status quo grid
status_quo_grid <- expand.grid(
  EIR = EIRs, 
  ft = fts, 
  p_test = p_tests,
  nd = 1,
  cycle_delays = -1
)
status_quo_grid$sequential <- -1
status_quo_grid$mft_flag <- FALSE
status_quo_grid$mft_prop <- "even"

# create our temporal cycling grid
temporal_cycle_grid <- expand.grid(
  EIR = EIRs, 
  ft = fts, 
  p_test = p_tests,
  nd = num_drugs,
  cycle_delays = cycle_delays
)
temporal_cycle_grid$sequential <- -1
temporal_cycle_grid$mft_flag <- FALSE
temporal_cycle_grid$mft_prop <- "even"

# create our sequential cycling grid
sequential_cycle_grid <- expand.grid(
  EIR = EIRs, 
  ft = fts, 
  p_test = p_tests,
  nd = num_drugs,
  cycle_delays = cycle_delays
)
sequential_cycle_grid$sequential <- 0.15
sequential_cycle_grid$mft_flag <- FALSE
sequential_cycle_grid$mft_prop <- "even"

# create our mft grid
mft_grid <- expand.grid(
  EIR = EIRs, 
  ft = fts, 
  p_test = p_tests,
  nd = num_drugs,
  mft_prop = mft_prop
)
mft_grid$sequential <- -1
mft_grid$mft_flag <- TRUE
mft_grid$cycle_delays <- -1

# combine into one
param_grid <- data.table::rbindlist(
  list(status_quo_grid, sequential_cycle_grid, temporal_cycle_grid, mft_grid),
                                    use.names = TRUE)

param_grid <- arrange(param_grid, EIR)

## ----------------------------------------------------o
## 3. Make cluster submissions --------------
## ----------------------------------------------------

## safe submission
try_fail_catch <- function(expr, attempts = 3){
  r <- NULL
  attempt <- 1
  while( is.null(r) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      r <- eval(expr)
    )
  }
  
}

for(i in seq_along(param_grid$EIR)){
  
  # get partner drug ratios
  if (param_grid$mft_prop[i] == "even") {
    pdr <- rep(1/param_grid$nd[i], param_grid$nd[i])
  } else {
    pdr <- al_heavy_ratios[[param_grid$nd[i]-1]]
  }
  
  drug_list <-
    magenta:::drug_list_create(
      resistance_flag = c(rep(FALSE, hd), rep(TRUE, tl)),
      mft_flag = param_grid$mft_flag[i],
      artemisinin_loci = 5,
      absolute_fitness_cost_flag = TRUE,
      partner_drug_ratios = pdr,
      drugs = dl[seq_len(param_grid$nd[i])],
      cost_of_resistance = rep(1 - 0.005, 6),
      sequential_cycling = param_grid$sequential[i],
      number_of_drugs = param_grid$nd[i],
      sequential_update = param_grid$cycle_delays[i],
      number_of_resistance_loci = nl
    )
  
  param_list_lmh <- list()
  
  for(j in 1:2){
    
    param_list_lmh[[j]] <- list(
      EIR = param_grid$EIR[i], 
      N = 20000, 
      years = tl+hd,
      itn_cov = 0,
      save_lineages = TRUE,
      irs_cov = 0,
      ft = param_grid$ft[i],
      num_loci = nl, 
      sample_reps = 1,
      mutation_flag=c(rep(FALSE,hd),rep(TRUE,tl)),
      mutation_treated_modifier = 100,
      mutation_rate=c(m,m,m,m,m,3*m,m*0.5),
      survival_percentage = 0.22,
      genetics_df_without_summarising = TRUE,
      spatial_incidence_matrix = c(rep(1,hd),rep(0,tl)),
      spatial_mosquitoFOI_matrix = c(rep(1,hd),rep(0,tl)),
      human_only_full_save=FALSE,
      spatial_type = "island",
      use_historic_interventions = TRUE, 
      seed = as.integer(runif(1, 1, 1000000000)),
      sample_size = c(100,1000),
      sample_states = c(1,2,4),
      genetics_df_without_summarising = TRUE,
      ibd_length = 1, update_length = 30,
      plaf=matrix(c(rep(c(crt_76T,mdr1_86Y,mdr1_184F,mdr1_CNV,k13_valid,pfpm23_CNV, 0),hd),rep(0,nl*tl)),ncol=nl,byrow=TRUE),
      update_save = TRUE, 
      human_update_save = TRUE,
      summary_saves_only = TRUE,
      housekeeping_list = magenta:::housekeeping_list_create(quiet = TRUE,cluster = TRUE),
      drug_list = drug_list,
      nmf_list = magenta:::nmf_list_create(nmf_flag = TRUE, prob_of_testing_nmf = param_grid$p_test[i]),
      island_imports_plaf_linked_flag = FALSE
    )
  }
  
  
  try_fail_catch(grp <- obj$lapply(
    X = param_list_lmh, 
    timeout=0, 
    FUN = function(x){
      return(magenta::pipeline(EIR=x$EIR,
                               seed=x$seed,
                               save_lineages = x$save_lineages,
                               N=x$N,
                               mutation_rate = x$mutation_rate, 
                               mutation_flag = x$mutation_flag,
                               sample_size = x$sample_size,
                               sample_reps = x$sample_reps,
                               years=x$years,
                               survival_percentage=x$survival_percentage,
                               itn_cov=x$itn_cov,
                               irs_cov=x$irs_cov,
                               ft=x$ft,
                               genetics_df_without_summarising = x$genetics_df_without_summarising,
                               spatial_incidence_matrix=x$spatial_incidence_matrix,
                               spatial_mosquitoFOI_matrix=x$spatial_mosquitoFOI_matrix,
                               spatial_type=x$spatial_type,
                               use_historic_interventions=x$use_historic_interventions,
                               human_only_full_save=x$human_only_full_save,
                               ibd_length=x$ibd_length,
                               num_loci=x$num_loci,
                               sample_states = x$sample_states,
                               update_length=x$update_length,
                               update_save=x$update_save,
                               human_update_save=x$human_update_save,
                               summary_saves_only=x$summary_saves_only,
                               housekeeping_list=x$housekeeping_list,
                               drug_list = x$drug_list, 
                               nmf_list = x$nmf_list, plaf = x$plaf,
                               island_imports_plaf_linked_flag = x$island_imports_plaf_linked_flag),
             name = paste0("test_pg_",i))
    }
    
  ))
  
}

# now submit workers
workers <- obj$submit_workers(200)

# and get the grp list here
grps <- lapply(obj$task_bundle_list(), function(x){obj$task_bundle_get(x)})

## ----------------------------------------------------o
## 4. Fetch Objects --------------
## ----------------------------------------------------o

# let's put our sim outputs here
dir.create(here::here("analysis/data-derived/sims_test"))

# first save the parms used to generate (X)
saveRDS(param_grid, here::here("analysis/data-derived/sims_test/X.rds"))

# Save our simulations across the parameter space
for(i in seq_along(param_grid$EIR)) {
  
  message(i)
  
  # create our results
  r_i <- map(seq_along(grps), function(x) {
    grps[[x]]$db$get_value(grps[[x]]$db$get_hash(grps[[x]]$tasks[[i]]$id, "task_results"), FALSE)
  })
  
  if(any(unlist(lapply(r_i, object.size))==0)){
    
  } else {
    
    # remove the final loggers that are very large and unneeded
    for(x in seq_along(r_i)) {
      to_rm <- which(names(r_i[[x]][[length(r_i[[x]])]]$Loggers) %in%
                       c("InfectionStates", "Ages", "IB", "ICA", "ICM", "ID") )
      r_i[[x]][[length(r_i[[x]])]]$Loggers[to_rm] <- NULL
    }
    
    # and save to file
    fn_i <- paste0("r_", i, ".rds")
    saveRDS(r_i, here::here("analysis/data-derived/sims_test", fn_i))
    gc()
    
  }
  
}

# These sims have not been pushed to Github due to Github memory/size constraints

# Analysis in next script uses these sims and the output of this is then saved
# in data-derived

