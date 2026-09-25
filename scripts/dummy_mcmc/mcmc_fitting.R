## MCMC FITTING - UNKNOWN REPORTING RATES ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
suppressMessages(require(BayesianTools))
suppressMessages(require(parallel))
options(dplyr.summarise.inform = FALSE) 

.args <- #if (interactive()) c(
  c(file.path("data", "inputs", "imd_age_pop.rds"),
    file.path("data", "inputs", "contact_matrix.rds"),
    file.path("data", "inputs", "subtype_years.rds"),
    file.path("data", "dummy_data", "dummy_surveillance.rds"),
    file.path("data", "dummy_data", "known_parameters.rds"),
    file.path("output", "data", "mcmc_samples.rds")
  ) #else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

set.seed(60)

#### LOAD DATA ####

## read in population data
imd_age_pop_reg <- readRDS(.args[1])
age_labels <- unique(imd_age_pop_reg$age_grp)
nage <- n_distinct(imd_age_pop_reg$age_grp)
nimd <- n_distinct(imd_age_pop_reg$imd_quintile)

## aggregate
imd_age_pop <- imd_age_pop_reg %>% 
  group_by(age_grp, imd_quintile) %>% 
  summarise(pop = sum(pop))

imd_age_pop$age_grp <- factor(imd_age_pop$age_grp, levels = age_labels)

## ENSURE ALL POP ORDERS ARE IMD THEN AGE
imd_age_pop <- imd_age_pop %>% 
  arrange(imd_quintile, age_grp)

## read in contact matrix
contact_matrix_1000 <- readRDS(.args[2])

## aggregate
contact_matrix <- contact_matrix_1000 %>% 
  group_by(p_imd_q, c_imd_q, p_age_group, c_age_group) %>% 
  summarise(n = mean(n))

contact_matrix$p_age_group <- factor(contact_matrix$p_age_group,
                                     levels = age_labels)
contact_matrix$c_age_group <- factor(contact_matrix$c_age_group,
                                     levels = age_labels)

## ENSURE ALL POP ORDERS ARE IMD THEN AGE
contact_matrix <- contact_matrix %>% 
  arrange(p_imd_q, c_imd_q,
          p_age_group, c_age_group,)

## into matrix
cm <- contact_matrix %>% ungroup() %>% 
  mutate(p_var = paste0(p_imd_q, '_', p_age_group),
         c_var = paste0(c_imd_q, '_', c_age_group)) %>% 
  select(p_var,c_var,n) %>% 
  pivot_wider(names_from = c_var, values_from = n) %>% 
  select(!p_var) %>% as.matrix()

## make matrix per capita
per_cap_matrix_45 <- t(t(cm)/imd_age_pop$pop)

## scale up for both risk groups 
## (assuming per capita contacts are independent of risk level)
pc_cm <- expand_contact_matrix(per_cap_matrix_45)

ng <- nrow(per_cap_matrix_45)
ndim <- nrow(pc_cm)
if(!all.equal(2*ng, ndim)){warning('dimensions not adding up')}
if(!all.equal(2*nimd*nage, ndim)){warning('dimensions not adding up')}

## SUBTYPE-SEASONS
subtype_seasons <- readRDS(.args[3])

## SURVEILLANCE DATA
surveillance_data <- readRDS(.args[4])

## KNOWN PARAMETERS
known_pars <- readRDS(.args[5])

years <- known_pars$years

delays <- c(known_pars$primary_care_delay, known_pars$secondary_care_delay)
names(delays) <- c('primary','secondary')

#### RUNNING MCMC ####

## MCMC pars
nchains <- 3
n_pop = 100
burn_in <- 200
thinning_value <- 1
n_samples <- 2000

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))  
if (is.na(n_cores) || n_cores < 1) n_cores <- 1            # safe fallback if run outside SLURM

mcmc_parallel <- function(i){
  
  txt_output <<- i
  year_i <- years[i]
  season_i <- paste0(year_i, '/', substr(year_i + 1, 3, 4))
  
  ## population data
  
  vaccinated_data_seasonal <- known_pars$vaccinated_data
  vaccinated_data_seasonal$age_grp <- factor(vaccinated_data_seasonal$age_grp, levels = age_labels)
  vaccinated_data_seasonal <- vaccinated_data_seasonal %>% 
    filter(start_of_season == year_i) %>% 
    arrange(subtype, desc(risk_level), imd_quintile, age_grp)
  vaccinated_data_seasonal_no_rep <- vaccinated_data_seasonal %>% 
    filter(subtype == vaccinated_data_seasonal$subtype[1])
  
  ## should be ordered by IMD then age
  if(vaccinated_data_seasonal$imd_quintile[2] != 1){warning('vaccinated_data_seasonal in wrong order')}
  if(vaccinated_data_seasonal$age_grp[2] != age_labels[2]){warning('vaccinated_data_seasonal in wrong order')}
  
  demography <- vaccinated_data_seasonal_no_rep %>% 
    mutate(population = pop) %>% 
    select(age_grp, imd_quintile, risk_level, population, risk_proportion) %>% 
    arrange(desc(risk_level), imd_quintile, age_grp)
  
  ## check population sum is correct
  tot_pop <- sum(imd_age_pop$pop)
  if(!all.equal(sum(demography$population), tot_pop)){warning('pop not adding up')}
  
  subtype_init_pars <- c(0.2, rep(0.5, 3), 2.5, 
                         rep(0.02, 6), rep(0.002, 6))
                         # c(transmissibility, 3x absolute susceptibility, log of initial infected, 
                         #   reporting rates for primary care, reporting rates for secondary care)
  
  run_mcmc_inference(
    demography_input = demography, 
    vaccinated_input = vaccinated_data_seasonal,
    subtype_season = subtype_seasons %>% filter(season == season_i),
    cm_input = pc_cm,
    epidemic_to_fit = surveillance_data %>% filter(index == i),
    epid_periods = known_pars$epid_periods,
    coverage_rates = known_pars$proportion_observed,
    care_delays = delays,
    initial_parameters = c(subtype_init_pars,
                           subtype_init_pars,
                           rep(0, 4)),
    #   subtype-specific parameters x2, IMD spline parameters x4
    n_samples = n_samples*nchains, 
    nburn = burn_in*nchains, 
    thinning = thinning_value,
    n_pop = n_pop,
    n_cores = n_cores,
    n_chains = 1 # the DEzs sampler produces three subchains, dealt with by
    # multiplying nburn and n_samples by 3
  )
}

mcmc_results <- mclapply(1:length(years), mcmc_parallel, mc.cores = length(years))

#### SAVE RESULTS ####

# save most recently run settings as a dummy save
write_rds(data.table(x=paste0(burn_in,'_',thinning_value,'_',n_samples),
                     HPC = F,
                     date = Sys.Date()), .args[6]) 

# save actual data
write_rds(mcmc_results, gsub('.rds',paste0('_', burn_in,'_',thinning_value,'_',n_samples,'_',Sys.Date(),'.rds'),
                             .args[6]))

