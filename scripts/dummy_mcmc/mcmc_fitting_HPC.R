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

.args <- c(file.path("data", "inputs", "imd_age_pop.rds"),
           file.path("data", "inputs", "contact_matrix.rds"),
           file.path("data", "dummy_data", "dummy_surveillance.rds"),
           file.path("data", "dummy_data", "known_parameters.rds"),
           file.path("data", "dummy_data", "mcmc_true_values.rds"),
           file.path("output", "data", "mcmc_samples.rds"))

bash_args <- if (interactive()) c(3, 1) else as.numeric(commandArgs(trailingOnly = TRUE))

#### RUN OPTIONS ####

i <- bash_args[1]
chain <- bash_args[2]
set.seed(chain)
message(paste0('Index: ', i, ', Chain: ', chain))

## Should high-risk IGPR and IHR be fixed in relation to low-risk?
RISK_RATIOS_FIXED <- T

## Will this run be merging with a previous run?
PREVIOUS_PARS <- F

###########################################################################################

## SOURCE SETUP SCRIPTS
source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

###########################################################################################

#### LOAD DATA ####

{
## read in population data
imd_age_pop_reg <- readRDS(.args[1])
age_labels <- unique(imd_age_pop_reg$age_grp)
nage <- n_distinct(imd_age_pop_reg$age_grp)
nimd <- n_distinct(imd_age_pop_reg$imd_quintile)

broad_ages <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages('children','adults','older_adults',age_labels)
)

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
}

## SURVEILLANCE DATA
surveillance_data <- readRDS(.args[3])[[1]]
subtype_seasons <- readRDS(.args[3])[[2]]

## KNOWN PARAMETERS
known_pars <- readRDS(.args[4])
years <- known_pars$years
delays <- c(known_pars$primary_care_delay, known_pars$secondary_care_delay)
names(delays) <- c('primary','secondary')
year_i <- years[i]
season_i <- paste0(year_i, '/', substr(year_i + 1, 3, 4))

n_subtypes <- n_distinct((subtype_seasons %>% filter(season==season_i))$subtype)
TWO_SUBTYPES <- (n_subtypes == 2)

## TRUE MCMC VALUES
mcmc_true_values <- readRDS(.args[5]) %>% 
  filter(season == season_i, name != 'Reff')

if(RISK_RATIOS_FIXED){
  mcmc_true_values <- mcmc_true_values %>% filter(!grepl('_high', name))
}

## POPULATION DATA
{
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
}

# crude estimate of IGPR and IHR
est_attack_rate <- 0.15
est_health_ratios <- surveillance_data %>% filter(index == i) %>%
  left_join(known_pars$proportion_observed,
            by = c('age_grp','imd_quintile','risk_level')) %>% 
  ungroup() %>% 
  summarise(primary = sum(primary_care/OS_COVERAGE),
            secondary = sum(secondary_care/OS_COVERAGE)) %>% 
  mutate(pop = sum(demography$population)) %>% 
  summarise(prim = sum(primary/(pop*est_attack_rate)),
            sec = sum(secondary/(pop*est_attack_rate)))

###########################################################################################

#### PARAMETERS TO FIT ####

# parameter names
subtype_parameter_names <- c('transmissibility',
                             'children_susceptibility','adults_susceptibility','older_adults_susceptibility',
                             'init_infected',
                             'primary_care_rate_children_low','primary_care_rate_adults_low','primary_care_rate_older_adults_low',
                             'primary_care_rate_children_high','primary_care_rate_adults_high','primary_care_rate_older_adults_high',
                             'secondary_care_rate_children_low','secondary_care_rate_adults_low','secondary_care_rate_older_adults_low',
                             'secondary_care_rate_children_high','secondary_care_rate_adults_high','secondary_care_rate_older_adults_high')
imd_spline_names <- c('imd_spline_primary_1','imd_spline_primary_2','imd_spline_secondary_1','imd_spline_secondary_2')

# if risk ratios are fixed, only 3 parameters to fit for each
N_IGPR_VARS <- N_IHR_VARS <- ifelse(RISK_RATIOS_FIXED, 3, 6)
if(RISK_RATIOS_FIXED){
  subtype_parameter_names <- subtype_parameter_names[!grepl('_high', subtype_parameter_names)]
}

# not necessarily all parameters have to be fitted in the MCMC
# (will be fixed at true values if not being fitted)

trans_fit <- F # transmissibility
susc_fit <- T # susceptibility (x3)
init_fit <- F # initial infected
primary_fit <- T # IGPRs (x3)
secondary_fit <- T # IHRs (x3)
imd_spline_fit <- T # IMD splines (x4)

message(paste0('Parameters fitted: ',
               ifelse(trans_fit,'transmissibility, ', ''),
               ifelse(susc_fit,'susceptibility, ', ''),
               ifelse(init_fit,'initial infected, ', ''),
               ifelse(primary_fit,'primary care rates, ', ''),
               ifelse(secondary_fit,'secondary care rates, ', ''),
               ifelse(imd_spline_fit,'IMD splines, ', ''),
               '\nParameters fixed:',
               ifelse(!trans_fit,', transmissibility', ''),
               ifelse(!susc_fit,', susceptibility', ''),
               ifelse(!init_fit,', initial infected', ''),
               ifelse(!primary_fit,', primary care rates', ''),
               ifelse(!secondary_fit,', secondary care rates', ''),
               ifelse(!imd_spline_fit,', IMD splines', '')))

subtype_pars_to_fit <- c(trans_fit,
                         rep(susc_fit, 3),
                         init_fit,
                         rep(primary_fit, N_IGPR_VARS),
                         rep(secondary_fit, N_IHR_VARS))

pars_to_fit <- subtype_pars_to_fit
names_pars_to_fit <- paste0(subtype_parameter_names, '_epid_1')

if(TWO_SUBTYPES){
  pars_to_fit <- c(pars_to_fit, subtype_pars_to_fit)
  names_pars_to_fit <- c(names_pars_to_fit,
                         paste0(subtype_parameter_names, '_epid_2'))
}

pars_to_fit <- c(pars_to_fit, rep(imd_spline_fit, 4))
names_pars_to_fit <- c(names_pars_to_fit, imd_spline_names)
names(pars_to_fit) <- names_pars_to_fit

## set true value vector
true_values_in <- mcmc_true_values$value 
if(length(pars_to_fit) != length(true_values_in)){stop('True value vector wrong length!')}
true_values_in[pars_to_fit] <- NA
names(true_values_in) <- names_pars_to_fit

###########################################################################################

#### INITIAL VALUES ####

# c(transmissibility, 3x absolute susceptibility, log of initial infected, 
#   reporting rates for primary care, reporting rates for secondary care)
subtype_init_pars <- c(0.2, rep(0.5, 3), 2.5, 
                       rep(est_health_ratios$prim[1], N_IGPR_VARS), rep(est_health_ratios$sec[1], N_IHR_VARS))

# filter to only those being fitted in this analysis 
subtype_init_pars <- subtype_init_pars[subtype_pars_to_fit]

# if in 2025/26, only using one set of epidemiological parameters (else using 2)
# (subtype-specific parameters x2, IMD spline parameters x4)

initial_parameters_vec <- subtype_init_pars
if(TWO_SUBTYPES){initial_parameters_vec <- c(initial_parameters_vec, subtype_init_pars)}
if(imd_spline_fit){initial_parameters_vec <- c(initial_parameters_vec, rep(0,4))}

initial_parameters_in <- matrix(initial_parameters_vec, nrow=1)

# check right number of columns
if(ncol(initial_parameters_in) != sum(pars_to_fit)){stop('initial_parameters_in wrong column count')}

###########################################################################################

#### SETUP MCMC ####

## MCMC pars
burn_in <- 0
thinning_value <- 1
n_samples <- 200

number_string <- paste0(burn_in,'_',thinning_value,'_',n_samples)

# make output folder
if(!file.exists(paste0('output/data/',number_string))){
  dir.create(paste0('output/data/',number_string))
}

cat('\n Total steps: ', burn_in + thinning_value*n_samples, ' (', 
    number_string, ')\n', sep = '')

# save all steps, or only those after burn in and thinning?
save_all_flag <- F

## OLD CODE FOR USING N_CORES > 1
# n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))  
# if (is.na(n_cores) || n_cores < 1) n_cores <- 1 # safe fallback if run outside SLURM
# cat('\nn_cores: ', n_cores, '\n')

if(PREVIOUS_PARS){
  # load in previous parameters 
  previous_run <- readRDS(.args[6])
  previous_par_files <- list.files(paste0("output/data/",previous_run$x))
  previous_par_files_filt <- previous_par_files[grepl(paste0('INDEX',i,'_CHAIN',chain,'_'),previous_par_files)]
  previous_pars <- readRDS(file.path('output','data',previous_run$x,previous_par_files_filt))
  
  parameters_vec <- sum(pars_to_fit) # TODO Fix when merging with not all parameters being fitted
  
  # make the end of the last run into the initial inputs for this run
  initial_parameters_in <- matrix(
    c(previous_pars[[1]][nrow(previous_pars[[1]]),1:parameters_vec],
      previous_pars[[2]][nrow(previous_pars[[2]]),1:parameters_vec],
      previous_pars[[3]][nrow(previous_pars[[3]]),1:parameters_vec]),
    nrow = parameters_vec)
  
  initial_parameters_in <- t(initial_parameters_in)
  
}

###########################################################################################

#### RUN MCMC ####

mcmc_results <- run_mcmc_inference(
  demography_input = demography,
  vaccinated_input = vaccinated_data_seasonal,
  subtype_season = subtype_seasons %>% filter(season == season_i),
  cm_input = pc_cm,
  epidemic_to_fit = surveillance_data %>% filter(index == i),
  epid_periods = known_pars$epid_periods,
  coverage_rates = known_pars$proportion_observed,
  care_delays = delays,
  initial_parameters = initial_parameters_in,
  save_all = save_all_flag,
  parameters_to_fit = pars_to_fit,
  true_values = true_values_in,
  n_samples = n_samples,
  nburn = burn_in,
  thinning = thinning_value,
  n_chains = 1,
  txt_output = paste0(i, '_', chain),
  manual_print_mod = 10
)

#### SAVE ####

# save most recently run settings as a dummy save
write_rds(list(x=number_string,
               HPC = T,
               save_all = save_all_flag,
               date = Sys.Date(),
               pars_to_fit = pars_to_fit,
               using_previous_pars = PREVIOUS_PARS,
               previous_pars = if(PREVIOUS_PARS){previous_run$x}else{NULL},
               previous_pars_date = if(PREVIOUS_PARS){previous_run$date}else{NULL},
               previous_pars_save_all = if(PREVIOUS_PARS){previous_run$save_all}else{NULL}), 
          .args[6])

# ALSO save as a differently named file in case multiple runs happening
write_rds(list(x=number_string,
               HPC = T,
               save_all = save_all_flag,
               date = Sys.Date(),
               pars_to_fit = pars_to_fit,
               using_previous_pars = PREVIOUS_PARS,
               previous_pars = if(PREVIOUS_PARS){previous_run$x}else{NULL},
               previous_pars_date = if(PREVIOUS_PARS){previous_run$date}else{NULL},
               previous_pars_save_all = if(PREVIOUS_PARS){previous_run$save_all}else{NULL}), 
          gsub('.rds',paste0('_', number_string,'_',paste(pars_to_fit, collapse = '_'),'.rds'),
               .args[6]))

# save actual data
write_rds(mcmc_results, gsub('.rds',paste0('_INDEX', i, '_CHAIN', chain, '_', 
                                           number_string,'_',Sys.Date(),'.rds'),
                             gsub('data/',paste0('data/',number_string,'/'),.args[6])))


