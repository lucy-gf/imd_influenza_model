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
    file.path("data", "dummy_data", "dummy_surveillance.rds"),
    file.path("data", "dummy_data", "known_parameters.rds"),
    file.path("output", "data", "mcmc_samples_rates_unknown.rds")
  ) #else commandArgs(trailingOnly = TRUE)

i <- as.numeric(commandArgs(trailingOnly = TRUE))

source(file.path('scripts','setup','colors.R'))
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

## SURVEILLANCE DATA
surveillance_data <- readRDS(.args[3])

## KNOWN PARAMETERS
known_pars <- readRDS(.args[4])

years <- known_pars$years

delays <- c(known_pars$primary_care_delay, known_pars$secondary_care_delay)
names(delays) <- c('primary','secondary')

risk_group_pop <- known_pars$risk_group_pop 
risk_group_pop$age_grp <- factor(risk_group_pop$age_grp, levels = age_labels)
risk_group_pop <- risk_group_pop %>% 
  arrange(imd_quintile, age_grp)

vaccinated_pop <- known_pars$vaccinated_pop
vaccinated_pop$age_grp <- factor(vaccinated_pop$age_grp, levels = age_labels)
vaccinated_pop <- vaccinated_pop %>% 
  left_join(known_pars$vaccination_efficacy, by = 'age_grp') %>% 
  mutate(effectively_vaccinated_population = VE*vaccinated_population) %>% 
  arrange(desc(risk_level), imd_quintile, age_grp)

demography <- rbind(risk_group_pop %>% mutate(risk_level = 'high'),
                    risk_group_pop %>% mutate(risk_level = 'low')) %>% 
  mutate(population = case_when(risk_level == 'high' ~ risk_population,
                                risk_level == 'low' ~ pop - risk_population)) %>% 
  select(!c(risk_population,pop)) %>% arrange(desc(risk_level), imd_quintile, age_grp)

## check population sum is correct
tot_pop <- sum(imd_age_pop$pop)
if(!all.equal(sum(demography$population), tot_pop)){warning('pop not adding up')}

#### RUNNING MCMC ####

nchains <- 3
burn_in <- 3000
thinning_value <- 1
n_samples <- 3000

mcmc_results <- run_mcmc_inference(
  demography_input = demography, 
  vaccinated_input = vaccinated_pop,
  cm_input = pc_cm, 
  epidemic_to_fit = surveillance_data %>% filter(index==i), 
  epid_periods = known_pars$epid_periods,
  coverage_rates = known_pars$proportion_observed,
  care_delays = delays,
  initial_parameters = c(0.07, rep(1, 2), 2, 
                         rep(0.02, 6), rep(0.002, 6),
                         rep(0, 4)),
  # c(transmissibility, 2x relative susceptibility, log of initial infected, 
  #   reporting rates for primary care, reporting rates for secondary care,
  #   IMD spline parameters x4)
  n_samples = n_samples*nchains, 
  nburn = burn_in*nchains, 
  thinning = thinning_value,
  n_chains = 1, # the DEzs sampler produces three subchains, dealt with by
  # multiplying nburn and n_samples by 3
  txt_output = i
)

write_rds(mcmc_results, gsub('.rds',paste0('_', i, '_', burn_in,'_',thinning_value,'_',n_samples,'.rds'),
                             .args[5])) # in case next save fails

mcmc_results <- c(mcmc_results,
                  list(burn_in=burn_in, thinning_value=thinning_value, n_samples=n_samples))

#### SAVE ####

write_rds(data.table(x=1), .args[5]) # dummy save
write_rds(mcmc_results, gsub('.rds',paste0('_', i, '_', burn_in,'_',thinning_value,'_',n_samples,'.rds'),
                             .args[5]))

