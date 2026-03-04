## CREATE DUMMY DATA (INFECTIONS) ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "contact_matrix.rds"),
  file.path("data", "inputs", "dummy_flu_data.rds")
) else commandArgs(trailingOnly = TRUE)
  
source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','dummy_seir_model.R'))

set.seed(60)

#### LOAD DATA ####

## number of years of data
years <- 2023:2025 # 2023-24 to 2025-26

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

#### SUPPLEMENT WITH DUMMY INPUTS ####

## vaccination coverage

vaccination_coverage_kids <- 0.5
vaccination_coverage_elderly <- 0.7
vaccination_coverage_risk <- 0.5

## risk group proportions

risk_group_kids <- 0.1
risk_group_adults <- 0.1
risk_group_elderly <- 0.4

#### PROPORTIONS IN RISK GROUPS ####

risk_group_age_vector <- c(
  rep(risk_group_kids, 3),
  rep(risk_group_adults, 4),
  rep(risk_group_elderly, 2)
)

## add some variation by IMD
risk_group_imd_variation <- seq(1.1, 0.9, by = -0.05)

risk_group_pop <- data.frame(rgav = rep(risk_group_age_vector, nimd), 
                             rgiv = rep(risk_group_imd_variation, each = nage)) %>% 
  mutate(age_grp = rep(age_labels, nimd),
         imd_quintile = rep(1:nimd, each = nage),
         risk_proportion = rgav*rgiv) %>% 
  left_join(imd_age_pop, by = c('age_grp','imd_quintile')) %>% 
  mutate(risk_population = risk_proportion*pop) %>% 
  select(!c(rgav, rgiv))

risk_group_pop %>% 
  ggplot() + geom_line(aes(age_grp, risk_population, col=imd_quintile,
                           group=imd_quintile))

risk_group_pop %>% 
  ggplot() + geom_line(aes(age_grp, risk_proportion, col=imd_quintile,
                           group=imd_quintile)) 

#### PROPORTIONS VACCINATED ####

vaccination_coverage_age_vector <- c(
  rep(vaccination_coverage_kids, 3),
  rep(vaccination_coverage_risk, 4),
  rep(vaccination_coverage_elderly, 2)
)

## add some variation by IMD
vaccination_imd_variation <- seq(0.9, 1.1, by = 0.05)

vaccinated_pop_1 <- data.table(vcav = rep(vaccination_coverage_age_vector, nimd), 
                             vciv = rep(vaccination_imd_variation, each = nage)) %>% 
  mutate(age_grp = rep(age_labels, nimd),
         imd_quintile = rep(1:nimd, each = nage),
         vaccinated_proportion = vcav*vciv) %>% 
  left_join(risk_group_pop, by = c('age_grp','imd_quintile')) 

# split by risk groups
vaccinated_pop <- rbind(vaccinated_pop_1 %>% mutate(risk_level = 'low', pop = pop - risk_population),
                        vaccinated_pop_1 %>% mutate(risk_level = 'high', pop = risk_population)) %>% 
  select(!c(risk_population, vcav, vciv))

vaccinated_pop <- vaccinated_pop %>% 
  mutate(vaccinated_population = case_when(
    age_grp %in% c('18-25', '26-34', '35-49', '50-69') & risk_level == 'low' ~ 0,
    T ~ vaccinated_proportion*pop
    )) %>% 
  ## ENSURE ALL POP ORDERS ARE IMD THEN AGE
  arrange(desc(risk_level), imd_quintile, age_grp)

vaccinated_pop %>% group_by(age_grp,imd_quintile) %>% 
  summarise(vaccinated_population=sum(vaccinated_population)) %>% 
  ggplot() + geom_line(aes(age_grp, vaccinated_population, col=imd_quintile,
                           group=imd_quintile))

vaccinated_pop %>% 
  ggplot() + geom_line(aes(age_grp, vaccinated_population, col=imd_quintile,
                           lty = risk_level, group=interaction(imd_quintile,risk_level)))

vaccinated_pop %>% 
  ggplot() + geom_line(aes(age_grp, vaccinated_population/pop, col=imd_quintile,
                           lty = risk_level, group=interaction(imd_quintile,risk_level)))


#### EPIDEMIOLOGICAL PARAMETERS ####

epid_periods <- c(2, 3) # latent and infectious periods

epid_parameters_s1 <- list(
  susceptibility = 0.45, # currently not age-dependent, one for each season
  transmissibility = 0.11,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[1], '35', '1'), "%Y%W%w")
)

epid_parameters_s2 <- list(
  susceptibility = 0.5, # currently not age-dependent, one for each season
  transmissibility = 0.13,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[2], '36', '1'), "%Y%W%w")
)

epid_parameters_s3 <- list(
  susceptibility = 0.55, # currently not age-dependent, one for each season
  transmissibility = 0.08,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[3], '33', '1'), "%Y%W%w")
)

## vaccination efficacy (age-dependent)
VE_pars <- c(0.70, 0.46) # currently just using the NGIV VE estimates with no mismatching

vaccination_efficacy <- c(
  rep(VE_pars[1], 7),
  rep(VE_pars[2], 2)
)

## check R0
# 
# for(i in 1:length(years)){
#   pars <- get(paste0('epid_parameters_s', i))
#   cat('R0 = ', R0_func(pars$susceptibility,
#                 pars$infectious_period,
#                 pars$transmissibility,
#                 cm), '\n', sep = '')
# }

#### EXPAND CONTACT MATRIX ####

expand_contact_matrix <- function(c45) {
  
  ndim <- nrow(c45)
  
  stopifnot(nrow(c45) == ncol(c45))
  
  # Create empty 90x90 matrix
  c90 <- matrix(0, nrow = 2*ndim, ncol = 2*ndim)
  
  # Fill all 4 quadrants with same matrix
  c90[1:ndim, 1:ndim]  <- c45   # low → low
  c90[1:ndim, (ndim+1):(2*ndim)]   <- c45   # low → high
  c90[(ndim+1):(2*ndim), 1:ndim]   <- c45   # high → low
  c90[(ndim+1):(2*ndim), (ndim+1):(2*ndim)] <- c45   # high → high
  
  return(c90)
}

## make matrix per capita
per_cap_matrix_45 <- t(t(cm)/imd_age_pop$pop)

## check calculation assumption is correct
# test_m<-matrix(c(1,2,3, 11,12,12), nrow = 2, ncol = 3, byrow = TRUE)
# should_be_111_1164 <- t(t(test_m)/c(1,2,3))

## scale up for both risk groups 
## (assuming per capita contacts are independent of risk level)
pc_cm <- expand_contact_matrix(per_cap_matrix_45)

ng <- nrow(per_cap_matrix_45)
ndim <- nrow(pc_cm)
if(!all.equal(2*ng, ndim)){warning('dimensions not adding up')}

#### SET MORE SEIR INPUTS ####

# population sizes
pop_stratified <- c(risk_group_pop$pop - risk_group_pop$risk_population,
                    risk_group_pop$risk_population)

pop_vaccinated <- c(vaccinated_pop$vaccinated_population)

tot_pop <- sum(imd_age_pop$pop)
if(!all.equal(sum(pop_stratified), tot_pop)){warning('pop not adding up')}

init_infected_num <- 1000
init_infected_vec <- (pop_stratified - pop_vaccinated)*init_infected_num/(tot_pop-sum(pop_vaccinated))
if(!all.equal(sum(init_infected_vec), init_infected_num)){warning('init infected not adding up')}

#### RUN EACH SEASON ####

seasonal_seir_outputs <- list()

for(i in 1:length(years)){
  
  pars <- get(paste0('epid_parameters_s', i))
  
  time_series <- run_model(
    pop = pop_stratified,
    I0 = init_infected_vec,
    vacc = pop_vaccinated,
    cm = pc_cm,
    trans = pars$transmissibility,
    lat_per = pars$latent_period,
    inf_per = pars$infectious_period
    ) 
    
  time_series <- time_series[vaccinated_pop %>% mutate(imd_quintile=as.character(imd_quintile)) %>% 
                               select(age_grp, imd_quintile, risk_level, pop),
                             on = c('age_grp','imd_quintile','risk_level')]
  
  time_series[, start_date := pars$start_date]
  
  seasonal_seir_outputs[[i]] <- time_series
  
}

names(seasonal_seir_outputs) <- years

#### PLOT EXAMPLES ####

seasonal_seir_outputs[[1]] %>% 
  filter(imd_quintile == 1) %>% 
  ggplot() +
  geom_line(aes(t, value, col = compartment)) +
  facet_grid(age_grp ~ risk_level, scales = 'free')

seasonal_seir_outputs[[1]] %>% 
  filter(imd_quintile == 1) %>% 
  ggplot() +
  geom_line(aes(t, value/pop, col = compartment)) +
  facet_grid(age_grp ~ risk_level, scales = 'free') +
  theme_bw()

#### FINAL SIZE ####

plot_final_size <- function(k){
  
  seasonal_seir_outputs[[k]][t == max(t) & compartment == 'cumI',] %>% 
    ggplot() + 
    geom_bar(aes(x = age_grp, y = 100*value/pop, fill = imd_quintile),
             stat = 'identity', position = 'dodge') + 
    theme_bw() + 
    scale_fill_manual(values = imd_quintile_colors) +
    facet_grid(.~risk_level) + ylim(c(0, 100)) + 
    labs(x = 'Age group', y = 'Final size (%)', fill = 'IMD quintile') +
    ggtitle(names(seasonal_seir_outputs)[k])

}

final_size_plots <- map(.x = 1:length(years), .f = plot_final_size)

patchwork::wrap_plots(final_size_plots, nrow = 3)

#### SAVE DUMMY DATA ####

saveRDS(seasonal_seir_outputs, .args[3])


