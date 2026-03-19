## CREATE DUMMY DATA (INFECTIONS) ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "contact_matrix.rds"),
  file.path("data", "dummy_data", "known_parameters.rds"),
  file.path("data", "dummy_data", "unknown_parameters.rds"),
  file.path("data", "dummy_data", "dummy_infections.rds")
) else commandArgs(trailingOnly = TRUE)
  
source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','seir_model.R'))

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

## KNOWN PARAMETERS
known_pars <- readRDS(.args[3])
years <- known_pars$years

risk_group_pop <- known_pars$risk_group_pop
vaccinated_pop <- known_pars$vaccinated_pop
 
v_p <- vaccinated_pop %>% group_by(age_grp, imd_quintile) %>% 
  summarise(vaccinated_population = sum(vaccinated_population), pop = sum(pop)) %>% 
  mutate(v_p = vaccinated_population/pop) %>% 
  arrange(imd_quintile, age_grp) %>% 
  left_join(known_pars$vaccination_efficacy, by = 'age_grp') %>% 
  mutate(eff_v_p = VE*v_p) ## effectively vaccinated individuals

## UNKNOWN PARAMETERS
unknown_pars <- readRDS(.args[4])

## check R0
cat('\n')
for(i in 1:length(years)){
  pars <- unknown_pars[[paste0('epid_parameters_s', i)]]
  cat('R0 = ', R0_func(pars$susceptibility,
                pars$infectious_period,
                pars$transmissibility,
                cm), '\n', sep = '')
}
## check Reff
cat('\n')
for(i in 1:length(years)){
  pars <- unknown_pars[[paste0('epid_parameters_s', i)]]
  cat('Reff = ', R0_func((1 - v_p$eff_v_p)*rep(pars$susceptibility, 5),
                       pars$infectious_period,
                       pars$transmissibility,
                       cm), '\n', sep = '')
}
cat('\n')

#### EXPAND CONTACT MATRIX ####

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

vaccinated_pop <- vaccinated_pop %>% 
  left_join(known_pars$vaccination_efficacy, by = 'age_grp') %>% 
  mutate(effectively_vaccinated_population = VE*vaccinated_population)
pop_vaccinated <- c(vaccinated_pop$effectively_vaccinated_population)

tot_pop <- sum(imd_age_pop$pop)
if(!all.equal(sum(pop_stratified), tot_pop)){warning('pop not adding up')}

#### RUN EACH SEASON ####

seasonal_seir_outputs <- list()

for(i in 1:length(years)){
  
  pars <- unknown_pars[[paste0('epid_parameters_s', i)]]
  
  init_infected_num <- pars$init_infected
  init_infected_vec <- (pop_stratified - pop_vaccinated)*init_infected_num/(tot_pop-sum(pop_vaccinated))
  if(!all.equal(sum(init_infected_vec), init_infected_num)){warning('init infected not adding up')}
  
  time_series <- run_model(
    pop = pop_stratified,
    I0 = init_infected_vec,
    vacc = pop_vaccinated,
    cm = pc_cm,
    trans = pars$transmissibility,
    susc = pars$susceptibility,
    lat_per = pars$latent_period,
    inf_per = pars$infectious_period
    ) 
    
  time_series <- time_series[vaccinated_pop %>% mutate(imd_quintile=as.character(imd_quintile)) %>% 
                               select(age_grp, imd_quintile, risk_level, pop),
                             on = c('age_grp','imd_quintile','risk_level')]
  
  time_series[, start_date := pars$start_date]
  
  time_series <- time_series[t %in% 0:365] ## take data from the start of each day
  
  seasonal_seir_outputs[[i]] <- time_series
  
}

names(seasonal_seir_outputs) <- years

#### PLOT EXAMPLES ####
  
seasonal_seir_outputs[[1]] %>% 
  ggplot() +
  geom_line(aes(t, infections/pop, col = imd_quintile)) +
  scale_color_manual(values = imd_quintile_colors) +
  facet_grid(age_grp ~ risk_level, scales = 'free')

#### FINAL SIZE ####

plot_final_size <- function(k){
  
  seasonal_seir_outputs[[k]] %>% 
    group_by(age_grp, imd_quintile, risk_level, pop) %>% 
    summarise(infections = sum(infections)) %>% 
    ggplot() + 
    geom_bar(aes(x = age_grp, y = 100*infections/pop, fill = imd_quintile),
             stat = 'identity', position = 'dodge') + 
    theme_bw() + 
    scale_fill_manual(values = imd_quintile_colors) +
    facet_grid(.~risk_level) +
    labs(x = 'Age group', y = 'Final size (%)', fill = 'IMD quintile') +
    ggtitle(names(seasonal_seir_outputs)[k])

}

final_size_plots <- map(.x = 1:length(years), .f = plot_final_size)

patchwork::wrap_plots(final_size_plots, nrow = 3)

#### SAVE DUMMY DATA ####

saveRDS(seasonal_seir_outputs, .args[5])


