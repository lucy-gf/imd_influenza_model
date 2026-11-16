## PRODUCE TRUE VALUES VECTOR ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  c(file.path("data", "inputs", "imd_age_pop.rds"),
    file.path("data", "inputs", "contact_matrix.rds"),
    file.path("data", "inputs", "subtype_years.rds"),
    file.path("data", "dummy_data", "known_parameters.rds"),
    file.path("data", "dummy_data", "unknown_parameters.rds"),
    file.path("data", "dummy_data", "mcmc_true_values.rds")
  )) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

set.seed(60)

figure_filename <- function(string){
  paste0('figures/dummy_mcmc/parameters/',string,'.png')
}

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

broad_ages <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages('children','adults','older_adults', age_labels)
)

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

#### SUBTYPE-SEASONS #### 
ukhsa_subtype_data <- readRDS(.args[3])
subtype_seasons <- ukhsa_subtype_data %>% 
  mutate(year = as.numeric(substr(season,1,4))) %>% 
  select(season, subtype, year) %>% unique() %>% 
  mutate(epidemic = year - min(year) + 1) %>% 
  group_by(epidemic) %>% mutate(epidemic_of_season = 1:n())

#### KNOWN PARAMETERS #### 
known_pars <- readRDS(.args[4])
years <- known_pars$years

#### UNKNOWN PARAMETERS #### 
unknown_pars <- readRDS(.args[5])

## true epid parameters being fitted
epid_pars <- data.frame()
for(i in 1:nrow(subtype_seasons)){
  
  up_ss <- unknown_pars$epid_parameters %>% filter(
    subtype == subtype_seasons$subtype[i],
    year == subtype_seasons$year[i]
  )
  up_care_base_ss <- unknown_pars$base_care_rates %>% filter(
    subtype == subtype_seasons$subtype[i],
    year == subtype_seasons$year[i]
  )
  up_imd_spline_ss <- unknown_pars$imd_spline_pars %>% filter(
    subtype == subtype_seasons$subtype[i],
    year == subtype_seasons$year[i]
  )
  
  epid_pars <- rbind(epid_pars,
                     data.frame(epidemic = subtype_seasons$epidemic[i],
                                season = subtype_seasons$season[i],
                                subtype = subtype_seasons$subtype[i],
                                year = subtype_seasons$year[i],
                                transmissibility = up_ss$transmissibility,
                                children_susceptibility = up_ss$children_susceptibility,
                                adults_susceptibility = up_ss$adults_susceptibility,
                                older_adults_susceptibility = up_ss$older_adults_susceptibility,
                                init_infected = log10(up_ss$init_infected),
                                primary_care_rate_children_low = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'low'))$gp_rate,
                                primary_care_rate_adults_low = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'low'))$gp_rate,
                                primary_care_rate_older_adults_low = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'low'))$gp_rate,
                                primary_care_rate_children_high = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'high'))$gp_rate,
                                primary_care_rate_adults_high = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'high'))$gp_rate,
                                primary_care_rate_older_adults_high = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'high'))$gp_rate,
                                secondary_care_rate_children_low = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'low'))$hosp_rate,
                                secondary_care_rate_adults_low = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'low'))$hosp_rate,
                                secondary_care_rate_older_adults_low = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'low'))$hosp_rate,
                                secondary_care_rate_children_high = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'high'))$hosp_rate,
                                secondary_care_rate_adults_high = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'high'))$hosp_rate,
                                secondary_care_rate_older_adults_high = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'high'))$hosp_rate,
                                imd_spline_primary_1 = up_imd_spline_ss$primary[1],
                                imd_spline_primary_2 = up_imd_spline_ss$primary[2],
                                imd_spline_secondary_1 = up_imd_spline_ss$secondary[1],
                                imd_spline_secondary_2 = up_imd_spline_ss$secondary[2],
                                Reff = R0_func(susceptibility = fcn_assign_ages(up_ss$children_susceptibility,
                                                                                up_ss$adults_susceptibility,
                                                                                up_ss$older_adults_susceptibility,
                                                                                age_labels),
                                               inf_period = known_pars$epid_periods[2],
                                               beta_in = up_ss$transmissibility,
                                               cm_in = cm)
                     ))
}

epid_pars <- epid_pars %>% pivot_longer(!c(epidemic,subtype,season,year)) 

epid_pars <- epid_pars %>% 
  left_join(subtype_seasons, by = c('epidemic','season','subtype','year')) %>% 
  mutate(par_name = case_when(
    grepl('imd_spline', name) ~ paste0(name, '_', year),
    T ~ paste0(name, '_', year, '_', subtype)),
    joining_name =  case_when(
      grepl('imd_spline', name) ~ name,
      T ~ paste0(name, '_epid_', epidemic_of_season)))

write_rds(epid_pars, .args[6])


