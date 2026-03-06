## KNOWN PARAMETERS ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "known_parameters.rds")
) else commandArgs(trailingOnly = TRUE)

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

imd_age_pop <- imd_age_pop %>% 
  arrange(imd_quintile, age_grp)

#### MAKE PARAMETERS ####

## number of years of data
years <- 2023:2025 # 2023-24 to 2025-26

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

## vaccinated proportions

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

#### EPI PERIODS ####

epid_periods <- c(2, 3) # latent and infectious periods

#### VACCINE EFFICACY ####
## (age-dependent)
VE_pars <- c(0.70, 0.46) # currently just using the NGIV VE estimates with no mismatching

vaccination_efficacy <- c(
  rep(VE_pars[1], 7),
  rep(VE_pars[2], 2)
)

init_infected <- 1000

#### DELAYS ####

primary_care_delay <- 1
secondary_care_delay <- 2

#### MAKE INTO LIST ####

known_pars <- list(
  years = years,
  risk_group_pop = risk_group_pop,
  vaccinated_pop = vaccinated_pop,
  epid_periods = epid_periods,
  vaccination_efficacy = vaccination_efficacy,
  init_infected = init_infected,
  primary_care_delay = primary_care_delay,
  secondary_care_delay = secondary_care_delay
)

#### SAVE KNOWN PARAMETERS ####

saveRDS(known_pars, .args[2])

