## KNOWN PARAMETERS ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "population", "risk_group_population_data.rds"),
  file.path("data", "dummy_data", "known_parameters.rds")
) else commandArgs(trailingOnly = TRUE)

if(!dir.exists(file.path("data", "dummy_data"))){dir.create(file.path("data", "dummy_data"))}

source(file.path('scripts','setup','colors.R'))

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

#### PROPORTIONS IN RISK GROUPS ####

## risk group proportions

USING_TRUE_RISK_DAT <- T

if(USING_TRUE_RISK_DAT){
  
  ## loading in true risk group data
  risk_group_dat <- readRDS(.args[2]) %>% 
    mutate(.by = interval,
           proportion = count/sum(count)) %>% 
    filter(risk_group %like% 'high') %>% arrange(lower)
  
  ## matching up age groups via midpoint
  age_labels_midpoints_list <- (str_split(gsub('[+]','-100',age_labels), pattern = '-'))
  age_labels_midpoints <- c(); corresponding_intervals <- c(); risk_group_age_vector <- c()
  for(interval in 1:nage){
    age_labels_midpoints[interval] <- mean(as.numeric(age_labels_midpoints_list[[interval]]))
    corresponding_intervals[interval] <- sum(age_labels_midpoints[interval] > risk_group_dat$lower)
    risk_group_age_vector[interval] <- risk_group_dat$proportion[corresponding_intervals[interval]]
  }
  
  ## add some variation by IMD
  rgiv_vec <- c(0.15, 0.25, 0.08) # children, adults, older adults
  risk_group_imd_variation_multipliers <- c(
    rep(seq(1 + 2*rgiv_vec[1], 1 - 2*rgiv_vec[1], by = -rgiv_vec[1]), 3),
    rep(seq(1 + 2*rgiv_vec[2], 1 - 2*rgiv_vec[2], by = -rgiv_vec[2]), 4),
    rep(seq(1 + 2*rgiv_vec[3], 1 - 2*rgiv_vec[3], by = -rgiv_vec[3]), 2)
  )
  
  risk_group_imd_variation <- c()
  for(age in 1:nage){
    for(imd in 1:nimd){
      risk_group_imd_variation[imd + nimd*(age - 1)] <- 
        risk_group_imd_variation_multipliers[imd + nimd*(age - 1)]*
        risk_group_age_vector[age]
      # cat(imd, ' ', age, ' ', imd + nimd*(age - 1),'\n')
    }
  }
  
}else{
  
  ## dummy data
  risk_group_kids <- 0.12
  risk_group_adults <- 0.08
  risk_group_elderly <- 0.4
  
  risk_group_age_vector <- c(
    rep(risk_group_kids, 3),
    rep(risk_group_adults, 4),
    rep(risk_group_elderly, 2)
  )
  
  ## add some variation by IMD
  rgiv_vec <- c(0.15, 0.25, 0.08) # children, adults, older adults
  risk_group_imd_variation <- c(
    rep(risk_group_kids*seq(1 + 2*rgiv_vec[1], 1 - 2*rgiv_vec[1], by = -rgiv_vec[1]), 3),
    rep(risk_group_adults*seq(1 + 2*rgiv_vec[2], 1 - 2*rgiv_vec[2], by = -rgiv_vec[2]), 4),
    rep(risk_group_elderly*seq(1 + 2*rgiv_vec[3], 1 - 2*rgiv_vec[3], by = -rgiv_vec[3]), 2)
  )
  
}

risk_group_pop <- data.frame(risk_proportion = risk_group_imd_variation) %>% 
  mutate(age_grp = rep(age_labels, each = nimd),
         imd_quintile = rep(1:nimd, nage)) %>% 
  left_join(imd_age_pop, by = c('age_grp','imd_quintile')) %>% 
  mutate(risk_population = round(risk_proportion*pop)) %>% 
  arrange(imd_quintile, age_grp)

risk_group_pop %>% 
  ggplot() + geom_line(aes(age_grp, risk_population, col=imd_quintile,
                           group=imd_quintile))

risk_group_pop %>% 
  ggplot() + geom_line(aes(age_grp, risk_proportion, col=imd_quintile,
                           group=imd_quintile)) 

#### VACCINATION COVERAGE #### 

# https://www.gov.uk/government/statistics/seasonal-influenza-vaccine-uptake-in-children-of-school-age-monthly-data-2025-to-2026
vaccination_coverage_kids <- 0.55
# https://www.gov.uk/government/statistics/seasonal-influenza-vaccine-uptake-in-gp-patients-monthly-data-2025-to-2026
vaccination_coverage_elderly <- 0.75
vaccination_coverage_risk <- 0.4

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
    T ~ round(vaccinated_proportion*pop)
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

vaccinated_pop %>% 
  mutate(key_group = case_when(
    age_grp %in% c('0-4','5-11','12-17') ~ 'Children',
    age_grp %notin% c('0-4','5-11','12-17','70-79','80+') & risk_level == 'high' ~ 'Risk group (18-69)',
    age_grp %in% c('70-79','80+') ~ 'Older adults',
    T~ NA
  )) %>% 
  group_by(key_group, imd_quintile) %>% 
  summarise(vaccinated_population = sum(vaccinated_population),
            pop = sum(pop)) %>% filter(!is.na(key_group)) %>% 
  ggplot() + 
  geom_bar(aes(x = key_group, group = imd_quintile, 
               fill = as.factor(imd_quintile), y = 100*vaccinated_population/pop),
           stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = imd_quintile_colors) + 
  theme_bw() + labs(x = '', fill = 'IMD quintile', 
                    y = 'Simulated vaccine uptake (%)') +
  theme(text = element_text(size = 12))

vaccinated_pop %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, group = as.factor(imd_quintile), 
               col = as.factor(imd_quintile), y = risk_proportion), lwd = 1) +
  geom_point(aes(x = age_grp, group = as.factor(imd_quintile), 
                 col = as.factor(imd_quintile), y = risk_proportion), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = as.factor(imd_quintile), 
                col = as.factor(imd_quintile), y = risk_proportion), 
             shape = 1, stroke=2, size = 3) +
  scale_color_manual(values = imd_quintile_colors) + ylim(c(0,NA)) +
  theme_bw() + labs(x = 'Age group', col = 'IMD quintile', 
                    y = 'Simulatedercentage in clinical risk group') +
  theme(text = element_text(size = 14))

#### EPI PERIODS ####

epid_periods <- c(2, 3) # latent and infectious periods

#### VACCINE EFFICACY ####
## (age-dependent, annual, eventually strain-specific)
# TODO Make strain-specific

#### AGAINST INFECTION ####
# 2023/24: https://onlinelibrary.wiley.com/doi/epdf/10.1111/irv.70194
# 2024/25: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#secondary-care-surveillance

vaccination_efficacy_infection <- CJ(
  age_grp = age_labels,    
  start_of_season = years,
  VE_INF = 0
)

children_ages <- c('0-4','5-11','12-17')
adult_ages <- c('18-25','26-34','35-49','50-69')
older_adult_ages <- c('70-79','80+')

## FOR NOW USING MADE UP DATA
## TODO UPDATE WHEN DATA AVAILABLE
vaccination_efficacy_infection[
  age_grp %in% children_ages, VE_INF := 0.45][
    age_grp %in% adult_ages, VE_INF := 0.35][
      age_grp %in% older_adult_ages, VE_INF := 0.3]

#### AGAINST HOSPITALISATION ####
# 2023/24: https://onlinelibrary.wiley.com/doi/epdf/10.1111/irv.70194
# 2024/25: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#secondary-care-surveillance
 
vaccination_efficacy_hospitalisation <- CJ(
  age_grp = age_labels,    
  start_of_season = years,
  VE_HOSP = 0
)

children_ages <- c('0-4','5-11','12-17')
adult_ages <- c('18-25','26-34','35-49','50-69')
older_adult_ages <- c('70-79','80+')

vaccination_efficacy_hospitalisation[
  age_grp %in% children_ages & start_of_season == 2023, VE := 0.56][
  age_grp %in% adult_ages & start_of_season == 2023, VE := 0.38][
  age_grp %in% older_adult_ages & start_of_season == 2023, VE := 0.18]

vaccination_efficacy_hospitalisation[
  age_grp %in% children_ages & start_of_season == 2024, VE := 0.62][
  age_grp %in% adult_ages & start_of_season == 2024, VE := 0.46][
  age_grp %in% older_adult_ages & start_of_season == 2024, VE := 0.40]

# for now, use 2024/25 values for 2025/26
vaccination_efficacy_hospitalisation[
  age_grp %in% children_ages & start_of_season == 2025, VE := 0.62][
  age_grp %in% adult_ages & start_of_season == 2025, VE := 0.46][
  age_grp %in% older_adult_ages & start_of_season == 2025, VE := 0.40]

#### DELAYS ####
## (in weeks)
primary_care_delay <- 1
secondary_care_delay <- 2

#### OPENSAFELY COVERAGE ####

## for now assuming that this is 42% everywhere,
## but in real-world model this will vary across subgroups

proportion_observed <- CJ(
  age_grp = age_labels,
  imd_quintile = 1:5,
  risk_level = c('high','low'),
  OS_COVERAGE = 0.42
) 
 
#### MAKE INTO LIST ####

known_pars <- list(
  years = years,
  risk_group_pop = risk_group_pop,
  vaccinated_pop = vaccinated_pop,
  epid_periods = epid_periods,
  vaccination_efficacy_infection = vaccination_efficacy_infection,
  vaccination_efficacy_hospitalisation = vaccination_efficacy_hospitalisation,
  primary_care_delay = primary_care_delay,
  secondary_care_delay = secondary_care_delay,
  proportion_observed = proportion_observed
)

#### SAVE KNOWN PARAMETERS ####

saveRDS(known_pars, .args[3])

