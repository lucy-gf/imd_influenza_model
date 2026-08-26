## KNOWN PARAMETERS ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(patchwork))
suppressMessages(require(tidyverse))
suppressMessages(require(readODS))
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
source(file.path('scripts','setup','age_grp_assignment.R'))

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
  
  risk_group_age_vector <- fcn_assign_ages(
    risk_group_kids, 
    risk_group_adults,
    risk_group_elderly,
    age_labels
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
  ggplot() + geom_line(aes(age_grp, risk_population, col=as.factor(imd_quintile),
                           group=imd_quintile), lwd = 0.8) + 
  scale_color_manual(values = imd_quintile_colors) + theme_bw() + 
  labs(x = '', y = 'Proportion in risk group', col = 'IMD')

risk_group_pop %>% 
  ggplot() + geom_line(aes(age_grp, risk_proportion, col=as.factor(imd_quintile),
                           group=imd_quintile), lwd = 0.8) + 
  scale_color_manual(values = imd_quintile_colors) + theme_bw() + 
  labs(x = '', y = 'Proportion in risk group', col = 'IMD')

#### VACCINATION COVERAGE #### 

# https://www.gov.uk/government/statistics/seasonal-influenza-vaccine-uptake-in-children-of-school-age-monthly-data-2025-to-2026
vaccination_coverage_kids <- 0.55
# https://www.gov.uk/government/statistics/seasonal-influenza-vaccine-uptake-in-gp-patients-monthly-data-2025-to-2026
vaccination_coverage_elderly <- 0.75
vaccination_coverage_risk <- 0.4

vaccination_coverage_age_vector <- fcn_assign_ages(
  vaccination_coverage_kids, 
  vaccination_coverage_risk,
  vaccination_coverage_elderly,
  age_labels
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
  select(!c(risk_population, vcav, vciv)) %>% 
  mutate(nrow = n()) %>% 
  ## duplicate across years for now
  slice(rep(row_number(), 3)) %>% 
  mutate(start_of_season = rep(years, each = nrow[1])) %>% select(!nrow) %>% 
  mutate(vaccinated_population = case_when(
    age_grp %in% c('18-25', '26-34', '35-49', '50-69') & risk_level == 'low' ~ 0,
    T ~ round(vaccinated_proportion*pop)
  )) %>% 
  ## ENSURE ALL POP ORDERS ARE IMD THEN AGE
  arrange(start_of_season, desc(risk_level), imd_quintile, age_grp)

vaccinated_pop %>% group_by(age_grp,imd_quintile) %>% 
  summarise(vaccinated_population=sum(vaccinated_population)) %>% 
  ggplot() + geom_line(aes(age_grp, vaccinated_population, col=imd_quintile,
                           group=imd_quintile))

vaccinated_pop %>% 
  ggplot() + geom_line(aes(age_grp, vaccinated_population, col=imd_quintile,
                           lty = risk_level, group=interaction(imd_quintile,risk_level)))

vacc_plot <- vaccinated_pop %>% 
  ggplot() +
  geom_line(aes(x = age_grp, group = interaction(as.factor(imd_quintile), risk_level), 
                col = as.factor(imd_quintile), y = vaccinated_population/pop), lwd = 1) +
  geom_point(aes(x = age_grp, group = interaction(as.factor(imd_quintile), risk_level), 
                 y = vaccinated_population/pop, shape = risk_level), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = interaction(as.factor(imd_quintile), risk_level), 
                 col = as.factor(imd_quintile), y = vaccinated_population/pop, 
                 shape = risk_level), stroke=1.5, size = 3) +
  scale_color_manual(values = imd_quintile_colors) + ylim(c(0,NA)) +
  scale_shape_manual(values = c(1, 2)) +
  theme_bw() + labs(x = 'Age group', col = 'IMD quintile', 
                    y = 'Simulated vaccination coverage',
                    shape = 'Risk level') +
  theme(text = element_text(size = 14)); vacc_plot

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
  theme(text = element_text(size = 14))

risk_plot <- vaccinated_pop %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, group = as.factor(imd_quintile), 
               col = as.factor(imd_quintile), y = risk_proportion), lwd = 1) +
  geom_point(aes(x = age_grp, group = as.factor(imd_quintile), y = risk_proportion), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = as.factor(imd_quintile), 
                col = as.factor(imd_quintile), y = risk_proportion), 
             shape = 1, stroke=2, size = 3) +
  scale_color_manual(values = imd_quintile_colors) + ylim(c(0,NA)) +
  theme_bw() + labs(x = 'Age group', col = 'IMD quintile', 
                    y = 'Simulated percentage in clinical risk group') +
  theme(text = element_text(size = 14),
        legend.position = 'none')

vacc_plot + risk_plot + plot_layout(nrow = 2, guides = 'collect')

ggsave(file.path('output','figures','dummy_infections','dummy_vacc_risk_props.png'),
       width = 10, height = 10)

#### VACCINE EFFICACY ####
## (age-dependent, annual, eventually strain-specific)
# TODO Make strain-specific

#### AGAINST INFECTION ####
# 2023/24: https://onlinelibrary.wiley.com/doi/epdf/10.1111/irv.70194
# 2024/25: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#secondary-care-surveillance

## FOR NOW USING MADE UP DATA
## TODO UPDATE WHEN DATA AVAILABLE

vaccination_efficacy_infection_A <- cross_join(
  data.table(age_grp = age_labels,  
             VE_INF = fcn_assign_ages(
               0.45, 
               0.35,
               0.3,
               age_labels
             )),
  data.table(start_of_season = years, strain = "A")
) 

vaccination_efficacy_infection_B <- cross_join(
  data.table(age_grp = age_labels,  
             VE_INF = fcn_assign_ages(
               0.45, 
               0.35,
               0.3,
               age_labels
             )),
  data.table(start_of_season = years, strain = "B")
) 

vaccination_efficacy_infection <- rbind(
  vaccination_efficacy_infection_A,
  vaccination_efficacy_infection_B
)

#### AGAINST HOSPITALISATION ####
# 2023/24: https://onlinelibrary.wiley.com/doi/epdf/10.1111/irv.70194
# 2024/25: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#secondary-care-surveillance
 
# from: https://doi.org/10.1111/irv.70194

VE_2_17_A_2022_2023 <- 0.60
VE_2_17_A_2023_2024 <- 0.50
VE_2_17_B_2022_2023 <- 0.87
VE_2_17_B_2023_2024 <- 0.84
VE_18_64_A_2022_2023 <- 0.28
VE_18_64_A_2023_2024 <- 0.34
VE_18_64_B_2022_2023 <- 0.53
VE_18_64_B_2023_2024 <- 0.70
VE_65_A_2022_2023 <- 0.25
VE_65_A_2023_2024 <- 0.17
VE_65_B_2022_2023 <- 0.28
VE_65_B_2023_2024 <- 0.39

# from: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#vaccination
# (Figure 56)

ukhsa_dat_2024_25 <- read_ods(file.path("data","ukhsa","annual_influenza_2024_2025.ods"),
                              sheet = 59, skip = 3)
colnames(ukhsa_dat_2024_25) <- c('age_group','flu_subtype','NA1','NA2','NA3','NA4','VE','VE_CIL','VE_CIU')
ukhsa_dat_2024_25 <- ukhsa_dat_2024_25 %>% mutate(VE = VE/100)

VE_2_17_A_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'Any Influenza A'))$VE
VE_2_17_B_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'Influenza B'))$VE
VE_18_64_A_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'Any Influenza A'))$VE
VE_18_64_B_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'Influenza B'))$VE
VE_65_A_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '65', flu_subtype %like% 'Any Influenza A'))$VE
VE_65_B_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '65', flu_subtype %like% 'Influenza B'))$VE

# from: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2025-to-2026/influenza-in-the-uk-annual-epidemiological-report-winter-2025-to-2026#vaccination
# (Figure 56)

ukhsa_dat_2025_26 <- read_ods(file.path("data","ukhsa","annual_influenza_2025_2026.ods"),
                              sheet = 59, skip = 3)
colnames(ukhsa_dat_2025_26) <- c('age_group','flu_subtype','NA1','NA2','NA3','NA4','VE')
ukhsa_dat_2025_26 <- ukhsa_dat_2025_26 %>% 
  mutate(VE = as.numeric(sub("%.*", "", VE))/100,
         flu_subtype = gsub(substr(ukhsa_dat_2025_26$flu_subtype[2], 4, 4), '', flu_subtype))

VE_2_17_A_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'AnyInfluenzaA'))$VE
VE_2_17_B_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'InfluenzaB'))$VE
VE_18_64_A_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'AnyInfluenzaA'))$VE
VE_18_64_B_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'InfluenzaB'))$VE
VE_65_A_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '65', flu_subtype %like% 'AnyInfluenzaA'))$VE
VE_65_B_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '65', flu_subtype %like% 'InfluenzaB'))$VE

vaccination_efficacy_hospitalisation <- CJ(strain = c('A','B'),
                                           start_of_season = years,
                                           age_grp = age_labels,
                                           VE_HOSP = 0)

for(strain_i in c('A','B')){
  for(season_i in years){
    vaccination_efficacy_hospitalisation[
      strain == strain_i & start_of_season == season_i, 
      VE_HOSP := fcn_assign_ages(
      get(paste0('VE_2_17_', strain_i, '_', season_i, '_', season_i + 1)), 
      get(paste0('VE_18_64_', strain_i, '_', season_i, '_', season_i + 1)),
      get(paste0('VE_65_', strain_i, '_', season_i, '_', season_i + 1)),
      age_labels
    )]
  }
}

vaccination_efficacy_hospitalisation %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, strain),
                col = as.factor(start_of_season)), lwd = 1) +
  geom_point(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, strain)), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, strain),
                 col = as.factor(start_of_season), shape = strain), 
             stroke=1.5, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values = season_colors) +
  theme_bw() + labs(x = 'Age group', col = 'Season start',
                    y = 'VE against hospitalisation') +
  facet_grid(strain ~.) + 
  theme(text = element_text(size = 14))

ggsave(file.path('output','figures','dummy_infections','ukhsa_ve_hospitalisation.png'),
       width = 10, height = 8)

## join together with vaccination coverage data

vaccinated_data <- rbind(vaccinated_pop %>% mutate(strain = 'A'),
                         vaccinated_pop %>% mutate(strain = 'B')) %>% 
  left_join(vaccination_efficacy_infection, by = c('age_grp','start_of_season', 'strain')) %>% 
  left_join(vaccination_efficacy_hospitalisation, by = c('age_grp','start_of_season', 'strain')) %>% 
  mutate(effectively_vaccinated_population = round(VE_INF*vaccinated_population))

#### EPI PERIODS ####

epid_periods <- c(2, 3) # latent and infectious periods

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
  vaccinated_data = vaccinated_data,
  epid_periods = epid_periods,
  primary_care_delay = primary_care_delay,
  secondary_care_delay = secondary_care_delay,
  proportion_observed = proportion_observed
)

#### SAVE KNOWN PARAMETERS ####

saveRDS(known_pars, .args[3])

