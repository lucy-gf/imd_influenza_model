## UNKNOWN PARAMETERS ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(patchwork))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "dummy_data", "unknown_parameters.rds")
) else commandArgs(trailingOnly = TRUE)

set.seed(60)

source(file.path("scripts","setup","colors.R"))

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

broad_ages <- data.table(
  age_grp = age_labels, 
  broad_age = c(rep('children', 3), rep('adults', 4), rep('older_adults', 2))
)

## number of years of data
years <- 2023:2025 # 2023-24 to 2025-26

#### EPIDEMIOLOGICAL PARAMETERS ####

epid_periods <- c(2, 3) # latent and infectious periods

susceptibility_matrix <- CJ(
  year = years,
  broad_age = c('children','adults','older_adults'),
  susceptibility = 0
)
susceptibility_matrix[broad_age == 'children', susceptibility := 0.6]
susceptibility_matrix[broad_age == 'adults', susceptibility := 0.3]
susceptibility_matrix[broad_age == 'older_adults', susceptibility := 0.45]
susceptibility_matrix[year == 2023, susceptibility := susceptibility*(0.95)]
susceptibility_matrix[year == 2023 & broad_age == 'children', susceptibility := susceptibility*(1.01)]
susceptibility_matrix[year == 2024, susceptibility := susceptibility*(1.05)]
susceptibility_matrix[year == 2023 & broad_age == 'older_adults', susceptibility := susceptibility*(1.01)]
susceptibility_matrix[year == 2025, susceptibility := susceptibility*(1.01)]
susceptibility_matrix[year == 2023 & broad_age == 'children', susceptibility := susceptibility*(0.98)]
if(nrow(susceptibility_matrix[susceptibility<0])>0){stop('Negative susceptibility')}

## make adults' susceptibility 1, everything else relative
susceptibility_adults <- susceptibility_matrix[broad_age=='adults']
setnames(susceptibility_adults,'susceptibility','adults_val')
susceptibility_adults[, broad_age := NULL]
susceptibility_matrix <- susceptibility_matrix[susceptibility_adults, on = c('year')]
susceptibility_matrix[, susceptibility := susceptibility/adults_val]
susceptibility_matrix[, adults_val := NULL]

## turn into age groups
susceptibility_long <- cross_join(susceptibility_matrix, broad_ages)[broad_age.x==broad_age.y]
susceptibility_long[, c('broad_age.x','broad_age.y') := NULL]
susceptibility_long$age_grp <- factor(susceptibility_long$age_grp, levels = age_labels)
setorder(susceptibility_long, year, age_grp)

epid_parameters_s1 <- list(
  susceptibility = susceptibility_long[year==years[1]]$susceptibility, # currently not age-dependent, one for each season
  transmissibility = 0.04,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[1]), "%d-%m-%Y"),
  init_infected = 300
)

epid_parameters_s2 <- list(
  susceptibility = susceptibility_long[year==years[2]]$susceptibility, # currently not age-dependent, one for each season
  transmissibility = 0.045,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[2]), "%d-%m-%Y"),
  init_infected = 400
)

epid_parameters_s3 <- list(
  susceptibility = susceptibility_long[year==years[3]]$susceptibility, # currently not age-dependent, one for each season
  transmissibility = 0.041,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[3]), "%d-%m-%Y"),
  init_infected = 200
)

#### REPORTING RATES ####

# TODO for now these are the same in each season

gp_rate <- c(1, 0.1, 2, 5, 3, 6)/50
hosp_rate <- c(2, 0.1, 2, 8, 3, 15)/500
## LOW RISK CHILDREN, ADULTS, OLDER ADULTS, 
## HIGH RISK CHILDREN, ADULTS, OLDER ADULTS

care_rate_df <- data.frame(
  broad_age = rep(unique(broad_ages$broad_age), 2),
  risk_level = rep(c('low','high'), each = 3),
  gp_rate = gp_rate,
  hosp_rate = hosp_rate
)

care_rate_age_df <- data.table(cross_join(
  care_rate_df,
  broad_ages))[broad_age.x == broad_age.y,]
care_rate_age_df[, c('broad_age.x','broad_age.y') := NULL]

imd_spline_pars <- data.table(
  primary = c(-0.8, 0.7),
  secondary = c(-0.3, 0.2)
)
rel_imd_rep_rates <- data.frame(imd_quintile = 1:5,
                                rel_primary_rates = imd_spline(imd_spline_pars$primary),
                                rel_secondary_rates = imd_spline(imd_spline_pars$secondary))

care_rate_imd_df <- cross_join(
  care_rate_age_df,
  rel_imd_rep_rates
  ) %>% 
  mutate(gp_rate = gp_rate*rel_primary_rates,
         hosp_rate = hosp_rate*rel_secondary_rates)

care_rate_imd_df$age_grp <- factor(care_rate_imd_df$age_grp, levels = age_labels) 

ratep1 <- care_rate_imd_df %>% 
  pivot_longer(c(gp_rate, hosp_rate)) %>% 
  ggplot() + 
  geom_line(aes(age_grp, value, group = interaction(name,imd_quintile), 
                col = as.factor(imd_quintile)), lwd = 0.8) +
  theme_bw() + ylim(c(0,NA)) + facet_grid(name ~ risk_level, scales = 'free') + 
  scale_color_manual(values = imd_quintile_colors) +
  labs(y='Rate', col = 'IMD quintile'); ratep1

ratep2 <- rel_imd_rep_rates %>% 
  pivot_longer(!imd_quintile) %>%
  ggplot() + 
  geom_line(aes(imd_quintile, value, group = name, lty = name), lwd = 0.8) +
  theme_bw() + ylim(c(0,NA)) +
  labs(y='Relative rate (baseline = IMD 3)', lty = 'Care setting'); ratep2

ratep1 + ratep2 + plot_layout(nrow = 1, widths = c(2,1))
ggsave(file.path('output','figures','dummy_mcmc','reporting_rates.png'), width = 16, height = 7)

#### MAKE INTO LIST ####

unknown_pars <- list(
  epid_parameters_s1 = epid_parameters_s1,
  epid_parameters_s2 = epid_parameters_s2,
  epid_parameters_s3 = epid_parameters_s3,
  care_rates = care_rate_imd_df,
  imd_spline_pars = imd_spline_pars,
  primary_care_rates = gp_rate,
  secondary_care_rates = hosp_rate
)

#### SAVE UNKNOWN PARAMETERS ####

saveRDS(unknown_pars, .args[2])


