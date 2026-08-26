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

broad_ages_susc <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages(
    'children', 
    'adults',
    'older_adults',
    age_labels
    )
)
broad_ages_care_rates <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages(
    'children', 
    'adults',
    'older_adults',
    age_labels
  )
)


## number of years of data
years <- 2023:2025 # 2023-24 to 2025-26

#### EPIDEMIOLOGICAL PARAMETERS ####

epid_periods <- c(2, 3) # latent and infectious periods

susceptibility_long <- cross_join(
  CJ(year = years, strain = c("A","B")),
  data.table(age_grp = age_labels, 
             broad_age = fcn_assign_ages('children','adults','older_adults', age_labels),
             susceptibility = fcn_assign_ages(0.6, 0.3, 0.45, age_labels))
)
susceptibility_long[year == 2023, susceptibility := susceptibility*(0.95)]
susceptibility_long[year == 2023 & broad_age == 'children', susceptibility := susceptibility*(1.01)]
susceptibility_long[year == 2023 & broad_age == 'older_adults', susceptibility := susceptibility*(0.96)]
susceptibility_long[year == 2024, susceptibility := susceptibility*(1.05)]
susceptibility_long[year == 2024 & broad_age == 'older_adults', susceptibility := susceptibility*(1.02)]
susceptibility_long[year == 2025, susceptibility := susceptibility*(1.01)]
susceptibility_long[year == 2025 & broad_age == 'children', susceptibility := susceptibility*(0.98)]
susceptibility_long[strain == 'B', susceptibility := susceptibility*(0.5)]
susceptibility_long[strain == 'B' & broad_age == 'older_adults', susceptibility := susceptibility*(0.8)]
if(nrow(susceptibility_long[susceptibility<0])>0){stop('Negative susceptibility')}

## make adults' susceptibility 1, everything else relative
susceptibility_adults <- susceptibility_long[broad_age=='adults']
setnames(susceptibility_adults,'susceptibility','adults_val')
susceptibility_adults[, c('broad_age', 'age_grp') := NULL]
susceptibility_long <- susceptibility_long[unique(susceptibility_adults), on = c('year', 'strain')]
susceptibility_long[, susceptibility := susceptibility/adults_val]
susceptibility_long[, c('broad_age','adults_val') := NULL]
susceptibility_long$age_grp <- factor(susceptibility_long$age_grp, levels = age_labels)
setorder(susceptibility_long, year, strain, age_grp)

susceptibility_long %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, y = susceptibility, group = interaction(year, strain),
                col = as.factor(year)), lwd = 1) +
  geom_point(aes(x = age_grp, y = susceptibility, group = interaction(year, strain)), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, y = susceptibility, group = interaction(year, strain),
                 col = as.factor(year), shape = strain), 
             stroke=1.5, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values = season_colors) +
  theme_bw() + labs(x = 'Age group', col = 'Season start',
                    y = 'VE against hospitalisation') +
  facet_grid(strain ~.) + 
  theme(text = element_text(size = 14))

epid_parameters_s1_A <- list(
  susceptibility = susceptibility_long[year==years[1] & strain == "A"]$susceptibility, 
  transmissibility = 0.041,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[1]), "%d-%m-%Y"),
  init_infected = 300
)

epid_parameters_s1_B <- list(
  susceptibility = susceptibility_long[year==years[1] & strain == "B"]$susceptibility, 
  transmissibility = 0.039,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[1]), "%d-%m-%Y"),
  init_infected = 300
)

epid_parameters_s2_A <- list(
  susceptibility = susceptibility_long[year==years[2] & strain == "A"]$susceptibility, 
  transmissibility = 0.041,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[2]), "%d-%m-%Y"),
  init_infected = 400
)

epid_parameters_s2_B <- list(
  susceptibility = susceptibility_long[year==years[2] & strain == "B"]$susceptibility, 
  transmissibility = 0.036,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[2]), "%d-%m-%Y"),
  init_infected = 400
)

epid_parameters_s3_A <- list(
  susceptibility = susceptibility_long[year==years[3] & strain == "A"]$susceptibility, 
  transmissibility = 0.041,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[3]), "%d-%m-%Y"),
  init_infected = 200
)

epid_parameters_s3_B <- list(
  susceptibility = susceptibility_long[year==years[3] & strain == "B"]$susceptibility, 
  transmissibility = 0.033,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0('01-09-', years[3]), "%d-%m-%Y"),
  init_infected = 200
)

#### REPORTING RATES ####

# TODO for now these are the same in each season and for each strain

gp_rate <- c(1, 0.1, 2, 5, 3, 6)/20
hosp_rate <- c(2, 0.1, 2, 8, 3, 15)/100
## LOW RISK CHILDREN, ADULTS, OLDER ADULTS, 
## HIGH RISK CHILDREN, ADULTS, OLDER ADULTS

care_rate_df <- data.frame(
  broad_age = rep(unique(broad_ages_care_rates$broad_age), 2),
  risk_level = rep(c('low','high'), each = 3),
  gp_rate = gp_rate,
  hosp_rate = hosp_rate
)

care_rate_age_df <- data.table(cross_join(
  care_rate_df,
  broad_ages_care_rates))[broad_age.x == broad_age.y,]
care_rate_age_df[, c('broad_age.x','broad_age.y') := NULL]

imd_spline_pars <- data.table(
  primary = c(-0.1, 0.12),
  secondary = c(-0.07, 0.05)
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
  geom_line(aes(age_grp, value, group = imd_quintile, 
                col = as.factor(imd_quintile)), lwd = 0.8) +
  geom_point(aes(x = age_grp, group = imd_quintile, 
                 y = value), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = imd_quintile, 
                 col = as.factor(imd_quintile), y = value), stroke=1.5, size = 3, shape = 1) +
  theme_bw() + ylim(c(0,NA)) + facet_grid(name ~ risk_level, scales = 'free') + 
  scale_color_manual(values = imd_quintile_colors) +
  labs(y='Rate', col = 'IMD quintile'); ratep1

ratep2 <- rel_imd_rep_rates %>% 
  pivot_longer(!imd_quintile) %>%
  ggplot() + 
  geom_line(aes(imd_quintile, value, group = name, lty = name), lwd = 0.8) +
  geom_point(aes(imd_quintile, value, group = name), 
             color = 'white', size = 3) +
  geom_point(aes(imd_quintile, value, group = name), 
             shape = 1, stroke = 2, size = 3) +
  theme_bw() + ylim(c(0.8,1.2)) +
  scale_linetype_manual(labels = c("Primary care", "Secondary care"),
                        values = c(1,2)) +
  theme(text = element_text(size = 12)) +
  labs(y='Relative rate (baseline = IMD 3)', lty = 'Care setting',
       x = 'IMD quintile'); ratep2

ratep1 + ratep2 + plot_layout(nrow = 1, widths = c(2,1))
ggsave(file.path('output','figures','dummy_infections','reporting_rates.png'), width = 16, height = 7)

#### MAKE INTO LIST ####

unknown_pars <- list(
  epid_parameters_s1_A = epid_parameters_s1_A,
  epid_parameters_s2_A = epid_parameters_s2_A,
  epid_parameters_s3_A = epid_parameters_s3_A,
  epid_parameters_s1_B = epid_parameters_s1_B,
  epid_parameters_s2_B = epid_parameters_s2_B,
  epid_parameters_s3_B = epid_parameters_s3_B,
  care_rates = care_rate_imd_df,
  imd_spline_pars = imd_spline_pars,
  primary_care_rates = gp_rate,
  secondary_care_rates = hosp_rate
)

#### SAVE UNKNOWN PARAMETERS ####

saveRDS(unknown_pars, .args[2])

