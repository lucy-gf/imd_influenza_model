## UNKNOWN PARAMETERS ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
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
susceptibility_matrix[year == 2024, susceptibility := susceptibility*(1.05)]
susceptibility_matrix[year == 2025, susceptibility := susceptibility*(1.01)]
# susceptibility_matrix[, susceptibility := (1 - runif(n=nrow(susceptibility_matrix),-0.05,0.05))*susceptibility]
if(nrow(susceptibility_matrix[susceptibility<0])>0){stop('Negative susceptibility')}
susceptibility_long <- cross_join(susceptibility_matrix, broad_ages)[broad_age.x==broad_age.y]
susceptibility_long[, c('broad_age.x','broad_age.y') := NULL]
susceptibility_long$age_grp <- factor(susceptibility_long$age_grp, levels = age_labels)
setorder(susceptibility_long, year, age_grp)

epid_parameters_s1 <- list(
  susceptibility = susceptibility_long[year==years[1]]$susceptibility, # currently not age-dependent, one for each season
  transmissibility = 0.15,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[1], '35', '1'), "%Y%W%w")
)

epid_parameters_s2 <- list(
  susceptibility = susceptibility_long[year==years[2]]$susceptibility, # currently not age-dependent, one for each season
  transmissibility = 0.13,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[2], '36', '1'), "%Y%W%w")
)

epid_parameters_s3 <- list(
  susceptibility = susceptibility_long[year==years[3]]$susceptibility, # currently not age-dependent, one for each season
  transmissibility = 0.16,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[3], '33', '1'), "%Y%W%w")
)

#### REPORTING RATES ####

gp_rate <- c(3,1,1,1,1,1,3,5,5)/50
hosp_rate <- c(8,3,1,1,1,1,5,10,20)/500
# TODO These should really vary by risk level too

care_rate_df <- data.frame(
  age_grp = age_labels,
  gp_rate = gp_rate,
  hosp_rate = hosp_rate
)

care_rate_imd_df <- cross_join(
  care_rate_df,
  data.frame(imd_quintile = 1:5)) %>% 
  mutate(gp_rate = gp_rate*(1 - (5 - imd_quintile)/20),
         hosp_rate = hosp_rate*(1 - (5 - imd_quintile)/20))

care_rate_imd_df$age_grp <- factor(care_rate_imd_df$age_grp, levels = age_labels) 

care_rate_imd_df %>% 
  pivot_longer(c(gp_rate, hosp_rate)) %>% 
  ggplot() + 
  geom_line(aes(age_grp, value, group = interaction(name,imd_quintile), 
                col = as.factor(imd_quintile), lty=name), lwd = 0.8) +
  theme_bw() + ylim(c(0,NA)) + 
  scale_color_manual(values = imd_quintile_colors) +
  labs(y='rate', lty = 'care setting', col = 'IMD quintile')

#### MAKE INTO LIST ####

unknown_pars <- list(
  epid_parameters_s1 = epid_parameters_s1,
  epid_parameters_s2 = epid_parameters_s2,
  epid_parameters_s3 = epid_parameters_s3,
  care_rates = care_rate_imd_df
)

#### SAVE UNKNOWN PARAMETERS ####

saveRDS(unknown_pars, .args[2])


