## UNKNOWN PARAMETERS ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "unknown_parameters.rds")
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

## number of years of data
years <- 2023:2025 # 2023-24 to 2025-26

#### EPIDEMIOLOGICAL PARAMETERS ####

epid_periods <- c(2, 3) # latent and infectious periods

epid_parameters_s1 <- list(
  susceptibility = 0.45, # currently not age-dependent, one for each season
  transmissibility = 0.1,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[1], '35', '1'), "%Y%W%w")
)

epid_parameters_s2 <- list(
  susceptibility = 0.5, # currently not age-dependent, one for each season
  transmissibility = 0.095,
  latent_period = epid_periods[1],
  infectious_period = epid_periods[2],
  start_date = as.Date(paste0(years[2], '36', '1'), "%Y%W%w")
)

epid_parameters_s3 <- list(
  susceptibility = 0.55, # currently not age-dependent, one for each season
  transmissibility = 0.09,
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
  ggplot() + 
  geom_line(aes(age_grp, gp_rate, group = imd_quintile, col = imd_quintile), lty=2) +
  geom_line(aes(age_grp, hosp_rate, group = imd_quintile, col = imd_quintile)) +
  theme_bw() + ylim(c(0,NA)) + ylab('')

#### MAKE INTO LIST ####

unknown_pars <- list(
  epid_parameters_s1 = epid_parameters_s1,
  epid_parameters_s2 = epid_parameters_s2,
  epid_parameters_s3 = epid_parameters_s3,
  care_rates = care_rate_imd_df
)

#### SAVE UNKNOWN PARAMETERS ####

saveRDS(unknown_pars, .args[2])


