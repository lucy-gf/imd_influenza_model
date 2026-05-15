## load population data ##

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(readxl))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "population", "imd_2025.xlsx"),
  file.path("data", "population", "lsoa_to_region.csv"),
  file.path("data", "inputs", "imd_age_pop.rds") 
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))

## load data

lsoa_imd <- data.table(read_xlsx(.args[1], sheet = 2) %>% 
                         select(starts_with('LSOA'),starts_with('Overall'),starts_with('Index')))
colnames(lsoa_imd) <- c('lsoa21cd','lsoa21nm','imd_rank','imd_decile')
lsoa_imd[, imd_quintile := ceiling(imd_decile/2)]

## load regions of England
# https://geoportal.statistics.gov.uk/datasets/ons::lsoa-2021-to-bua-to-lad-to-region-december-2022-best-fit-lookup-in-ew-v2/about
lsoa_to_reg <- data.table(read_csv(.args[2], show_col_types = F))[, c(1,11)]
colnames(lsoa_to_reg) <- c('lsoa21cd','p_engreg')

## load age-specific population in each LSOA
# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates
# using 2024 population, but can go back to 2019 with these datasets (vary `years` and `sheet_input`)
years <- c('20192022','20222024')[2]
sheet_input <- ifelse(years == '20192022', 8, 7)
lsoa_pop <- data.table(read_xlsx(file.path('data','population',paste0('sapelsoasyoa', years, '.xlsx')), 
                                 sheet = sheet_input, skip = 3))
colnames(lsoa_pop) <- c('lad21cd','lad21nm','lsoa21cd','lsoa21nm','total',paste0('F_', 0:90), paste0('M_', 0:90))

## merge:
lsoa_dat <- lsoa_pop %>% 
  select(!starts_with('lad')) %>% 
  left_join(lsoa_imd %>% select(lsoa21cd, imd_quintile), by = 'lsoa21cd') %>% 
  left_join(lsoa_to_reg, by = 'lsoa21cd') %>% 
  filter(substr(lsoa21cd, 1, 1) == 'E')

## aggregate:
imd_dat <- data.table(lsoa_dat %>% select(!starts_with('lsoa')))
imd_dat <- imd_dat[, lapply(.SD, sum), by = c('p_engreg','imd_quintile')]

## age groupings:
ages_1 <- seq(0, 75, 5)
ages_1_names <- paste0(ages_1, '-', lead(ages_1) - 1)
ages_1_names[length(ages_1_names)] <- '75+'

ages_2 <- c(0,5,12,18,26,35,50,70,80)
ages_2_names <- paste0(ages_2, '-', lead(ages_2) - 1)
ages_2_names[length(ages_2_names)] <- '80+'

## make

imd_dat_long <- imd_dat %>% 
  select(!total) %>% 
  pivot_longer(cols = !c('p_engreg', 'imd_quintile')) %>% 
  mutate(age = as.numeric(substr(name, 3, 4))) %>% 
  group_by(p_engreg, imd_quintile, age) %>% 
  summarise(pop = sum(value))

imd_age_1 <- imd_dat_long %>% 
  mutate(age_grp = cut(age, c(ages_1, Inf), right = F, labels = ages_1_names)) %>% 
  group_by(p_engreg, imd_quintile, age_grp) %>% 
  summarise(pop = sum(pop))

imd_age_2 <- imd_dat_long %>% 
  mutate(age_grp = cut(age, c(ages_2, Inf), right = F, labels = ages_2_names)) %>% 
  group_by(p_engreg, imd_quintile, age_grp) %>% 
  summarise(pop = sum(pop))


## check sums
if(sum(imd_age_1$pop) - sum(imd_age_2$pop) != 0){stop('Ages not equal')}

## save 
write_rds(imd_age_2, .args[3])
