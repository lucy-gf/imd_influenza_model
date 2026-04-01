## loading in IMD- and age-stratified contact matrices ##

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "contact_matrix", "fitted_matrs_balanced.csv"),
  file.path("data", "inputs", "contact_matrix.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))

## in reality this should have a national/regional sensitivity analysis flag

## read
cm <- read_csv(.args[1])

## summarise
group_vars <- c('p_age_group','c_age_group','p_imd_q','c_imd_q')

summ_matr <- cm %>% 
  group_by(!!!syms(group_vars)) %>% 
  summarise(n = mean(n))

## save
write_rds(summ_matr, .args[2])

