## TURN DUMMY INFECTIONS INTO SURVEILLANCE DATA ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "dummy_infections.rds"),
  file.path("data", "inputs", "dummy_surveillance.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))

set.seed(60)

#### LOAD DATA ####





