## TURN DUMMY INFECTIONS INTO SURVEILLANCE DATA ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "dummy_infections.rds"),
  file.path("data", "dummy_data", "known_parameters.rds"),
  file.path("data", "dummy_data", "unknown_parameters.rds"),
  file.path("data", "dummy_data", "dummy_surveillance.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))

set.seed(60)

#### LOAD DATA ####

infections <- readRDS(.args[1])

## KNOWN PARAMETERS
known_pars <- readRDS(.args[2])
years <- known_pars$years

risk_group_pop <- known_pars$risk_group_pop
vaccinated_pop <- known_pars$vaccinated_pop

## UNKNOWN PARAMETERS
unknown_pars <- readRDS(.args[3])

#### SAMPLE SURVEILLANCE DATA ####

full_df <- data.frame()

for(k in 1:length(infections)){
  
  if(k==1){cat('Year: ')}
  
  year_i <- names(infections)[k]
  
  cat(year_i, ', ', sep = '')
  
  ## turn into new infections (from compartment cumI)
  infections_df <- infections[[k]] %>% 
    filter(compartment == 'cumI') %>% 
    group_by(age_grp, imd_quintile, risk_level, pop, start_date) %>% 
    mutate(infections = value - lag(value, default = 0)) %>% ungroup() %>% 
    mutate(imd_quintile = as.numeric(imd_quintile),
           infections = round(infections)) %>% ## round to nearest integer
    select(!c(compartment,value)) %>% 
    left_join(unknown_pars$care_rates, by = c('age_grp','imd_quintile'))
  
  #### SAMPLE PRIMARY CARE ####
  
  jitter_width <- 0.01
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)
  
  primary_care <- c()
  for(i in 1:nrow(infections_df)){
    if(infections_df$infections[i] < 1){
      primary_care <- c(primary_care, 0)
    }else{
      prob_i <- (infections_df$gp_rate[i] + unif_jitter[i])
      if(prob_i <= 0){prob_i <- 1e-10}
      prc_i <- rbinom(size = infections_df$infections[i], n = 1,
                      prob = prob_i)
      primary_care <- c(primary_care, prc_i)
    }
  }
  primary_care <- c(rep(0, known_pars$primary_care_delay), primary_care)
  primary_care <- primary_care[1:nrow(infections_df)]
  
  infections_df$primary_care <- primary_care
  
  #### SAMPLE SECONDARY CARE ####
  
  jitter_width <- 0.001
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)
  
  secondary_care <- c()
  for(i in 1:nrow(infections_df)){
    if(infections_df$infections[i] < 1){
      secondary_care <- c(secondary_care, 0)
    }else{
      prob_i <- (infections_df$hosp_rate[i] + unif_jitter[i])
      if(prob_i <= 0){prob_i <- 1e-10}
      sec_i <- rbinom(size = infections_df$infections[i], n = 1,
                      prob = prob_i)
      secondary_care <- c(secondary_care, sec_i)
    }
  }
  secondary_care <- c(rep(0, known_pars$secondary_care_delay), secondary_care)
  secondary_care <- secondary_care[1:nrow(infections_df)]
  
  infections_df$secondary_care <- secondary_care
  
  #### MAKE INTO TIME SERIES ####
  
  infections_df <- infections_df %>% 
    mutate(date = start_date + t)
  
  #### SELECT KEY VARIABLES ####
  
  key_vars <- c('date', 'age_grp', 'imd_quintile', 'risk_level', 'pop')
  
  infections_filtered <- infections_df %>% 
    select(!!!syms(key_vars), infections, primary_care, secondary_care) %>% 
    mutate(index = k)
  
  full_df <- rbind(full_df,
                   infections_filtered)
  
}

full_df_agg <- full_df %>% 
  group_by(!!!syms(key_vars), index) %>% 
  summarise(infections = sum(infections), 
            primary_care = sum(primary_care), 
            secondary_care = sum(secondary_care))
 
full_df_agg %>% group_by(date, imd_quintile, risk_level) %>%
  summarise(infections = sum(infections),
            primary_care = sum(primary_care),
            secondary_care = sum(secondary_care)) %>% 
  filter(year(date)=='2024') %>% 
  ggplot() +
  geom_line(aes(date, primary_care, col = as.factor(imd_quintile), group = imd_quintile)) +
  geom_line(aes(date, secondary_care, col = as.factor(imd_quintile), group = imd_quintile),lty=2) +
  facet_grid(risk_level ~ ., scales = 'free') +
  theme_bw() + labs(color='IMD') +
  scale_color_manual(values = imd_quintile_colors)

full_df_agg %>% group_by(date, imd_quintile) %>%
  summarise(infections = sum(infections),
            primary_care = sum(primary_care),
            secondary_care = sum(secondary_care)) %>% 
  filter(year(date)=='2024') %>% 
  ggplot() +
  geom_line(aes(date, primary_care, col = as.factor(imd_quintile), group = imd_quintile)) +
  geom_line(aes(date, secondary_care, col = as.factor(imd_quintile), group = imd_quintile),lty=2) +
  theme_bw() + labs(color='IMD') +
  scale_color_manual(values = imd_quintile_colors)

## print outcomes

cat('Outcomes per season:\n')

full_df_agg %>% group_by(index) %>%
  summarise(infections = sum(infections),
            primary_care = sum(primary_care),
            secondary_care = sum(secondary_care)) 

#### SAVE DATA ####

write_rds(full_df_agg, .args[4])
