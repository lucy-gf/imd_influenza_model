## TURN DUMMY INFECTIONS INTO SURVEILLANCE DATA ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "dummy_data", "dummy_infections.rds"),
  file.path("data", "dummy_data", "known_parameters.rds"),
  file.path("data", "dummy_data", "unknown_parameters.rds"),
  file.path("data", "dummy_data", "dummy_surveillance.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','seir_model.R'))

set.seed(60)

#### LOAD DATA ####

infections <- readRDS(.args[1])

## KNOWN PARAMETERS
known_pars <- readRDS(.args[2])
years <- known_pars$years

risk_group_pop <- known_pars$risk_group_pop
vaccinated_pop <- known_pars$vaccinated_pop

opensafely_coverage <- known_pars$proportion_observed

## UNKNOWN PARAMETERS
unknown_pars <- readRDS(.args[3])

#### SAMPLE SURVEILLANCE DATA ####

full_df <- data.frame()

for(k in 1:length(infections)){
  
  if(k==1){cat('Year: ')}
  
  year_i <- names(infections)[k]
  
  cat(year_i, ', ', sep = '')
  
  ## merge with coverage rates and attendance rates
  infections_df <- infections[[k]] %>% 
    mutate(imd_quintile = as.numeric(imd_quintile)) %>% 
    left_join(unknown_pars$care_rates, by = c('age_grp','imd_quintile')) %>% 
    left_join(opensafely_coverage, by = c('age_grp','imd_quintile','risk_level')) %>% 
    mutate(observed_infections = round(OS_COVERAGE*infections)) 
  ## round to nearest integer, when considering only infections in OpenSAFELY population
  
  #### SAMPLE PRIMARY CARE ####
  
  jitter_width <- 0 # 0.01 # no jitter for now
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)
  
  primary_care <- c()
  for(i in 1:nrow(infections_df)){
    if(infections_df$observed_infections[i] < 1){
      primary_care <- c(primary_care, 0)
    }else{
      prob_i <- (infections_df$gp_rate[i] + unif_jitter[i])
      if(prob_i <= 0){prob_i <- 1e-10}
      prc_i <- rbinom(size = infections_df$observed_infections[i], n = 1,
                      prob = prob_i)
      primary_care <- c(primary_care, prc_i)
    }
  }
  primary_care <- c(rep(0, 7*known_pars$primary_care_delay), primary_care)
  primary_care <- primary_care[1:nrow(infections_df)]
  
  infections_df$primary_care <- primary_care
  
  #### SAMPLE SECONDARY CARE ####
  
  jitter_width <- 0 # 0.001 # no jitter for now
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)
  
  secondary_care <- c()
  for(i in 1:nrow(infections_df)){
    if(infections_df$observed_infections[i] < 1){
      secondary_care <- c(secondary_care, 0)
    }else{
      prob_i <- (infections_df$hosp_rate[i] + unif_jitter[i])
      if(prob_i <= 0){prob_i <- 1e-10}
      sec_i <- rbinom(size = infections_df$observed_infections[i], n = 1,
                      prob = prob_i)
      secondary_care <- c(secondary_care, sec_i)
    }
  }
  secondary_care <- c(rep(0, 7*known_pars$secondary_care_delay), secondary_care)
  secondary_care <- secondary_care[1:nrow(infections_df)]
  
  infections_df$secondary_care <- secondary_care
  
  #### MAKE INTO TIME SERIES ####
  
  infections_df <- infections_df %>% 
    mutate(date = start_date + t)
  
  #### SELECT KEY VARIABLES ####
  
  key_vars <- c('date', 'age_grp', 'imd_quintile', 'risk_level', 'pop')
  
  infections_filtered <- infections_df %>% 
    select(!!!syms(key_vars), infections, observed_infections, primary_care, secondary_care) %>% 
    mutate(index = k)
  
  full_df <- rbind(full_df,
                   infections_filtered)
  
}

full_df_agg <- full_df %>% 
  group_by(!!!syms(key_vars), index) %>% 
  summarise(infections = sum(infections),
            observed_infections = sum(observed_infections),
            primary_care = sum(primary_care), 
            secondary_care = sum(secondary_care))

full_df_agg %>% mutate(date = last_monday(date)) %>% 
  group_by(date, imd_quintile) %>%
  summarise(primary_care = sum(primary_care),
            secondary_care = sum(secondary_care)) %>% 
  pivot_longer(c(primary_care,secondary_care)) %>% 
  ggplot() +
  geom_point(aes(date, value, col = as.factor(imd_quintile), group = imd_quintile, shape = name)) +
  scale_shape_manual(values = c(1,2)) + 
  theme_bw() + labs(color='IMD') + facet_grid(name ~ imd_quintile, scales = 'free') + 
  scale_color_manual(values = imd_quintile_colors)

## print outcomes

cat('Outcomes per season:\n')

full_df_agg %>% group_by(index) %>%
  summarise(infections = sum(infections),
            observed_infections = sum(observed_infections),
            proportion_observed = paste0(round(100*sum(observed_infections)/sum(infections), 2), '%'),
            primary_care = sum(primary_care),
            secondary_care = sum(secondary_care)) 

#### SAVE DATA ####

# actual surveillance data won't include infections, 
# and will only be weekly
surveillance_data <- full_df_agg %>% 
  group_by(age_grp, imd_quintile, risk_level) %>% 
  complete(date = seq.Date(from = as.Date(paste0('01-01-',year(full_df_agg$date[1])), format = '%d-%m-%Y') + 7,
                           to = max(full_df_agg$date), by = 1),
           fill = list(index = 1, infections = 0, primary_care = 0, secondary_care = 0)) %>% 
  mutate(week_start = last_monday(date)) %>% 
  group_by(week_start, age_grp, imd_quintile, index, risk_level) %>% 
  summarise(primary_care = sum(primary_care),
            secondary_care = sum(secondary_care))

surveillance_data %>% 
  ggplot() +
  geom_line(aes(week_start, primary_care, col = as.factor(imd_quintile), group = imd_quintile)) +
  geom_line(aes(week_start, secondary_care, col = as.factor(imd_quintile), group = imd_quintile),lty=2) +
  theme_bw() + labs(color='IMD') + facet_grid(risk_level ~ age_grp, scales = 'free') +
  scale_color_manual(values = imd_quintile_colors)

write_rds(surveillance_data, .args[4])

