## TURN DUMMY INFECTIONS INTO SURVEILLANCE DATA ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
suppressMessages(require(scales))
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
    left_join(unknown_pars$care_rates, by = c('age_grp','imd_quintile', 'risk_level')) %>% 
    left_join(opensafely_coverage, by = c('age_grp','imd_quintile','risk_level')) %>% 
    mutate(observed_infections = floor(OS_COVERAGE*infections)) 
  ## round to nearest integer, when considering only infections in OpenSAFELY population
  
  #### SAMPLE PRIMARY CARE ####
  
  primary_care <- c()
  for(i in 1:nrow(infections_df)){
    if(infections_df$observed_infections[i] < 1){
      primary_care <- c(primary_care, 0)
    }else{
      prob_i <- (infections_df$gp_rate[i])
      if(prob_i <= 0){prob_i <- 1e-10}
      prc_i <- rbinom(size = infections_df$observed_infections[i], n = 1,
                      prob = prob_i)
      primary_care <- c(primary_care, prc_i)
    }
  }
  
  infections_df$primary_care <- primary_care
  
  #### SAMPLE SECONDARY CARE ####
  
  secondary_care <- c()
  for(i in 1:nrow(infections_df)){
    if(infections_df$observed_infections[i] < 1){
      secondary_care <- c(secondary_care, 0)
    }else{
      prob_i <- (infections_df$hosp_rate[i])
      if(prob_i <= 0){prob_i <- 1e-10}
      sec_i <- rbinom(size = infections_df$observed_infections[i], n = 1,
                      prob = prob_i)
      secondary_care <- c(secondary_care, sec_i)
    }
  }
  
  infections_df$secondary_care <- secondary_care
  
  #### MAKE INTO TIME SERIES ####
  
  infections_df <- infections_df %>% 
    mutate(date = start_date + t)
  
  #### ADD REPORTING DELAYS ####
  key_vars <- c('date', 'age_grp', 'imd_quintile', 'risk_level', 'pop')
  key_vars_no_date <- key_vars[key_vars != 'date']
  
  infections_df <- infections_df %>% 
    group_by(!!!syms(key_vars_no_date)) %>% 
    mutate(primary_care = shift(primary_care, n = 7*known_pars$primary_care_delay, fill = 0),
           secondary_care = shift(secondary_care, n = 7*known_pars$secondary_care_delay, fill = 0))
  
  #### SELECT KEY VARIABLES ####
  
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
  geom_line(aes(date, value, col = as.factor(imd_quintile), group = imd_quintile)) +
  geom_point(aes(date, value, col = as.factor(imd_quintile), group = imd_quintile, shape = name)) +
  scale_shape_manual(values = c(1,2)) + 
  theme_bw() + labs(color='IMD') + facet_grid(name ~ imd_quintile, scales = 'free') + 
  scale_color_manual(values = imd_quintile_colors) +
  scale_x_date(breaks = "1 year", labels=date_format("%Y"))
ggsave(gsub('data','output/figures',gsub('.rds','.png',gsub('dummy_data','dummy_infections',.args[4]))), width = 16, height = 8)

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

sd_plot <- surveillance_data %>% 
  left_join(vaccinated_pop %>% select(age_grp, imd_quintile, risk_level, pop),
            by = c('age_grp','imd_quintile','risk_level')) %>% 
  left_join(opensafely_coverage,
            by = c('age_grp','imd_quintile','risk_level')) %>% 
  mutate(pop = pop*OS_COVERAGE) %>% select(!OS_COVERAGE) %>% 
  filter(index == 3) %>% 
  mutate(age_grp = case_when(
    age_grp %in% c('18-25','26-34','35-49','50-69') ~ '18-69',
    T ~ age_grp)) %>% 
  group_by(week_start, age_grp, imd_quintile) %>% 
  summarise(secondary_care = sum(secondary_care), 
            pop = sum(pop)) %>% 
  mutate(secondary_care = 100000*secondary_care/pop)

sd_plot$age_grp <- factor(sd_plot$age_grp,
                          levels = c('0-4','5-11','12-17','18-69','70-79','80+'))
sd_plot %>% 
  ggplot() +
  geom_line(aes(week_start, secondary_care, col = as.factor(imd_quintile), group = imd_quintile)) +
  geom_point(data = sd_plot %>% filter(secondary_care > 0),
             aes(week_start, secondary_care, col = as.factor(imd_quintile), group = imd_quintile),
             alpha = 1) +
  theme_bw() + labs(color='IMD') + facet_wrap(. ~ age_grp, scales = 'free') +
  scale_color_manual(values = imd_quintile_colors) +
  labs(y = 'Hospital attendance per 100,000', x='')

write_rds(surveillance_data, .args[4])

