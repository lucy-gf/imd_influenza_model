## CREATE DUMMY DATA FOR MODEL FITTING ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "dummy_flu_data.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))

set.seed(60)

#### LOAD DATA ####

## number of years of data
years <- 2023:2026 # 2023-24 to 2025-26

## read in population data

imd_age_pop <- readRDS(.args[1])
age_labels <- unique(imd_age_pop$age_grp)

## rough epidemiological parameters

attack_rate <- c(3,4,5,4,3,2,2,1,1)/10
gp_rate <- c(3,1,1,1,1,1,3,5,5)/10
hosp_rate <- c(8,3,1,1,1,1,5,10,20)/100

epi_pars <- data.frame(
  age_grp = age_labels,
  attack_rate = attack_rate,
  gp_rate = gp_rate,
  hosp_rate = hosp_rate
)

epi_pars %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, y = attack_rate, group=1), col = 1) +
  geom_line(aes(x = age_grp, y = attack_rate*gp_rate, group=1), col = 2) +
  geom_line(aes(x = age_grp, y = attack_rate*hosp_rate, group=1), col = 3) +
  theme_bw() + ylim(c(0,NA))

#### MAKE DUMMY DATA #### 

dummy_data <- cross_join(CJ(year = years,week = 1:52),
                         imd_age_pop %>% select(!pop))

dummy_data <- dummy_data %>% 
  mutate(infections = 0, primary_care = 0, secondary_care = 0)

for(year_i in years[1:(length(years)-1)]){
  
  cat('Year: ', year_i, ', ', sep = '')
  
  #### SAMPLE INFECTIONS #### 
  
  infections_df <- imd_age_pop %>% 
    left_join(epi_pars, by = 'age_grp') 
  
  jitter_width <- 0.03
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)
  
  infections <- c()
  for(i in 1:nrow(infections_df)){
    inf_i <- rbinom(size = infections_df$pop[i], n = 1,
                    prob = (infections_df$attack_rate[i] + 
                              unif_jitter[i] - infections_df$imd_quintile[i]/100))
    infections <- c(infections, inf_i)
  }
  
  infections_df$total_infections <- infections
  
  #### SAMPLE PRIMARY CARE ####

  jitter_width <- 0.01
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)

  primary_care <- c()
  for(i in 1:nrow(infections_df)){
    prc_i <- rbinom(size = infections_df$total_infections[i], n = 1,
                    prob = (infections_df$gp_rate[i] +
                              unif_jitter[i] - infections_df$imd_quintile[i]/100))
    primary_care <- c(primary_care, prc_i)
  }

  infections_df$primary_care_rate <- primary_care/infections_df$total_infections

  #### SAMPLE SECONDARY CARE ####

  jitter_width <- 0.001
  unif_jitter <- runif(nrow(infections_df), min = -jitter_width, max = jitter_width)

  secondary_care <- c()
  for(i in 1:nrow(infections_df)){
    sec_i <- rbinom(size = infections_df$total_infections[i], n = 1,
                    prob = (infections_df$hosp_rate[i] +
                              unif_jitter[i]))
    secondary_care <- c(secondary_care, sec_i)
  }

  infections_df$secondary_care_rate <- secondary_care/infections_df$total_infections
  
  #### MAKE INTO TIME SERIES ####
  
  ## baseline weekly incidence
  n_week <- sample(x = 35:45, size = 1) # number of weeks of epidemic, must be <= 52
  x <- 1:52
  weekly_time_series <- (x^5)*((n_week-x)^2)
  weekly_time_series[(n_week + 1):52] <- 0 # set to 0 outside of epidemic
  weekly_time_series <- weekly_time_series/sum(weekly_time_series) # normalise
  
  epidemic_start <- sample(25:35, size = 1) # week of the year of start
  
  weekly_time_series <- c(rep(0,epidemic_start - 1), weekly_time_series, rep(0,53 - epidemic_start))
  
  # time delays for primary (1wk) and secondary care (2wk)
  prim_delay <- 1; sec_delay <- 2
  primary_time_series <- c(rep(0, prim_delay), weekly_time_series[1:(length(weekly_time_series) - prim_delay)])
  secondary_time_series <- c(rep(0, sec_delay), weekly_time_series[1:(length(weekly_time_series) - sec_delay)])
  
  ## scale weekly time series by attack rates etc.
  join_vec <- c('p_engreg','imd_quintile','age_grp')
  
  dummy_data_filt <- dummy_data %>% 
    filter(year %in% year_i:(year_i + 1)) %>% 
    left_join(infections_df %>% 
                select(!!!syms(join_vec), 
                       total_infections, primary_care_rate, secondary_care_rate), 
              by = join_vec) %>% 
    arrange(!!!syms(join_vec)) %>% 
    mutate(inf_ts = rep(weekly_time_series, nrow(dummy_data)/(length(years)*52)),
           weekly_inf = 0, weekly_primary = 0, weekly_secondary = 0)
  
  for(k in 1:nrow(dummy_data_filt)){
    if(dummy_data_filt$inf_ts[k]==0){
      dummy_data_filt$weekly_inf[k] <- 0
      dummy_data_filt$weekly_primary[k] <- 0
      dummy_data_filt$weekly_secondary[k] <- 0
    }else{
      jitters <- runif(3, 0.8, 1)
      sample_inf <- rbinom(n = 1, size = dummy_data_filt$total_infections[k], 
                           prob = jitters[1]*dummy_data_filt$inf_ts[k])
      sample_prim <- rbinom(n = 1, size = sample_inf, 
                            prob = jitters[2]*dummy_data_filt$primary_care_rate[k])
      sample_sec <- rbinom(n = 1, size = sample_inf, 
                           prob = jitters[3]*dummy_data_filt$secondary_care_rate[k])
      
      if(sample_inf < sample_prim){warning(cat('Primary > Infections (k = ', k, ')', sep = ''))}
      if(sample_inf < sample_sec){warning(cat('Secondary > Infections (k = ', k, ')', sep = ''))}
      
      dummy_data_filt$weekly_inf[k] <- sample_inf
      dummy_data_filt$weekly_primary[k] <- sample_prim
      dummy_data_filt$weekly_secondary[k] <- sample_sec
    }
    
    if(k == 1){cat('Rows done: ', sep='')}
    if(k %% 5000 == 0){cat(k, ', ', sep='')}
    if(k == nrow(dummy_data_filt)){cat(k,'\n', sep='')}
  }
  
  dummy_data <- dummy_data %>% 
    left_join(dummy_data_filt %>% 
                select(!!!syms(join_vec), year, week,
                       weekly_inf, weekly_primary, weekly_secondary),
              by = c(join_vec, 'year','week')) %>% 
    mutate(across(c(weekly_inf, weekly_primary, weekly_secondary), na_to_0)) %>% 
    mutate(infections = infections + weekly_inf,
           primary_care = primary_care + weekly_primary,
           secondary_care = secondary_care + weekly_secondary) %>% 
    select(!c(weekly_inf, weekly_primary, weekly_secondary))
  
}

## add week as date
dummy_data <- dummy_data %>% 
  mutate(week_char = case_when(week < 10 ~ paste0('0',week), T ~ paste0(week)),
         date = as.Date(paste0(year,week_char,'1'), "%Y%W%w"))

## add population back in 
dummy_data <- dummy_data %>% 
  left_join(imd_age_pop,
            by = c('imd_quintile','age_grp','p_engreg'))

plot_rates_func <- function(char){
  
  dummy_data %>% 
    filter(age_grp==char) %>% 
    pivot_longer(!c(year,week,p_engreg,imd_quintile,age_grp,week_char,date,pop)) %>% 
    ggplot() + 
    geom_line(aes(x = date, y = value/pop, col = as.factor(imd_quintile), 
                  group = interaction(as.factor(imd_quintile),name), lty = name)) + 
    facet_wrap(p_engreg ~ .) + labs(col = 'IMD', lty = 'Type', y = 'Rate') + 
    scale_color_manual(values = imd_quintile_colors) + 
    theme_bw() +
    ggtitle(paste0('Infections, primary care, and secondary care (age group ', char, ')'))

  }

plot_rates_func(age_labels[5])

dummy_data %>% 
  group_by(p_engreg, age_grp, imd_quintile, pop) %>% 
  summarise(infections = sum(infections),
            primary_care = sum(primary_care),
            secondary_care = sum(secondary_care)) %>% 
  ggplot() + 
  geom_bar(aes(x = age_grp, y = secondary_care/(pop*3), fill = as.factor(imd_quintile), 
                group = as.factor(imd_quintile)),
           position = 'dodge', stat='identity') + 
  facet_wrap(p_engreg ~ .) + labs(fill = 'IMD') + 
  scale_fill_manual(values = imd_quintile_colors) + 
  theme_bw() + ggtitle('Mean annual secondary care by region, age, and IMD quintile')

dummy_data %>% 
  pivot_longer(!c(year,week,p_engreg,imd_quintile,age_grp,week_char,date,pop)) %>% 
  group_by(date, imd_quintile, age_grp, name) %>% 
  summarise(value = sum(value), pop = sum(pop)) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = value/pop, col = as.factor(imd_quintile), 
                group = interaction(as.factor(imd_quintile),name), lty = name)) + 
  facet_wrap(age_grp ~ ., scales = 'free') + labs(col = 'IMD', lty = 'Type', y = 'Rate') + 
  scale_color_manual(values = imd_quintile_colors) + 
  theme_bw() 

dummy_data %>% 
  pivot_longer(!c(year,week,p_engreg,imd_quintile,age_grp,week_char,date,pop)) %>% 
  group_by(date, imd_quintile, age_grp, name) %>% 
  summarise(value = sum(value), pop = sum(pop)) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = value, col = as.factor(imd_quintile), 
                group = interaction(as.factor(imd_quintile),name), lty = name)) + 
  facet_wrap(age_grp ~ ., scales = 'free') + labs(col = 'IMD', lty = 'Type', y = 'Rate') + 
  scale_color_manual(values = imd_quintile_colors) + 
  theme_bw() 


#### TURN INTO LINELIST #### 





#### SAVE #### 




