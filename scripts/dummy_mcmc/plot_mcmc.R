## MCMC PLOTTING ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
suppressMessages(require(BayesianTools))
suppressMessages(require(patchwork))
suppressMessages(require(GGally))
suppressMessages(require(parallel))
suppressMessages(require(scales))
options(dplyr.summarise.inform = FALSE) 

.args <- #if (interactive()) c(
  c(file.path("data", "inputs", "imd_age_pop.rds"),
    file.path("data", "inputs", "contact_matrix.rds"),
    file.path("data", "inputs", "subtype_years.rds"),
    file.path("data", "dummy_data", "dummy_infections.rds"),
    file.path("data", "dummy_data", "dummy_surveillance.rds"),
    file.path("data", "dummy_data", "known_parameters.rds"),
    file.path("data", "dummy_data", "unknown_parameters.rds"),
    file.path("output", "data", "mcmc_samples.rds"),
    file.path("output", "data", "mcmc_posteriors.rds")
  ) #else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

set.seed(60)

figure_filename <- function(string){
  paste0('figures/dummy_mcmc/parameters/',string,'.png')
}

#### LOAD DATA ####

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

## ENSURE ALL POP ORDERS ARE IMD THEN AGE
imd_age_pop <- imd_age_pop %>% 
  arrange(imd_quintile, age_grp)

broad_ages <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages('children','adults','older_adults', age_labels)
)

## read in contact matrix
contact_matrix_1000 <- readRDS(.args[2])

## aggregate
contact_matrix <- contact_matrix_1000 %>% 
  group_by(p_imd_q, c_imd_q, p_age_group, c_age_group) %>% 
  summarise(n = mean(n))

contact_matrix$p_age_group <- factor(contact_matrix$p_age_group,
                                     levels = age_labels)
contact_matrix$c_age_group <- factor(contact_matrix$c_age_group,
                                     levels = age_labels)

## ENSURE ALL POP ORDERS ARE IMD THEN AGE
contact_matrix <- contact_matrix %>% 
  arrange(p_imd_q, c_imd_q,
          p_age_group, c_age_group,)

## into matrix
cm <- contact_matrix %>% ungroup() %>% 
  mutate(p_var = paste0(p_imd_q, '_', p_age_group),
         c_var = paste0(c_imd_q, '_', c_age_group)) %>% 
  select(p_var,c_var,n) %>% 
  pivot_wider(names_from = c_var, values_from = n) %>% 
  select(!p_var) %>% as.matrix()

## make matrix per capita
per_cap_matrix_45 <- t(t(cm)/imd_age_pop$pop)

## scale up for both risk groups 
## (assuming per capita contacts are independent of risk level)
pc_cm <- expand_contact_matrix(per_cap_matrix_45)

#### SUBTYPE-SEASONS #### 
ukhsa_subtype_data <- readRDS(.args[3])
subtype_seasons <- ukhsa_subtype_data %>% 
  mutate(year = as.numeric(substr(season,1,4))) %>% 
  select(season, subtype, year) %>% unique() %>% 
  mutate(epidemic = year - min(year) + 1) %>% 
  group_by(epidemic) %>% mutate(epidemic_of_season = 1:n())

#### TRUE INFECTIONS DATA #### 
true_infections_list <- readRDS(.args[4])

#### SURVEILLANCE DATA #### 
surveillance_data <- readRDS(.args[5])
if(nrow(surveillance_data) != nrow(unique(surveillance_data %>% ungroup() %>% 
                                          select(age_grp, imd_quintile, risk_level, week_start)))){
  warning('Surveillance data wrongly indexed')
}

#### KNOWN PARAMETERS #### 
known_pars <- readRDS(.args[6])
years <- known_pars$years

#### UNKNOWN PARAMETERS #### 
unknown_pars <- readRDS(.args[7])

## true epid parameters being fitted
epid_pars <- data.frame()
for(i in 1:nrow(subtype_seasons)){
  
  up_ss <- unknown_pars$epid_parameters %>% filter(
    subtype == subtype_seasons$subtype[i],
    year == subtype_seasons$year[i]
  )
  up_care_base_ss <- unknown_pars$base_care_rates %>% filter(
    subtype == subtype_seasons$subtype[i],
    year == subtype_seasons$year[i]
  )
  up_imd_spline_ss <- unknown_pars$imd_spline_pars %>% filter(
    subtype == subtype_seasons$subtype[i],
    year == subtype_seasons$year[i]
  )
  
  epid_pars <- rbind(epid_pars,
    data.frame(epidemic = subtype_seasons$epidemic[i],
               season = subtype_seasons$season[i],
               subtype = subtype_seasons$subtype[i],
               year = subtype_seasons$year[i],
               transmissibility = up_ss$transmissibility,
               children_susceptibility = up_ss$children_susceptibility,
               adults_susceptibility = up_ss$adults_susceptibility,
               older_adults_susceptibility = up_ss$older_adults_susceptibility,
               init_infected = log10(up_ss$init_infected),
               primary_care_rate_children_low = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'low'))$gp_rate,
               primary_care_rate_adults_low = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'low'))$gp_rate,
               primary_care_rate_older_adults_low = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'low'))$gp_rate,
               primary_care_rate_children_high = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'high'))$gp_rate,
               primary_care_rate_adults_high = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'high'))$gp_rate,
               primary_care_rate_older_adults_high = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'high'))$gp_rate,
               secondary_care_rate_children_low = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'low'))$hosp_rate,
               secondary_care_rate_adults_low = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'low'))$hosp_rate,
               secondary_care_rate_older_adults_low = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'low'))$hosp_rate,
               secondary_care_rate_children_high = (up_care_base_ss %>% filter(broad_age == 'children', risk_level == 'high'))$hosp_rate,
               secondary_care_rate_adults_high = (up_care_base_ss %>% filter(broad_age == 'adults', risk_level == 'high'))$hosp_rate,
               secondary_care_rate_older_adults_high = (up_care_base_ss %>% filter(broad_age == 'older_adults', risk_level == 'high'))$hosp_rate,
               imd_spline_primary_1 = up_imd_spline_ss$primary[1],
               imd_spline_primary_2 = up_imd_spline_ss$primary[2],
               imd_spline_secondary_1 = up_imd_spline_ss$secondary[1],
               imd_spline_secondary_2 = up_imd_spline_ss$secondary[2],
               R0 = R0_func(susceptibility = fcn_assign_ages(up_ss$children_susceptibility,
                                                             up_ss$adults_susceptibility,
                                                             up_ss$older_adults_susceptibility,
                                                             age_labels),
                            inf_period = known_pars$epid_periods[2],
                            beta_in = up_ss$transmissibility,
                            cm_in = cm)
    ))
}

epid_pars <- epid_pars %>% pivot_longer(!c(epidemic,subtype,season,year)) 

epid_pars <- epid_pars %>% 
  left_join(subtype_seasons, by = c('epidemic','season','subtype','year')) %>% 
  mutate(par_name = case_when(
    grepl('imd_spline', name) ~ paste0(name, '_', year),
    T ~ paste0(name, '_', year, '_', subtype)),
    joining_name =  case_when(
      grepl('imd_spline', name) ~ name,
      T ~ paste0(name, '_epid_', epidemic_of_season)))

#### MCMC SAMPLES #### 
get_samples <- function(i){
  samp <- getSample(mcmc_samples_list[[i]],coda=T,parametersOnly=F)
  chain_list <- lapply(seq_along(samp), function(i) {
    cbind(as.matrix(samp[[i]]), chain = i, iteration = 1:nrow(as.matrix(samp[[i]])))
  })
  data.table(do.call(rbind, chain_list))[, epidemic := i]
}

read_and_get_samples <- function(i){
  get_samples_parallel <- function(k){
    samp <- data.table(dat$chain[[k]])
    samp[, chain := k][, iteration := 1:nrow(samp)]
    samp
  }
  
  # if file name doesn't exist, may have finished a day earlier
  CHANGE_DATE <- !file.exists(gsub('.rds',paste0('_', i, '_', number_date_str,'.rds'),.args[8]))
  if(CHANGE_DATE){
    n_char <- nchar(number_date_str)
    number_date_str <- paste0(substr(number_date_str, 1, n_char - 1),
                              as.numeric(substr(number_date_str, n_char, n_char)) - 1)
  }
  
  dat <- readRDS(gsub('.rds',paste0('_', i, '_', number_date_str,'.rds'),.args[8]))
  list_samples <- mclapply(1:length(dat$chain), get_samples_parallel)
  samples_out <- rbindlist(list_samples)
  samples_out[, epidemic := i]
  samples_out
}

# load most recently run settings (burn in, thinning, samples)
output_details_file <- readRDS(.args[8])
number_str <- output_details_file$x[1]
run_date <- output_details_file$date
message('\n------------\nDate run: ',as.character(run_date),'\n------------',sep='')
message('\n------------\nSettings: ',number_str,'\n------------\n',sep='')
number_date_str <- paste0(number_str, '_', run_date)

# was it run on the HPC? The files saved differently
WAS_HPC <- output_details_file$HPC

if(WAS_HPC){
  
  dat_example <- readRDS(gsub('.rds',paste0('_1_', number_date_str,'.rds'),.args[8]))
  mcmc_samples <- rbindlist(lapply(1:3, read_and_get_samples))
  mcmc_samples[, c('LP','LPr') := NULL]
  
}else{
  
  mcmc_samples_list <- readRDS(gsub('.rds',paste0('_', number_date_str,'.rds'),.args[8]))
  mcmc_samples <- rbindlist(lapply(1:3, get_samples))
  mcmc_samples[, c('Lposterior','Lprior') := NULL]
  
}

total_fitted_pars <- unique(epid_pars$par_name)
total_fitted_pars <- total_fitted_pars[substr(total_fitted_pars, 1, 2) != 'R0']
n_subtype_seasons <- nrow(epid_pars %>% select(subtype, season) %>% unique())
n_seasons <- n_distinct(epid_pars$season)
message('\n', length(total_fitted_pars),' fitted parameters, ', 
        n_subtype_seasons,' unique epidemics, over ', n_seasons,
        ' seasons.\n', sep = '')

epids_actual <- n_distinct(mcmc_samples$epidemic)
if(3 %in% unique(mcmc_samples$epidemic)){epids_actual <- epids_actual - 17/38}
total_actual_fitted_pars <- (ncol(mcmc_samples) - 4)*epids_actual
message('\n', total_actual_fitted_pars,' posteriors in the mcmc_samples file.\n', sep = '')

fitted_pars <- unique(epid_pars$name)
fitted_pars <- fitted_pars[substr(fitted_pars, 1, 2) != 'R0']
fitted_pars_not_imd <- fitted_pars[!grepl('imd_spline', fitted_pars)]
fitted_pars_imd <- fitted_pars[grepl('imd_spline', fitted_pars)]
posterior_cols <- c(paste0(fitted_pars_not_imd, '_epid_1'),
                    paste0(fitted_pars_not_imd, '_epid_2'), 
                    fitted_pars_imd)
admin_cols <- c('likelihood', 'chain', 'iteration', 'epidemic')
colnames(mcmc_samples) <- c(posterior_cols,
                            admin_cols)

#### ADD R0 #### 
cat('Adding R0: ')
unique_df <- unique(mcmc_samples[, ..posterior_cols])

pb <- txtProgressBar(min = 1, max = nrow(unique_df), style = 3)

for(i in 1:nrow(unique_df)){
  
  row <- unique_df[i, ]
  
  mcmc_samples[transmissibility_epid_1 == row$transmissibility_epid_1 &
                 children_susceptibility_epid_1 == row$children_susceptibility_epid_1 & 
                 adults_susceptibility_epid_1 == row$adults_susceptibility_epid_1 &
                 older_adults_susceptibility_epid_1 == row$older_adults_susceptibility_epid_1,
               R0_epid_1 := R0_func(susceptibility = fcn_assign_ages(row$children_susceptibility_epid_1,
                                                                     row$adults_susceptibility_epid_1,
                                                                     row$older_adults_susceptibility_epid_1,
                                                                     age_labels),
                             inf_period = known_pars$epid_periods[2],
                             beta_in = row$transmissibility_epid_1,
                             cm_in = cm)
               ]
  mcmc_samples[transmissibility_epid_2 == row$transmissibility_epid_2 &
                 children_susceptibility_epid_2 == row$children_susceptibility_epid_2 & 
                 adults_susceptibility_epid_2 == row$adults_susceptibility_epid_2 &
                 older_adults_susceptibility_epid_2 == row$older_adults_susceptibility_epid_2,
               R0_epid_2 := R0_func(susceptibility = fcn_assign_ages(row$children_susceptibility_epid_2,
                                                                     row$adults_susceptibility_epid_2,
                                                                     row$older_adults_susceptibility_epid_2,
                                                                     age_labels),
                                    inf_period = known_pars$epid_periods[2],
                                    beta_in = row$transmissibility_epid_2,
                                    cm_in = cm)
  ]
  
  # Print progress
  setTxtProgressBar(pb, i)
  
}
close(pb)

posterior_cols <- c(posterior_cols, 'R0_epid_1', 'R0_epid_2')
plotting_cols <- unique(gsub('_epid_1|_epid_2', '', posterior_cols))

#### FILTER #### 
burn_in <- as.numeric(strsplit(number_str, split = '_')[[1]][1])
thinning_value <- as.numeric(strsplit(number_str, split = '_')[[1]][2])
n_samples <- (max(mcmc_samples$iteration) - burn_in)/thinning_value

mcmc_samples_filtered <- mcmc_samples[iteration > burn_in & iteration %% thinning_value == 0,]
mcmc_samples_filtered[, iteration := 1:n_samples, .(chain, epidemic)]

#### PLOT DENSITY #### 
densities <- map(.x = plotting_cols, .f = plot_density)
patchwork::wrap_plots(densities, nrow = 6)
ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('fitted_densities'),.args[length(.args)]), width = 30, height = 14)

#### PLOT TRACE #### 
traces <- map(.x = plotting_cols, .f = plot_trace)
patchwork::wrap_plots(traces, nrow = 6) 
ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('fitted_traces'),.args[length(.args)]), width = 30, height = 14)
traces_filtered <- map(.x = plotting_cols, .f = ~{plot_trace(var=.x, filtered=T)})
patchwork::wrap_plots(traces_filtered, nrow = 6) 
ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('fitted_filtered_traces'),.args[length(.args)]), width = 30, height = 14)
plot_trace(var="transmissibility", filtered=T)
ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('example_filtered_trace'),.args[length(.args)]), width = 10, height = 5)
log_likelihood_plot <- plot_trace('likelihood'); log_likelihood_plot
ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('fitted_likelihood'),.args[length(.args)]), width = 8, height = 6)

## PAIRWISE PLOTS

mcmc_samples_filtered[, chain_it_id := paste0(chain, '_', iteration)]

pairs_cols <- c(posterior_cols, 'epidemic', 'chain_it_id')

pairs_data <- if(nrow(mcmc_samples_filtered) >= 10000){
  # only taking 1% of the mcmc_samples_filtered dataset as it takes too long otherwise!
  mcmc_samples_filtered[seq(1, nrow(mcmc_samples_filtered), by = 100), ..pairs_cols]
}else{
  mcmc_samples_filtered[1:nrow(mcmc_samples_filtered), ..pairs_cols]
}

pairs_long <- melt(pairs_data, id.vars = c('epidemic', 'chain_it_id'))
pairs_long[, variable := as.character(variable)]
pairs_long[, epidemic_of_season := as.numeric(substr(variable, nchar(variable), nchar(variable)))]
pairs_long[, variable := substr(variable, 1, nchar(variable) - 7)]

pairs_long <- pairs_long[subtype_seasons, on = c('epidemic','epidemic_of_season')]
pairs_long[, subtype_season := paste0(season, ': ', subtype)]
pairs_long[, c('epidemic','epidemic_of_season','subtype','season', 'year') := NULL]

pairs_wide <- dcast(pairs_long, 
                    chain_it_id + subtype_season ~ variable, value.var = 'value')

colnames(pairs_wide) <- gsub('_rate_','_rate\n_', colnames(pairs_wide))
colnames(pairs_wide) <- gsub('_spline_','_spline\n_', colnames(pairs_wide))

pairs_wide[,chain_it_id := NULL]

p_FULL <- ggpairs(pairs_wide, columns = 2:ncol(pairs_wide), aes(color = as.factor(subtype_season), alpha = 0.5))
ggsave(filename = gsub('data/mcmc_posteriors.rds',figure_filename('fitted_pairwise'),.args[length(.args)]),
       plot = p_FULL, width = 40, height = 40)

#### SAVE DATA ####

write_rds(mcmc_samples_filtered, .args[length(.args)])






