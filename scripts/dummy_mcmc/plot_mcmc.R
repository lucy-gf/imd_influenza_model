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
surveillance_data <- readRDS(.args[5])[[1]]
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
epid_pars <- readRDS(.args[8])

#### MCMC SAMPLES #### 
get_samples <- function(i){
  samp <- getSample(mcmc_samples_list[[i]],coda=T,parametersOnly=F)
  chain_list <- lapply(seq_along(samp), function(i) {
    cbind(as.matrix(samp[[i]]), chain = i, iteration = 1:nrow(as.matrix(samp[[i]])))
  })
  data.table(do.call(rbind, chain_list))[, epidemic := i]
}

read_and_get_samples <- function(i,
                                 run_date_in = run_date,
                                 number_str_in = number_str){
  
  get_samples_parallel <- function(k){
    samp <- data.table(dat$chain[[k]])
    llcomps <- dat$LLcomponents[[k]][, c('LL1','LL2','LL3')]
    samp <- cbind(samp, llcomps)
    samp[, job := chain_job]
    samp[, chain := k + length(dat$chain)*(chain_job - 1)]
    samp[, iteration := 1:nrow(samp)]
    samp
  }
  
  extract_data <- function(k){
    
    samp <- data.table(dat[[k]])
    
    # if in 2025/26, add blank parameters
    if(ncol(samp) == 27){
      samp <- cbind(samp[, 1:17], 
                    matrix(0, ncol = 17, nrow = nrow(samp)),
                    samp[, 18:27])
      colnames(samp) <- c(paste0('par ', 1:38), colnames(samp)[39:44])
    }
    
    samp[, job := chain_job]
    samp[, chain := k + length(dat)*(chain_job - 1)]
    samp[, iteration := 1:nrow(samp)]
    samp
    
  }
  
  # Code could have finished running in the days before
  dates <- as.character(run_date_in - 0:6)
  
  # files in directory
  all_files <- list.files(path = gsub('mcmc_samples.rds',number_str_in,.args[9]))
  
  # files finishing on relevant dates
  date_files <- c()
  for(date in dates){date_files <- c(date_files, all_files[grepl(date, all_files)])}
  
  # files with correct burn-in, thinning etc.
  date_files <- date_files[grepl(number_str_in, date_files)]
  
  # files with correct index
  date_files <- date_files[grepl(paste0('INDEX',i), date_files)]
  
  # how many chain jobs for a given index
  chain_jobs <- suppressWarnings(unique(c(as.numeric(substr(date_files, 26, 26)),
                                          as.numeric(substr(date_files, 26, 27)))))
  if(sum(grepl('CHAIN1_', date_files)) == 0){chain_jobs <- chain_jobs[chain_jobs != 1]}
  chain_jobs <- chain_jobs[!is.na(chain_jobs)]
  total_chains <- length(chain_jobs)
  
  samples_out <- data.table()
  
  for(chain_job in chain_jobs){
    
    date_file <- date_files[grepl(paste0("CHAIN", chain_job, "_"), date_files)]
    if(length(date_file) != 1){stop(paste0('Too many eligible data files (', chain_job, ': ', date_file[1], ')'))}
    
    dat <- tryCatch(expr = suppressWarnings(readRDS(gsub('mcmc_samples.rds',paste0(number_str_in,'/',date_file),.args[9]))),
                    error = NA)
    
    if(length(dat) == 1){
      if(is.na(dat)){
        dat <- tryCatch(expr = suppressWarnings(readRDS(gsub('mcmc_samples.rds',paste0(number_str_in,'/',date_file),.args[9]))))
      }
    }
    
    list_samples <- if(save_all_flag){
      mclapply(1:length(dat), get_samples_parallel)
    }else{
        mclapply(1:length(dat), extract_data)
      }
    
    samples_out_CJ <- rbindlist(list_samples)
    
    samples_out <- rbind(samples_out, 
                         samples_out_CJ)
    
  }
  
  samples_out[, epidemic := i]
  
  return(samples_out)
}

# load most recently run settings (burn in, thinning, samples)
output_details_file <- readRDS(.args[9])
number_str <- output_details_file$x[1]
run_date <- output_details_file$date
save_all_flag <- output_details_file$save_all
merge_flag <- output_details_file$using_previous_pars
if(merge_flag){
  merging_run <- output_details_file$previous_pars
}

message('\n------------\nDate run: ',as.character(run_date),'\n------------',sep='')
message('\n------------\nSettings: ',number_str,'\n------------\n',sep='')
message('\n------------\nSaved: ',ifelse(save_all_flag,'all','minimal'),'\n------------\n',sep='')
if(merge_flag){message('\n------------\nMerge with previous run (', merging_run,')\n------------\n',sep='')}

# was it run on the HPC? The files saved differently
WAS_HPC <- output_details_file$HPC

if(WAS_HPC){
  
  m1 <- read_and_get_samples(1)
  m2 <- read_and_get_samples(2)
  m3 <- read_and_get_samples(3)
  
  if(merge_flag){ # merge with previous iterations, bump up iteration count
    
    m1_prev <- read_and_get_samples(1,
                                    output_details_file$previous_pars_date,
                                    merging_run)
    m2_prev <- read_and_get_samples(2,
                                    output_details_file$previous_pars_date,
                                    merging_run)
    m3_prev <- read_and_get_samples(3,
                                    output_details_file$previous_pars_date,
                                    merging_run)
    
    m1 <- rbind(m1_prev, m1[, iteration := iteration + max(m1_prev$iteration)])
    m2 <- rbind(m2_prev, m2[, iteration := iteration + max(m2_prev$iteration)])
    m3 <- rbind(m3_prev, m3[, iteration := iteration + max(m3_prev$iteration)])
    
  }
  
  # bind any non-empty data
  mcmc_samples <- data.table()
  for(t in 1:3){
    dt <- get(paste0('m',t))
    if(ncol(dt) > 1){
      mcmc_samples <- rbind(mcmc_samples, dt)
    }
  }
  
  mcmc_samples[, c('LP','LPr') := NULL]
  
}else{
  
  number_date_str <- paste0(number_str, '_', run_date)
  
  mcmc_samples_list <- readRDS(gsub('.rds',paste0('_', number_date_str,'.rds'),.args[9]))
  mcmc_samples <- rbindlist(lapply(1:3, get_samples))
  mcmc_samples[, c('Lposterior','Lprior') := NULL]
  
}

total_fitted_pars <- unique(epid_pars$par_name)
total_fitted_pars <- total_fitted_pars[substr(total_fitted_pars, 1, 4) != 'Reff']
n_subtype_seasons <- nrow(epid_pars %>% select(subtype, season) %>% unique())
n_seasons <- n_distinct(epid_pars$season)
message('\n', length(total_fitted_pars),' fitted parameters, ', 
        n_subtype_seasons,' unique epidemics, over ', n_seasons,
        ' seasons.\n', sep = '')

epids_actual <- n_distinct(mcmc_samples$epidemic)
if(3 %in% unique(mcmc_samples$epidemic)){epids_actual <- epids_actual - 17/38}
total_actual_fitted_pars <- (ncol(mcmc_samples) - 8)*epids_actual
message('\n', total_actual_fitted_pars,' posteriors in the mcmc_samples file.\n', sep = '')

message('\n', n_distinct(mcmc_samples$epidemic),' seasons in the mcmc_samples file.\n', sep = '')
message('\n', mcmc_samples %>% group_by(epidemic) %>% 
                     summarise(chains = n_distinct(chain)) %>% select(chains),' chains per season.\n', sep = '')

fitted_pars <- unique(epid_pars$name)
fitted_pars <- fitted_pars[substr(fitted_pars, 1, 4) != 'Reff']
fitted_pars_not_imd <- fitted_pars[!grepl('imd_spline', fitted_pars)]
fitted_pars_imd <- fitted_pars[grepl('imd_spline', fitted_pars)]
posterior_cols <- c(paste0(fitted_pars_not_imd, '_epid_1'),
                    paste0(fitted_pars_not_imd, '_epid_2'), 
                    fitted_pars_imd)

admin_cols <- c('likelihood',
                'gp_likelihood', 'hospital_likelihood', 'subtype_likelihood',
                'job', 'chain', 'iteration', 'epidemic')
colnames(mcmc_samples) <- c(posterior_cols,
                            admin_cols)

#### FILTER #### 
burn_in <- as.numeric(strsplit(number_str, split = '_')[[1]][1])
thinning_value <- as.numeric(strsplit(number_str, split = '_')[[1]][2])
n_samples <- (max(mcmc_samples$iteration) - burn_in)/thinning_value

## print acceptance rates
for(epid in unique(mcmc_samples$epidemic)){
  filt <- mcmc_samples[epidemic == epid,]
  accept_rate <- n_distinct(filt$likelihood)/nrow(filt)
  filt_burned <- filt[iteration > 0.8*max(filt$iteration),]
  accept_rate_burned <- n_distinct(filt_burned$likelihood)/nrow(filt_burned)
  message('Epidemic: ', epid, '\nOverall acceptance rate: ', round(100*accept_rate, 1), 
          '%\nAcceptance rate in last 20% of steps: ', round(100*accept_rate_burned, 1), '%')
}

mcmc_samples_filtered <- if(save_all_flag){
  # reduce whole chain by thinning value
  mcmc_samples[iteration %% thinning_value == 0 &
                 iteration >= burn_in,] 
}else{
  copy(mcmc_samples)[, iteration := burn_in + iteration*thinning_value]
}

#### ADD Reff #### 
message('Adding Reff')
unique_df <- unique(mcmc_samples_filtered[, ..posterior_cols])

# if(save_all_flag){
  # only do every 100th unique row of mcmc_samples_filtered
  # Reff_filter <- 100
  # unique_df <- unique_df[Reff_filter*(1:(nrow(unique_df)/Reff_filter)),]
# }

pb <- txtProgressBar(min = 1, max = nrow(unique_df), style = 3)

RVACC <- c()

for(i in 1:nrow(unique_df)){
  
  row <- unique_df[i, ]
  
  R1 <- R0_func(susceptibility = fcn_assign_ages(row$children_susceptibility_epid_1,
                                                 row$adults_susceptibility_epid_1,
                                                 row$older_adults_susceptibility_epid_1,
                                                 age_labels),
                inf_period = known_pars$epid_periods[2],
                beta_in = row$transmissibility_epid_1,
                cm_in = cm)
  
  RVACC <- c(RVACC, R0_func(susceptibility = fcn_assign_ages((1-0.248)*row$children_susceptibility_epid_1,
                                           (1-0.14)*row$adults_susceptibility_epid_1,
                                           (1-0.225)*row$older_adults_susceptibility_epid_1,
                                           age_labels),
          inf_period = known_pars$epid_periods[2],
          beta_in = row$transmissibility_epid_1,
          cm_in = cm))
  
  R2 <- R0_func(susceptibility = fcn_assign_ages(row$children_susceptibility_epid_2,
                                                 row$adults_susceptibility_epid_2,
                                                 row$older_adults_susceptibility_epid_2,
                                                 age_labels),
                inf_period = known_pars$epid_periods[2],
                beta_in = row$transmissibility_epid_2,
                cm_in = cm)
  
  mcmc_samples_filtered[transmissibility_epid_1 == row$transmissibility_epid_1 &
                 children_susceptibility_epid_1 == row$children_susceptibility_epid_1 & 
                 adults_susceptibility_epid_1 == row$adults_susceptibility_epid_1 &
                 older_adults_susceptibility_epid_1 == row$older_adults_susceptibility_epid_1,
               Reff_epid_1 := R1]
  
  mcmc_samples_filtered[transmissibility_epid_2 == row$transmissibility_epid_2 &
                 children_susceptibility_epid_2 == row$children_susceptibility_epid_2 &
                 adults_susceptibility_epid_2 == row$adults_susceptibility_epid_2 &
                 older_adults_susceptibility_epid_2 == row$older_adults_susceptibility_epid_2,
               Reff_epid_2 := R2]

  # Print progress
  setTxtProgressBar(pb, i)
  
}
close(pb)

## approximate R after accounting for effective vaccination
hist(RVACC)

posterior_cols <- c(posterior_cols, 'Reff_epid_1', 'Reff_epid_2')
plotting_cols <- unique(gsub('_epid_1|_epid_2', '', posterior_cols))

mcmc_samples_filtered_burned <- if(save_all_flag){
  ## filter past burn-in
  mcmc_samples_filtered[iteration > burn_in,]#[, iteration := 1:n_samples, .(chain, epidemic)]
}else{
  copy(mcmc_samples_filtered)
  # mcmc_samples_filtered[iteration > 3e5,] # manually add 100,000 burn in
}

#### PLOT TRACE #### 

message('Plotting filtered traces')
traces_filtered <- map(.x = plotting_cols, .f = ~{plot_trace(var=.x, filtered=T)})
map(.x = 1:length(traces_filtered), .f = ~{
  ggsave(filename = gsub('data/mcmc_posteriors.rds',figure_filename(plotting_cols[.x]),.args[length(.args)]), 
         plot = traces_filtered[.x],
         width = 14, height = 7)
})
tf <- patchwork::wrap_plots(traces_filtered, nrow = 6) 
ggsave(filename = gsub('data/mcmc_posteriors.rds',figure_filename('fitted_filtered_traces'),.args[length(.args)]), 
       plot = tf,
       width = 45, height = 22)

message('Plotting log-likelihood')
log_likelihood_plot <- plot_trace('likelihood'); suppressWarnings(log_likelihood_plot)
ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('fitted_likelihood'),.args[length(.args)]), 
       width = 12, height = 16)

risk_ratios <- mcmc_samples_filtered_burned %>% 
  select(job, chain, iteration, epidemic, contains('care_rate')) %>% 
  pivot_longer(!c(job, chain, iteration, epidemic)) %>% 
  mutate(setting = case_when(grepl('primary', name) ~ 'primary', grepl('secondary', name) ~ 'secondary'),
         subtype = case_when(grepl('epid_1', name) ~ '1', grepl('epid_2', name) ~ '2'),
         broad_age = case_when(grepl('children', name) ~ 'children', grepl('older_adults', name) ~ 'older_adults', T ~ 'adults'),
         risk = case_when(grepl('high', name) ~ 'high', T ~ 'low')) %>% select(!name) %>% 
  pivot_wider(names_from = risk, values_from = value) %>% 
  mutate(high_over_low = high/low) %>% mutate(subtype_season = paste0(years[epidemic], ' (subtype ', subtype, ')'))

true_ratios <- epid_pars %>% filter(grepl('care_rate', joining_name)) %>%  
  mutate(setting = case_when(grepl('primary', joining_name) ~ 'primary', grepl('secondary', joining_name) ~ 'secondary'),
         subtype = case_when(grepl('epid_1', joining_name) ~ '1', grepl('epid_2', joining_name) ~ '2'),
         broad_age = case_when(grepl('children', joining_name) ~ 'children', grepl('older_adults', joining_name) ~ 'older_adults', T ~ 'adults'),
         risk = case_when(grepl('high', joining_name) ~ 'high', T ~ 'low'), 
         subtype_season = paste0(year, ' (subtype ', subtype, ')')) %>% 
  select(subtype_season, setting, subtype, broad_age, risk, value) %>% 
  pivot_wider(names_from = risk, values_from = value) %>% 
  mutate(true_high_over_low = high/low)

risk_ratios <- risk_ratios %>% 
  left_join(true_ratios %>% select(!c(low,high)), by = c('subtype_season','setting','subtype','broad_age')) %>% 
  filter(!subtype_season == '2025 (subtype 2)')

risk_ratios$broad_age <- factor(risk_ratios$broad_age, levels = c('children','adults','older_adults')) 

risk_ratios %>% 
  ggplot() + 
  geom_line(aes(iteration, high_over_low, color = as.factor(job), group = interaction(job, chain))) + 
  geom_hline(aes(yintercept = true_high_over_low), lty = 2) +
  facet_grid(interaction(broad_age,setting) ~ subtype_season, scales = 'free') + theme_bw() +
  labs(y = 'IGPR or IHR (high risk / low risk)', col = 'job')

#### PLOT DENSITY #### 
# message('Plotting densities')
# densities <- map(.x = plotting_cols, .f = plot_density)
# patchwork::wrap_plots(densities, nrow = 6)
# ggsave(gsub('data/mcmc_posteriors.rds',figure_filename('fitted_densities'),.args[length(.args)]), width = 30, height = 14)

if(F){
## PAIRWISE PLOTS
message('Plotting pairwise correlations')
mcmc_samples_filtered_burned[, chain_it_id := paste0(chain, '_', iteration)]

pairs_cols <- c(posterior_cols, 'epidemic', 'chain_it_id')

pairs_data <- if(nrow(mcmc_samples_filtered_burned) >= 10000){
  # only taking 1% of the mcmc_samples_filtered dataset as it takes too long otherwise!
  mcmc_samples_filtered_burned[seq(1, nrow(mcmc_samples_filtered_burned), by = 100), ..pairs_cols]
}else{
  mcmc_samples_filtered_burned[1:nrow(mcmc_samples_filtered_burned), ..pairs_cols]
}

pairs_long <- melt(pairs_data, id.vars = c('epidemic', 'chain_it_id'))
pairs_long[, variable := as.character(variable)]
pairs_long[, epidemic_of_season := as.numeric(substr(variable, nchar(variable), nchar(variable)))]
pairs_long[substr(variable, 1, 3) == 'imd', epidemic_of_season := 1]
add_imd_data <- pairs_long[substr(variable, 1, 3) == 'imd',]
add_imd_data[, epidemic_of_season := 2]
pairs_long <- rbind(pairs_long,
                    add_imd_data)
pairs_long[substr(variable, 1, 3) != 'imd', variable := substr(variable, 1, nchar(variable) - 7)]

pairs_long <- pairs_long[subtype_seasons[subtype_seasons$epidemic %in% unique(pairs_long$epidemic),], on = c('epidemic','epidemic_of_season')]
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
}

#### SAVE DATA ####

write_rds(mcmc_samples_filtered_burned, .args[length(.args)])






