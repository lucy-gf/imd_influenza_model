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
    file.path("data", "dummy_data", "dummy_infections.rds"),
    file.path("data", "dummy_data", "dummy_surveillance.rds"),
    file.path("data", "dummy_data", "known_parameters.rds"),
    file.path("data", "dummy_data", "unknown_parameters.rds"),
    file.path("output", "data", "mcmc_samples_rates_unknown.rds"),
    file.path("output", "figures", "dummy_mcmc", "fitted_epidemics.png")
  ) #else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

if(!dir.exists(gsub('fitted_epidemics.png','',.args[length(.args)]))){
  dir.create(gsub('fitted_epidemics.png','',.args[length(.args)]), recursive=T)
}

set.seed(60)

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
  broad_age = c(rep('children', 3), rep('adults', 4), rep('older_adults', 2))
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

#### TRUE INFECTIONS DATA #### 
true_infections_list <- readRDS(.args[3])

#### SURVEILLANCE DATA #### 
surveillance_data <- readRDS(.args[4])

#### KNOWN PARAMETERS #### 
known_pars <- readRDS(.args[5])
years <- known_pars$years

risk_group_pop <- known_pars$risk_group_pop 
risk_group_pop$age_grp <- factor(risk_group_pop$age_grp, levels = age_labels)
risk_group_pop <- risk_group_pop %>% 
  arrange(imd_quintile, age_grp)

vaccinated_pop <- known_pars$vaccinated_pop
vaccinated_pop$age_grp <- factor(vaccinated_pop$age_grp, levels = age_labels)
vaccinated_pop <- vaccinated_pop %>% 
  left_join(known_pars$vaccination_efficacy, by = 'age_grp') %>% 
  mutate(effectively_vaccinated_population = VE*vaccinated_population) %>% 
  arrange(desc(risk_level), imd_quintile, age_grp)

demography <- rbind(risk_group_pop %>% mutate(risk_level = 'high'),
                    risk_group_pop %>% mutate(risk_level = 'low')) %>% 
  mutate(population = case_when(risk_level == 'high' ~ risk_population,
                                risk_level == 'low' ~ pop - risk_population)) %>% 
  select(!c(risk_population,pop)) %>% arrange(desc(risk_level), imd_quintile, age_grp)

## check population sum is correct
tot_pop <- sum(imd_age_pop$pop)
if(!all.equal(sum(demography$population), tot_pop)){warning('pop not adding up')}

#### UNKNOWN PARAMETERS #### 
unknown_pars <- readRDS(.args[6])
care_rates <- unknown_pars$care_rates %>% 
  rename(primary_care = gp_rate,
         secondary_care = hosp_rate)

## true epid parameters being fitted
epid_pars <- data.frame()
for(i in 1:3){
  epid_pars <- rbind(epid_pars,
    data.frame(epidemic = i,
               transmissibility = unknown_pars[[paste0('epid_parameters_s',i)]]$transmissibility,
               susceptibility_1 = unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[1],
               susceptibility_2 = unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[3],
               init_infected = log10(unknown_pars[[paste0('epid_parameters_s',i)]]$init_infected),
               primary_care_rate_children_low = unknown_pars$primary_care_rates[1],
               primary_care_rate_adults_low = unknown_pars$primary_care_rates[2],
               primary_care_rate_older_adults_low = unknown_pars$primary_care_rates[3],
               primary_care_rate_children_high = unknown_pars$primary_care_rates[4],
               primary_care_rate_adults_high = unknown_pars$primary_care_rates[5],
               primary_care_rate_older_adults_high = unknown_pars$primary_care_rates[6],
               secondary_care_rate_children_low = unknown_pars$secondary_care_rates[1],
               secondary_care_rate_adults_low = unknown_pars$secondary_care_rates[2],
               secondary_care_rate_older_adults_low = unknown_pars$secondary_care_rates[3],
               secondary_care_rate_children_high = unknown_pars$secondary_care_rates[4],
               secondary_care_rate_adults_high = unknown_pars$secondary_care_rates[5],
               secondary_care_rate_older_adults_high = unknown_pars$secondary_care_rates[6],
               imd_spline_primary_1 = unknown_pars$imd_spline_pars$primary[1],
               imd_spline_primary_2 = unknown_pars$imd_spline_pars$primary[2],
               imd_spline_secondary_1 = unknown_pars$imd_spline_pars$secondary[1],
               imd_spline_secondary_2 = unknown_pars$imd_spline_pars$secondary[2],
               R0 = R0_func(susceptibility = c(rep(unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[1],3),
                                               rep(unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[2],4),
                                               rep(unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[3],2)),
                            inf_period = known_pars$epid_periods[2],
                            beta_in = unknown_pars[[paste0('epid_parameters_s',i)]]$transmissibility,
                            cm_in = cm)
    ))
}
epid_pars <- epid_pars %>% pivot_longer(!epidemic)

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
  dat <- readRDS(gsub('.rds',paste0('_', i, '_', number_str,'.rds'),.args[7]))
  list_samples <- mclapply(1:3, get_samples_parallel)
  samples_out <- rbindlist(list_samples)
  samples_out[, epidemic := i]
  samples_out
}

# load most recently run settings (burn in, thinning, samples)
output_details_file <- readRDS(.args[7])
number_str <- output_details_file$x[1]
run_date <- output_details_file$date
cat('\n------------\n',run_date,'\n------------\n',sep='')
cat('\n------------\n',number_str,'\n------------\n',sep='')
number_date_str <- paste0(number_str, '_', date)

# was it run on the HPC? The files saved differently
WAS_HPC <- output_details_file$HPC

if(WAS_HPC){
  
  dat_example <- readRDS(gsub('.rds',paste0('_1_', number_date_str,'.rds'),.args[7]))
  mcmc_samples <- rbindlist(lapply(1:3, read_and_get_samples))
  mcmc_samples[, c('LP','LPr') := NULL]
  
}else{
  
  mcmc_samples_list <- readRDS(gsub('.rds',paste0('_', number_date_str,'.rds'),.args[7]))
  mcmc_samples <- rbindlist(lapply(1:3, get_samples))
  mcmc_samples[, c('Lposterior','Lprior') := NULL]
  
}

fitted_pars <- unique(epid_pars$name)
fitted_pars <- fitted_pars[fitted_pars %notin% c('epidemic', 'R0')]
cat('\n',length(fitted_pars),' fitted parameters\n', sep = '')
colnames(mcmc_samples) <- c(fitted_pars, 'likelihood', 'chain', 'iteration', 'epidemic')

#### ADD R0 #### 
cat('Adding R0: ')
unique_df <- unique(mcmc_samples[, ..fitted_pars])
mod_val <- if(nrow(unique_df) > 50000){1000}else{
  ifelse(nrow(unique_df) > 5000,100,10)}
for(i in 1:nrow(unique_df)){
  row <- unique_df[i, ]
  mcmc_samples[transmissibility == row$transmissibility &
                 susceptibility_1 == row$susceptibility_1 &
                 susceptibility_2 == row$susceptibility_2,
               R0 := R0_func(susceptibility = susc_vector(c(row$susceptibility_1,
                                                            row$susceptibility_2)),
                             inf_period = known_pars$epid_periods[2],
                             beta_in = row$transmissibility,
                             cm_in = cm)
               ]
  if(i %% 100 == mod_val){cat(round(100*i/nrow(unique_df)), '%, ', sep='')}
}

fitted_pars <- c(fitted_pars, 'R0')

#### FILTER #### 
burn_in <- as.numeric(strsplit(number_str, split = '_')[[1]][1])
thinning_value <- as.numeric(strsplit(number_str, split = '_')[[1]][2])
n_samples <- (max(mcmc_samples$iteration) - burn_in)/thinning_value

mcmc_samples_filtered <- mcmc_samples[iteration > burn_in & iteration %% thinning_value == 0,]
mcmc_samples_filtered[, iteration := 1:n_samples, .(chain, epidemic)]

#### PLOT DENSITY #### 
densities <- map(.x = fitted_pars, .f = plot_density)
patchwork::wrap_plots(densities, nrow = 3)
ggsave(gsub('epidemics','densities',.args[length(.args)]), width = 30, height = 14)

#### PLOT TRACE #### 
traces <- map(.x = fitted_pars, .f = plot_trace)
patchwork::wrap_plots(traces, nrow = 3) 
ggsave(gsub('epidemics','traces',.args[length(.args)]), width = 30, height = 14)
traces_filtered <- map(.x = fitted_pars, .f = ~{plot_trace(var=.x, filtered=T)})
patchwork::wrap_plots(traces_filtered, nrow = 3) 
ggsave(gsub('epidemics','filtered_traces',.args[length(.args)]), width = 30, height = 14)
log_likelihood_plot <- plot_trace('likelihood'); log_likelihood_plot
ggsave(gsub('epidemics','likelihood',.args[length(.args)]), width = 8, height = 6)

## PAIRWISE PLOTS
pairs_data <- if(nrow(mcmc_samples_filtered) >= 10000){
  # only taking 1% of the mcmc_samples_filtered dataset as it takes too long otherwise!
  mcmc_samples_filtered[seq(1, nrow(mcmc_samples_filtered), by = 100), c(1:20, 24)]
}else{
  mcmc_samples_filtered[1:nrow(mcmc_samples_filtered), c(1:20, 24)]
}

colnames(pairs_data) <- gsub('_rate_','_rate\n_', colnames(pairs_data))
colnames(pairs_data) <- gsub('_spline_','_spline\n_', colnames(pairs_data))
pairs_data$iteration <- rep(1:(nrow(pairs_data)/n_distinct(pairs_data$epidemic)), n_distinct(pairs_data$epidemic))
p <- ggpairs(pairs_data, columns = 1:4, aes(color = as.factor(epidemic), alpha = 0.5))
print(p)
p_FULL <- ggpairs(pairs_data, columns = 1:20, aes(color = as.factor(epidemic), alpha = 0.5))
print(p_FULL)
ggsave(gsub('epidemics','pairwise',.args[length(.args)]), width = 30, height = 30)

#### RUN FITTED EPIDEMICS #### 

pop_vaccinated <- vaccinated_pop$effectively_vaccinated_population

## only do for every tenth/hundredth/thousandth epidemic to save time
mod_val <- if(nrow(mcmc_samples_filtered) > 100000){1000}else{
  ifelse(nrow(mcmc_samples_filtered) > 10000,100,10)}

cat('\nRunning fitted epidemics (', nrow(mcmc_samples_filtered)/mod_val, ' total): ', sep = '')
fitted_epidemics <- data.table()
for(bs in 1:nrow(mcmc_samples_filtered)){ 
  
  if(bs %% mod_val != 0){next} 
  
  init_infected_num <- 10^(mcmc_samples_filtered$init_infected[bs])
  init_infected_vec <- (demography$population - pop_vaccinated)*init_infected_num/
    (sum(demography$population)-sum(pop_vaccinated))
  
  time_series <- run_model(
    pop = demography$population,
    I0 = init_infected_vec,
    vacc = pop_vaccinated,
    cm = pc_cm,
    trans = mcmc_samples_filtered$transmissibility[bs],
    susc = susc_vector(mcmc_samples_filtered[bs, paste0('susceptibility_',1:2)]),
    lat_per = known_pars$epid_periods[1],
    inf_per = known_pars$epid_periods[2]
  )
  
  fitted_epidemics <- rbind(
    fitted_epidemics,
    time_series[, chain := mcmc_samples_filtered$chain[bs]][, iteration := mcmc_samples_filtered$iteration[bs]][, epidemic := mcmc_samples_filtered$epidemic[bs]]
  )
  
  if(bs %% (10*mod_val) == 0){cat(bs/mod_val, ', ', sep = '')}
  
}

## make long reporting rates data table

## column names
primary_names <- paste0('primary_imd_', 1:5)
secondary_names <- paste0('secondary_imd_', 1:5)

## apply imd_spline function
rep_rates <- mcmc_samples_filtered %>% 
  bind_cols(pmap(list(mcmc_samples_filtered$imd_spline_primary_1, mcmc_samples_filtered$imd_spline_primary_2),
      function(x, y) {
        result <- imd_spline(c(x, y))
        setNames(as.list(result), primary_names)
      }) %>% bind_rows()) %>% 
  bind_cols(pmap(list(mcmc_samples_filtered$imd_spline_secondary_1, mcmc_samples_filtered$imd_spline_secondary_2),
                 function(x, y) {
                   result <- imd_spline(c(x, y))
                   setNames(as.list(result), secondary_names)
                 }) %>% bind_rows())

age_risk_rates <- rep_rates %>% 
  select(contains('rate'), chain, iteration, epidemic) %>%
  pivot_longer(!c(chain,iteration,epidemic)) %>% 
  mutate(care_setting = case_when(substr(name, 1, 1) == 'p' ~ 'primary_rate', T ~ 'secondary_rate'),
         risk_level = case_when(grepl('low', name) ~ 'low', T ~ 'high'),
         broad_age = case_when(grepl('children', name) ~ 'children', 
                               grepl('older_adults', name) ~ 'older_adults', 
                               T ~ 'adults')) %>% 
  select(!name)

imd_rates <- rep_rates %>% 
  select(contains('ary_imd'), chain, iteration, epidemic) %>%
  pivot_longer(!c(chain,iteration,epidemic)) %>% 
  mutate(care_setting = case_when(substr(name, 1, 1) == 'p' ~ 'primary_rate', T ~ 'secondary_rate'),
         imd_quintile = as.factor(gsub('primary_imd_|secondary_imd_','',name))) %>% 
  select(!name)

mcmc_surveillance_rates <- full_join(
  age_risk_rates, imd_rates, by = c('chain','iteration','epidemic','care_setting'), relationship = "many-to-many"
) %>% mutate(value = value.x*value.y) %>% select(!c(value.x, value.y))

mcmc_surveillance_rates_w <- mcmc_surveillance_rates %>% 
  pivot_wider(names_from = 'care_setting', values_from = value)

#### PLOT HEALTHCARE RATE ASCERTAINMENT #### 

mcmc_surveillance_rates_w$broad_age <- factor(
  mcmc_surveillance_rates_w$broad_age, 
  levels = c('children','adults','older_adults')
)

hc_rates_fitted <- mcmc_surveillance_rates_w %>% 
  pivot_longer(c(primary_rate, secondary_rate)) %>% 
  group_by(risk_level, broad_age, imd_quintile, name) %>% 
  summarise(median = median(value),
            l = quantile(value, 0.025),
            u = quantile(value, 0.975)) %>% 
  left_join(care_rates %>% 
              select(!contains('rel_')) %>% 
              mutate(broad_age = case_when(
                age_grp %in% c('0-4','5-11','12-17') ~ 'children',
                age_grp %in% c('70-79','80+') ~ 'older_adults', 
                T ~ 'adults'
              )) %>% select(!age_grp) %>% unique() %>% 
              rename(primary_rate = primary_care,
                     secondary_rate = secondary_care) %>% 
              pivot_longer(c(primary_rate, secondary_rate)) %>% 
              rename(true_value = value) %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('risk_level', 'imd_quintile', 'broad_age', 'name'))

hc_rates_fitted$broad_age <- factor(
  hc_rates_fitted$broad_age, 
  levels = c('children','adults','older_adults')
)

hc_rates_fitted %>% 
  mutate(name = gsub('_rate', ' care', name),
         risk_level = paste0(risk_level, ' risk')) %>% 
  ggplot(aes(x = broad_age, group = imd_quintile, col = imd_quintile)) + 
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.4,
                position = position_dodge(width = 0.4), alpha=1) +
  geom_point(aes(y = median), #shape = 1, 
             position = position_dodge(width = 0.4)) +
  geom_point(aes(y = true_value), shape = 4, size = 3, 
             position = position_dodge(width = 0.4), stroke = 0.8) +
  theme_bw() +
  scale_color_manual(values = imd_quintile_colors) + 
  facet_wrap(risk_level ~ name, scales = 'free', ncol = 2) + 
  theme(legend.position = 'none') +
  scale_y_continuous(labels = scales::percent, limits = c(0,NA)) +
  labs(x = '', y = 'Healthcare attendance upon infection')

# hc_rates_fitted_2 <- mcmc_surveillance_rates_w %>% 
#   pivot_longer(c(primary_rate, secondary_rate)) %>% 
#   left_join(care_rates %>% 
#               select(!contains('rel_')) %>% 
#               mutate(broad_age = case_when(
#                 age_grp %in% c('0-4','5-11','12-17') ~ 'children',
#                 age_grp %in% c('70-79','80+') ~ 'older_adults', 
#                 T ~ 'adults'
#               )) %>% select(!age_grp) %>% unique() %>% 
#               rename(primary_rate = primary_care,
#                      secondary_rate = secondary_care) %>% 
#               pivot_longer(c(primary_rate, secondary_rate)) %>% 
#               rename(true_value = value) %>% mutate(imd_quintile = as.factor(imd_quintile)), 
#             by = c('risk_level', 'imd_quintile', 'broad_age', 'name')) %>% 
#   mutate(name = gsub('_rate', ' care', name),
#          risk_level = paste0(risk_level, ' risk')) 
# 
# hc_rates_fitted_2$broad_age <- factor(
#   hc_rates_fitted_2$broad_age, 
#   levels = c('children','adults','older_adults')
# )
# 
# hc_rates_fitted_2 %>% 
#   ggplot() + 
#   geom_violin(aes(x = broad_age, y = value, group = interaction(imd_quintile, broad_age),
#                   col = imd_quintile, fill = imd_quintile), 
#               alpha = 0.4,
#               position = position_dodge(width = 0.4)) +
#   geom_point(aes(y = true_value, x = broad_age, group = imd_quintile, col = imd_quintile), 
#              shape = 4, size = 3, 
#              position = position_dodge(width = 0.4), stroke = 0.8) +
#   theme_bw() +
#   scale_color_manual(values = imd_quintile_colors) +
#   scale_fill_manual(values = imd_quintile_colors) +
#   facet_wrap(risk_level ~ name, scales = 'free', ncol = 2) + 
#   theme(legend.position = 'none') +
#   scale_y_continuous(labels = scales::percent, limits = c(0,NA)) +
#   labs(x = '', y = 'Healthcare attendance upon infection')

## add in start_date, reporting rates
fitted_epidemics[, start_of_epidemic := as.Date(paste0('01-09-', years[1] + epidemic - 1), format = '%d-%m-%Y')]
fitted_epidemics[, date := start_of_epidemic + t] ## add date
fitted_epidemics[, c('t','start_of_epidemic') := NULL]

# aggregate week
fitted_epidemics[, date := last_monday(date)]
fitted_epidemics_agg <- fitted_epidemics[, .(infections = ceiling(sum(infections))),
                                             by = .(date, epidemic, chain, iteration, age_grp, imd_quintile, risk_level)]

fitted_epidemics_surv <- copy(fitted_epidemics_agg) ## copy for surveillance data

fitted_epidemics_agg[, c('epidemic','iteration','chain') := NULL]
fitted_epidemics_agg_m <- rbind(
  fitted_epidemics_agg[, lapply(.SD, median), by = c('date', 'age_grp', 'imd_quintile', 'risk_level')][, measure := 'median'],
  fitted_epidemics_agg[, lapply(.SD, max), by = c('date', 'age_grp', 'imd_quintile', 'risk_level')][, measure := 'u'],
  fitted_epidemics_agg[, lapply(.SD, min), by = c('date', 'age_grp', 'imd_quintile', 'risk_level')][, measure := 'l'])

fitted_epidemics_agg_l <- dcast(fitted_epidemics_agg_m, 
                                date + age_grp + imd_quintile + risk_level ~ measure,
                                value.var = 'infections')
# fitted_epidemics_agg_l <- fitted_epidemics_agg_m %>% 
#   pivot_wider(names_from = measure, values_from = infections)

## true infections
true_infections <- rbindlist(true_infections_list)[, date := start_date + t][, c('date','age_grp','imd_quintile','infections','risk_level')]
true_infections <- true_infections[, lapply(.SD, sum), by = c('date','age_grp','imd_quintile','risk_level')]
 
## make weekly
true_infections <- true_infections[, date := last_monday(date)][, lapply(.SD, sum), by = c('date','age_grp','imd_quintile','risk_level')]

## combine
fitted_and_obs <- fitted_epidemics_agg_l %>% 
  left_join(true_infections, by = c('date','age_grp','imd_quintile','risk_level')) 

#### PLOT FITTED EPIDEMICS #### 
fitted_and_obs %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=interaction(imd_quintile,risk_level), fill=imd_quintile), alpha=0.4) +
  geom_line(aes(date, median, group=interaction(imd_quintile,risk_level), col=imd_quintile)) +
  geom_point(aes(date, infections, group=interaction(imd_quintile,risk_level), col=imd_quintile)) +
  facet_wrap(age_grp~., scales='free') + theme_bw() +
  scale_color_manual(values = imd_quintile_colors) +
  scale_fill_manual(values = imd_quintile_colors)

fitted_and_obs %>% filter(imd_quintile == 3, age_grp=='5-11', risk_level == 'low') %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=age_grp, fill=age_grp), alpha=0.4) +
  geom_line(aes(date, median, group=age_grp, col=age_grp)) +
  geom_point(aes(date, infections, group=age_grp, col=age_grp), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  theme(legend.position = 'none',
        text = element_text(size=14)) +
  labs(y = 'infections')

fitted_and_obs %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=risk_level, fill=age_grp), alpha=0.4) +
  geom_line(aes(date, median, group=risk_level, col=age_grp, lty = risk_level)) +
  geom_point(aes(date, infections, group=risk_level, col=age_grp), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  scale_color_manual(values = age_colors) +
  scale_linetype_manual(values = c(2,1)) +
  scale_fill_manual(values = age_colors) +
  theme(legend.position = 'none',
        text = element_text(size=14)) +
  labs(y = 'infections')

fitted_and_obs %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  filter(imd_quintile == 'IMD 1', age_grp == '5-11') %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=interaction(age_grp,risk_level), fill=risk_level), 
              alpha=0.4) +
  geom_line(aes(date, median, group=interaction(age_grp,risk_level), col=risk_level), lwd = 0.8) +
  geom_point(aes(date, infections, group=interaction(age_grp,risk_level), col=risk_level), 
             shape = 1, stroke = 1) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  theme(text = element_text(size=12)) +
  scale_x_date(breaks = "1 year", labels=date_format("%Y")) +
  labs(y = 'Infections (originally unobserved)', x = '', col = 'Clinical risk', fill = 'Clinical risk')

## SAVE PNG
ggsave(.args[length(.args)], width = 16, height = 12)

## CRUDE ATTACK RATE ESTIMATES
## (doing sum of median instead of median of sum, for ease of calculation (preliminary))
fitted_and_obs %>% group_by(age_grp, imd_quintile) %>% 
  summarise(med = sum(median)/3, 
            l = sum(l)/3,
            u = sum(u)/3,
            inf = sum(infections/3)) %>% 
  left_join(imd_age_pop %>% mutate(imd_quintile=as.character(imd_quintile)), 
            by = c('age_grp', 'imd_quintile')) %>% 
  ggplot() + 
  geom_bar(aes(x=age_grp, y=med/pop, fill=imd_quintile),
           stat='identity',position='dodge') +
  geom_errorbar(aes(x=age_grp, ymin=l/pop, ymax=u/pop, group=imd_quintile),
                position = position_dodge(width = 0.9), width=0.4, alpha=0.7) +
  theme_bw() + scale_fill_manual(values = imd_quintile_colors) +
  labs(x='', fill='IMD quintile', y='Attack rate')

## fitted surveillance data

fitted_surv_dat <- fitted_epidemics_surv %>% 
  left_join(broad_ages, by = 'age_grp') %>% 
  left_join(surveillance_data %>% 
              rename(date = week_start, epidemic = index) %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('date','age_grp','imd_quintile','risk_level', 'epidemic')) %>% 
  left_join(mcmc_surveillance_rates_w,
            by = c('broad_age','imd_quintile','risk_level','chain','iteration','epidemic')) %>% 
  left_join(known_pars$proportion_observed %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('age_grp','imd_quintile','risk_level')) %>% 
  mutate(observed_primary = infections*primary_rate*OS_COVERAGE,
         observed_secondary = infections*secondary_rate*OS_COVERAGE) %>% 
  group_by(iteration, chain, date, age_grp, imd_quintile, 
           risk_level, primary_care, secondary_care) %>% 
  summarise(observed_primary = sum(observed_primary),
            observed_secondary = sum(observed_secondary))
fitted_surv_dat <- data.table(fitted_surv_dat)

fitted_surv_dat[, c('iteration', 'chain') := NULL]

fitted_surv_agg <- rbind(
  fitted_surv_dat[, lapply(.SD, median), by = c('date', 'age_grp', 'imd_quintile', 'risk_level', 'primary_care', 'secondary_care')][, measure := 'median'],
  fitted_surv_dat[, lapply(.SD, max), by = c('date', 'age_grp', 'imd_quintile', 'risk_level', 'primary_care', 'secondary_care')][, measure := 'u'],
  fitted_surv_dat[, lapply(.SD, min), by = c('date', 'age_grp', 'imd_quintile', 'risk_level', 'primary_care', 'secondary_care')][, measure := 'l'])

fitted_surv_agg <- fitted_surv_agg %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  group_by(age_grp, imd_quintile, risk_level) %>% 
  mutate(observed_primary = lag(observed_primary, default=0),
         observed_secondary = lag(lag(observed_secondary, default=0), default=0)) 

fitted_primary_plot <- fitted_surv_agg %>%
  select(!observed_secondary) %>% 
  pivot_wider(names_from = measure, values_from = observed_primary) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=interaction(age_grp,risk_level), fill=risk_level), alpha=0.4) +
  geom_line(aes(date, median, group=interaction(age_grp,risk_level), col=risk_level)) +
  geom_point(aes(date, primary_care, group=interaction(age_grp,risk_level), col=risk_level), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  theme(text = element_text(size=14)) +
  scale_x_date(breaks = "1 year", labels=date_format("%Y")) +
  labs(y = 'primary care', x = ''); fitted_primary_plot

fitted_secondary_plot <- fitted_surv_agg %>%
  select(!observed_primary) %>% 
  pivot_wider(names_from = measure, values_from = observed_secondary) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=interaction(age_grp,risk_level), fill=risk_level), alpha=0.4) +
  geom_line(aes(date, median, group=interaction(age_grp,risk_level), col=risk_level)) +
  geom_point(aes(date, secondary_care, group=interaction(age_grp,risk_level), col=risk_level), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  theme(text = element_text(size=14)) +
  scale_x_date(breaks = "1 year", labels=date_format("%Y")) +
  labs(y = 'primary care', x = ''); fitted_secondary_plot

fitted_surv_agg %>%
  select(!observed_primary) %>% 
  filter(imd_quintile == 'IMD 1', age_grp == '5-11') %>% 
  pivot_wider(names_from = measure, values_from = observed_secondary) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=interaction(age_grp,risk_level), fill=risk_level), 
              alpha=0.4) +
  geom_line(aes(date, median, group=interaction(age_grp,risk_level), col=risk_level), lwd = 0.8) +
  geom_point(aes(date, secondary_care, group=interaction(age_grp,risk_level), col=risk_level), 
             shape = 1, stroke = 1) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  theme(text = element_text(size=12)) +
  scale_x_date(breaks = "1 year", labels=date_format("%Y")) +
  labs(y = 'Hospitalisations', x = '', col = 'Clinical risk', fill = 'Clinical risk')

fitted_surv_agg %>%
  select(!observed_secondary) %>% 
  filter(imd_quintile == 'IMD 1', age_grp == '5-11') %>% 
  pivot_wider(names_from = measure, values_from = observed_primary) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=interaction(age_grp,risk_level), fill=risk_level), 
              alpha=0.4) +
  geom_line(aes(date, median, group=interaction(age_grp,risk_level), col=risk_level), lwd = 0.8) +
  geom_point(aes(date, primary_care, group=interaction(age_grp,risk_level), col=risk_level), 
             shape = 1, stroke = 1) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  theme(text = element_text(size=12)) +
  scale_x_date(breaks = "1 year", labels=date_format("%Y")) +
  labs(y = 'Primary care', x = '', col = 'Clinical risk', fill = 'Clinical risk')

## SAVE PNGs
fitted_primary_plot
ggsave(gsub('epidemics','primary_data',.args[length(.args)]), width = 16, height = 12)
fitted_secondary_plot
ggsave(gsub('epidemics','secondary_data',.args[length(.args)]), width = 16, height = 12)

