## MCMC PLOTTING ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(readr))
suppressMessages(require(BayesianTools))
suppressMessages(require(patchwork))
suppressMessages(require(parallel))
options(dplyr.summarise.inform = FALSE) 

.args <- #if (interactive()) c(
  c(file.path("data", "inputs", "imd_age_pop.rds"),
    file.path("data", "inputs", "contact_matrix.rds"),
    file.path("data", "dummy_data", "dummy_infections.rds"),
    file.path("data", "dummy_data", "dummy_surveillance.rds"),
    file.path("data", "dummy_data", "known_parameters.rds"),
    file.path("data", "dummy_data", "unknown_parameters.rds"),
    file.path("output", "data", "mcmc_samples.rds"),
    file.path("output", "figures", "dummy_mcmc", "fitted_epidemics.png")
  ) #else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

if(!dir.exists(gsub('fitted_epidemics.png','',.args[length(.args)]))){
  dir.create(gsub('fitted_epidemics.png','',.args[length(.args)]))
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
               susceptibility_1 = unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[1],
               susceptibility_2 = unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[2],
               susceptibility_3 = unique(unknown_pars[[paste0('epid_parameters_s',i)]]$susceptibility)[3],
               transmissibility = unknown_pars[[paste0('epid_parameters_s',i)]]$transmissibility,
               start_date = as.numeric(unknown_pars[[paste0('epid_parameters_s',i)]]$start_date - 
                                         as.Date(paste0('01-01-',year(unknown_pars[[paste0('epid_parameters_s',i)]]$start_date)),
                                                 format = '%d-%m-%Y')),
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

mcmc_samples_list <- readRDS(gsub('.rds','_1000_10_1000.rds',.args[7]))
mcmc_samples <- rbindlist(lapply(1:3, get_samples))

mcmc_samples[, c('Lprior','Lposterior') := NULL]

fitted_pars <- c('transmissibility', paste0('susceptibility_', 1:3), 'start_date')
colnames(mcmc_samples) <- c(fitted_pars, 'llikelihood', 'chain', 'iteration', 'epidemic')

#### ADD R0 #### 
cat('Adding R0: ')
unique_df <- unique(mcmc_samples[, ..fitted_pars])
for(i in 1:nrow(unique_df)){
  row <- unique_df[i, ]
  mcmc_samples[transmissibility == row$transmissibility &
                 susceptibility_1 == row$susceptibility_1 &
                 susceptibility_2 == row$susceptibility_2 &
                 susceptibility_3 == row$susceptibility_3, 
               R0 := R0_func(susceptibility = c(rep(row$susceptibility_1,3),
                                                rep(row$susceptibility_2,4),
                                                rep(row$susceptibility_3,2)),
                             inf_period = known_pars$epid_periods[2],
                             beta_in = row$transmissibility,
                             cm_in = cm)
               ]
  if(i %% 100 == 0){cat(i, ', ', sep='')}
}

fitted_pars <- c(fitted_pars, 'R0')

#### FILTER #### 
burn_in <- mcmc_samples_list$burn_in
thinning_value <- mcmc_samples_list$thinning_value
n_samples <- mcmc_samples_list$n_samples

mcmc_samples_filtered <- mcmc_samples[iteration > burn_in & iteration %% thinning_value == 0,]
mcmc_samples_filtered[, iteration := 1:n_samples, .(chain, epidemic)]

var_cols <- c('#9B7EDE','#832161','#F7B801','#52050A','#BCD2EE','#F35B04')
names(var_cols) <- fitted_pars

densities <- map(.x = fitted_pars, .f = plot_density)
patchwork::wrap_plots(densities, nrow = 2)
ggsave(gsub('epidemics','densities',.args[length(.args)]), width = 12, height = 8)

traces <- map(.x = fitted_pars, .f = plot_trace)
patchwork::wrap_plots(traces, nrow = 2) 
ggsave(gsub('epidemics','traces',.args[length(.args)]), width = 15, height = 8)
log_likelihood_plot <- plot_trace('llikelihood'); log_likelihood_plot
ggsave(gsub('epidemics','llikelihood',.args[length(.args)]), width = 8, height = 6)

mcmc_samples_filtered %>%
  ggplot() +
  geom_point(aes(x = transmissibility, y = susceptibility_1, col = as.factor(epidemic))) +
  geom_point(aes(x = transmissibility, y = susceptibility_2, col = as.factor(epidemic)), shape = 1) +
  geom_point(aes(x = transmissibility, y = susceptibility_3, col = as.factor(epidemic)), shape = 2) +
  theme_bw() + labs(col='Epidemic')

#### RUN FITTED EPIDEMICS #### 

pop_vaccinated <- vaccinated_pop$effectively_vaccinated_population
init_infected_num <- known_pars$init_infected
init_infected_vec <- (demography$population - pop_vaccinated)*init_infected_num/
  (sum(demography$population)-sum(pop_vaccinated))

fitted_epidemics <- data.table()

for(bs in 1:nrow(mcmc_samples_filtered)){ 
  
  if(bs %% 10 != 0){next} ## only do for every tenth epidemic to save time
  
  time_series <- run_model(
    pop = demography$population,
    I0 = init_infected_vec,
    vacc = pop_vaccinated,
    cm = pc_cm,
    trans = mcmc_samples_filtered$transmissibility[bs],
    susc = susc_vector(mcmc_samples_filtered[bs, paste0('susceptibility_',1:3)]),
    lat_per = known_pars$epid_periods[1],
    inf_per = known_pars$epid_periods[2]
  )
  
  fitted_epidemics <- rbind(
    fitted_epidemics,
    time_series[, chain := mcmc_samples_filtered$chain[bs]][, iteration := mcmc_samples_filtered$iteration[bs]][, epidemic := mcmc_samples_filtered$epidemic[bs]]
  )
  
  if(bs %% 500 == 0){cat(bs, ', ', sep = '')}
  
}

## save
write_rds(fitted_epidemics, gsub('mcmc','epidemic',.args[7]))

## add in start_date
fitted_epidemics_agg <- fitted_epidemics[mcmc_samples_filtered[iteration %in% unique(fitted_epidemics$iteration), 
                                                               c('iteration','chain','start_date','epidemic')], on = c('iteration','chain','epidemic')]
fitted_epidemics_agg[, start_of_year := as.Date(paste0('01-01-', years[1] + epidemic - 1), format = '%d-%m-%Y')]
fitted_epidemics_agg[, date := start_of_year + (floor(start_date) - 1) + t] ## add date
fitted_epidemics_agg[, c('t','start_of_year','start_date','risk_level','epidemic') := NULL]

# aggregate over risk_level and week
fitted_epidemics_agg[, date := last_monday(date)]
fitted_epidemics_agg <- fitted_epidemics_agg[, .(infections = ceiling(sum(infections))),
                                             by = .(date, age_grp, chain, iteration, imd_quintile)]

fitted_epidemics_agg[, c('iteration','chain') := NULL]
fitted_epidemics_agg_m <- rbind(
  fitted_epidemics_agg[, lapply(.SD, median), by = c('date', 'age_grp', 'imd_quintile')][, measure := 'median'],
  fitted_epidemics_agg[, lapply(.SD, max), by = c('date', 'age_grp', 'imd_quintile')][, measure := 'u'],
  fitted_epidemics_agg[, lapply(.SD, min), by = c('date', 'age_grp', 'imd_quintile')][, measure := 'l'])

fitted_epidemics_agg_l <- dcast(fitted_epidemics_agg_m, 
                                date + age_grp + imd_quintile ~ measure,
                                value.var = 'infections')
# fitted_epidemics_agg_l <- fitted_epidemics_agg_m %>% 
#   pivot_wider(names_from = measure, values_from = infections)

## true infections
true_infections <- rbindlist(true_infections_list)[, date := start_date + t][, c('date','age_grp','imd_quintile','infections')]
true_infections <- true_infections[, lapply(.SD, sum), by = c('date','age_grp','imd_quintile')]
 
## make weekly
true_infections <- true_infections[, date := last_monday(date)][, lapply(.SD, sum), by = c('date','age_grp','imd_quintile')]

#### PLOT FITTED EPIDEMICS #### 
fitted_epidemics_agg_l %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=imd_quintile, fill=imd_quintile), alpha=0.4) +
  geom_line(aes(date, median, group=imd_quintile, col=imd_quintile)) +
  geom_point(data = true_infections, aes(date, infections, group=imd_quintile, col=imd_quintile)) +
  facet_wrap(age_grp~., scales='free') + theme_bw() +
  scale_color_manual(values = imd_quintile_colors) +
  scale_fill_manual(values = imd_quintile_colors)

fitted_epidemics_agg_l %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=l, ymax=u, group=age_grp, fill=age_grp), alpha=0.4) +
  geom_line(aes(date, median, group=age_grp, col=age_grp)) +
  geom_point(data = true_infections %>% mutate(imd_quintile := paste0('IMD ', imd_quintile)), 
             aes(date, infections, group=age_grp, col=age_grp), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  theme(legend.position = 'none',
        text = element_text(size=14)) +
  labs(y = 'infections')

## SAVE PNG
ggsave(.args[length(.args)], width = 16, height = 12)

## fitted surveillance data

fitted_primary_plot <- fitted_epidemics_agg_l %>% 
  left_join(surveillance_data %>% 
              rename(date = week_start) %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('date','age_grp','imd_quintile')) %>% 
  left_join(care_rates %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('age_grp','imd_quintile'), suffix = c('','_rate')) %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  group_by(age_grp, imd_quintile) %>% 
  mutate(fitted_pc = lag(median, default=0)*primary_care_rate,
         fitted_pc_l = lag(l, default=0)*primary_care_rate,
         fitted_pc_u = lag(u, default=0)*primary_care_rate) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=fitted_pc_l, ymax=fitted_pc_u, group=age_grp, fill=age_grp), alpha=0.4) +
  geom_line(aes(date, fitted_pc, group=age_grp, col=age_grp)) +
  geom_point(aes(date, primary_care, group=age_grp, col=age_grp), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  theme(legend.position = 'none',
        text = element_text(size=14)) +
  labs(y = 'primary care'); fitted_primary_plot

fitted_secondary_plot <- fitted_epidemics_agg_l %>% 
  left_join(surveillance_data %>% 
              rename(date = week_start) %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('date','age_grp','imd_quintile')) %>% 
  left_join(care_rates %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('age_grp','imd_quintile'), suffix = c('','_rate')) %>% 
  mutate(imd_quintile := paste0('IMD ', imd_quintile)) %>% 
  group_by(age_grp, imd_quintile) %>% 
  mutate(fitted_sc = lag(lag(median, default=0), default=0)*secondary_care_rate,
         fitted_sc_l = lag(lag(l, default=0), default=0)*secondary_care_rate,
         fitted_sc_u = lag(lag(u, default=0), default=0)*secondary_care_rate) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin=fitted_sc_l, ymax=fitted_sc_u, group=age_grp, fill=age_grp), alpha=0.4) +
  geom_line(aes(date, fitted_sc, group=age_grp, col=age_grp)) +
  geom_point(aes(date, secondary_care, group=age_grp, col=age_grp), shape = 1, alpha = 0.6) +
  facet_grid(age_grp~imd_quintile, scales='free') + theme_bw() +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  theme(legend.position = 'none',
        text = element_text(size=14)) +
  labs(y = 'secondary care'); fitted_secondary_plot

## SAVE PNGs
fitted_primary_plot
ggsave(gsub('epidemics','primary_data',.args[length(.args)]), width = 16, height = 12)
fitted_secondary_plot
ggsave(gsub('epidemics','secondary_data',.args[length(.args)]), width = 16, height = 12)

