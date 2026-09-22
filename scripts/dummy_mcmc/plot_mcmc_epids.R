#### RUN FITTED EPIDEMICS #### 

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
    file.path("output", "data", "mcmc_posteriors.rds"),
    file.path("output", "data", "posterior_epidemics.rds")
  ) #else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))
source(file.path('scripts','seir_model.R'))
source(file.path('scripts','dummy_mcmc','mcmc_functions.R'))

set.seed(60)

figure_filename <- function(string){
  paste0('figures/dummy_mcmc/epidemics/',string,'.png')
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

vaccinated_data <- known_pars$vaccinated_data
vaccinated_data$age_grp <- factor(vaccinated_data$age_grp, levels = age_labels)
vaccinated_data <- vaccinated_data %>% 
  arrange(desc(risk_level), imd_quintile, age_grp)

demography <- vaccinated_data %>% 
  mutate(population = pop) %>% 
  select(age_grp, imd_quintile, risk_level, population, risk_proportion) %>% 
  arrange(desc(risk_level), imd_quintile, age_grp) %>% 
  unique()

## check population sum is correct
tot_pop <- sum(imd_age_pop$pop)
if(!all.equal(sum(demography$population), tot_pop)){warning('pop not adding up')}

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
                                adult_susceptibility = up_ss$susceptibility,
                                rel_children_susceptibility = up_ss$rel_susc_children,
                                rel_older_adults_susceptibility = up_ss$rel_susc_older_adults,
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
                                R0 = R0_func(susceptibility = fcn_assign_ages(up_ss$susceptibility*up_ss$rel_susc_children,
                                                                              up_ss$susceptibility,
                                                                              up_ss$susceptibility*up_ss$rel_susc_older_adults,
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
      T ~ paste0(name, '_epid_', epidemic_of_season)),
    subtype = case_when(
      grepl('imd_spline', name) ~ NA, 
      T ~ subtype),
    epidemic_of_season = case_when(
      grepl('imd_spline', name) ~ NA, 
      T ~ epidemic_of_season)
    ) %>% unique()

par_name_df <- epid_pars %>% 
  select(epidemic, season, subtype, year, name, par_name, joining_name) %>% unique() %>% 
  rename(variable = joining_name)

#### LOAD IN MCMC POSTERIOR DATA ####

mcmc_samples_filtered <- data.table(readRDS(.args[8]))

admin_cols <- c('likelihood','chain','iteration','epidemic')
posterior_cols <- colnames(mcmc_samples_filtered)[colnames(mcmc_samples_filtered) %notin% admin_cols]

## only do for every tenth/hundredth/thousandth epidemic to save time
mod_val <- if(nrow(mcmc_samples_filtered) > 100000){1000}else{
  ifelse(nrow(mcmc_samples_filtered) > 10000,100,10)}
mcmc_samples_f_f <- mcmc_samples_filtered[iteration %% mod_val == 0, ]

# make unique IDs across multiple seasons
mcmc_samples_f_f[, chain_it_id := paste0(chain, '_', iteration)]
mcmc_samples_f_f[, c('chain','iteration') := NULL]

admin_cols <- c('likelihood','chain_it_id','epidemic')

# removing any subtype-seasons which don't align with true data
# (i.e. don't care about epidemic 2 in 2025)

post_long <- melt.data.table(mcmc_samples_f_f, 
                             id.vars = admin_cols)[par_name_df, on = c('epidemic','variable')]

# make wide again
post_wide <- dcast.data.table(post_long,
                              likelihood + chain_it_id + epidemic + season + year + subtype ~ name, value.var = 'value')

cat('\nRunning fitted epidemics (', nrow(post_wide)*2/3, ' total): ', sep = '')

fitted_epidemics <- data.table()

for(season_i in unique(post_long$season)){ 
  
  cat('\n\n', season_i)
  
  start_of_season_i <- as.numeric(substr(season_i, 1, 4))
  
  # data to be run
  epids_to_run <- post_wide[season == season_i & !is.na(subtype)]
  
  # how many epidemics (subtypes) in that season
  subtypes_to_run <- unique(epids_to_run$subtype[!is.na(epids_to_run$subtype)])
  n_epids_to_run <- length(subtypes_to_run)
  
  for(subtype_i in subtypes_to_run){
    
    cat('\n', subtype_i, '\n')
    
    vaccinated_pop_seasonal <- vaccinated_data %>% 
      filter(start_of_season == start_of_season_i,
             subtype == subtype_i) %>% 
      arrange(start_of_season, desc(risk_level), imd_quintile, age_grp)
    
    ## should be ordered by IMD then age
    if(vaccinated_pop_seasonal$imd_quintile[2] != 1){warning('vaccinated_data in wrong order')}
    if(vaccinated_pop_seasonal$age_grp[2] != age_labels[2]){warning('vaccinated_data in wrong order')}
    
    # population sizes (done here as may become season-specific)
    pop_stratified <- vaccinated_pop_seasonal$pop 
    pop_vaccinated <- vaccinated_pop_seasonal$vaccinated_population
    VE_INF <- vaccinated_pop_seasonal$VE_INF
    
    tot_pop <- sum(imd_age_pop$pop)
    if(!all.equal(sum(pop_stratified), tot_pop)){warning('pop not adding up')}
    
    epids_to_run_subtype <- epids_to_run[subtype == subtype_i, ]
    
    for(i in 1:nrow(epids_to_run_subtype)){
      
      data <- epids_to_run_subtype[i, ]
      
      init_infected_num <- 10^(data$init_infected)
      init_infected_vec <- (pop_stratified - pop_vaccinated)*init_infected_num/(tot_pop-sum(pop_vaccinated))
      if(!all.equal(sum(init_infected_vec), init_infected_num)){warning('init infected not adding up')}
      
      susceptibility_vec <- fcn_assign_ages(
        data$adult_susceptibility*data$rel_children_susceptibility,
        data$adult_susceptibility,
        data$adult_susceptibility*data$rel_older_adults_susceptibility,
        age_labels
      )
      
      time_series <- run_model(
        pop = pop_stratified,
        I0 = init_infected_vec,
        vacc_cov = pop_vaccinated,
        ve_inf = VE_INF,
        cm = pc_cm,
        trans = data$transmissibility,
        susc = susceptibility_vec,
        lat_per = known_pars$epid_periods[1],
        inf_per = known_pars$epid_periods[2]
      )
      
      time_series <- time_series[, lapply(.SD, sum), 
                                 by = c('t', 'age_grp', 'imd_quintile', 'risk_level', 'vaccinated')]
      
      time_series[, epidemic := data$epidemic]
      
      data_join <- data[, c('likelihood','chain_it_id','epidemic','season','year','subtype')]
      time_series <- time_series[data_join, on = 'epidemic']
      
      fitted_epidemics <- rbind(
        fitted_epidemics,
        time_series
        
      )
      
      if(i %% 10 == 0){cat(i, ', ', sep = '')}
      
    }
  }
  
}

#### IHRs AND IGPRs ####

## column names
primary_names <- paste0('primary_imd_', 1:5)
secondary_names <- paste0('secondary_imd_', 1:5)

## apply imd_spline function
rep_rates <- mcmc_samples_f_f %>% 
  bind_cols(pmap(list(mcmc_samples_f_f$imd_spline_primary_1, mcmc_samples_f_f$imd_spline_primary_2),
                 function(x, y) {
                   result <- imd_spline(c(x, y))
                   setNames(as.list(result), primary_names)
                 }) %>% bind_rows()) %>% 
  bind_cols(pmap(list(mcmc_samples_f_f$imd_spline_secondary_1, mcmc_samples_f_f$imd_spline_secondary_2),
                 function(x, y) {
                   result <- imd_spline(c(x, y))
                   setNames(as.list(result), secondary_names)
                 }) %>% bind_rows())

age_risk_rates <- rep_rates %>% 
  select(contains('rate'), chain_it_id, epidemic) %>%
  pivot_longer(!c(chain_it_id,epidemic)) %>% 
  mutate(care_setting = case_when(substr(name, 1, 1) == 'p' ~ 'primary_rate', T ~ 'secondary_rate'),
         risk_level = case_when(grepl('low', name) ~ 'low', T ~ 'high'),
         broad_age = case_when(grepl('children', name) ~ 'children', 
                               grepl('older_adults', name) ~ 'older_adults', 
                               T ~ 'adults'),
         epidemic_of_season = as.numeric(substr(name, nchar(name), nchar(name)))) %>% 
  select(!name) %>% 
  left_join(subtype_seasons, by = c('epidemic','epidemic_of_season')) %>% 
  filter(!is.na(subtype)) # e.g. removing epidemic 2 of 2025

imd_rates <- rep_rates %>% 
  select(contains('ary_imd'), chain_it_id, epidemic) %>%
  pivot_longer(!c(chain_it_id,epidemic)) %>% 
  mutate(care_setting = case_when(substr(name, 1, 1) == 'p' ~ 'primary_rate', T ~ 'secondary_rate'),
         imd_quintile = as.factor(gsub('primary_imd_|secondary_imd_','',name))) %>% 
  select(!name)

mcmc_surveillance_rates <- full_join(
  age_risk_rates, imd_rates, by = c('chain_it_id','epidemic','care_setting'), relationship = "many-to-many"
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
  group_by(season, risk_level, broad_age, imd_quintile, name, epidemic_of_season, subtype) %>% 
  summarise(median = median(value),
            l = quantile(value, 0.025),
            u = quantile(value, 0.975)) %>% 
  left_join(unknown_pars$care_rates %>% 
              select(!contains('rel_')) %>% 
              left_join(broad_ages, by = 'age_grp') %>% 
              select(!age_grp) %>% unique() %>% 
              rename(primary_rate = gp_rate,
                     secondary_rate = hosp_rate) %>% 
              pivot_longer(c(primary_rate, secondary_rate)) %>% 
              rename(true_value = value) %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('risk_level', 'imd_quintile', 'broad_age', 'name', 'subtype', 'season'))

hc_rates_fitted$broad_age <- factor(
  hc_rates_fitted$broad_age, 
  levels = c('children','adults','older_adults')
)

plot_healthcare_rates <- function(i){
  
  hc_rates_fitted %>% 
    filter(subtype == subtype_seasons$subtype[i],
           season == subtype_seasons$season[i]) %>% 
    mutate(name = gsub('_rate', ' care', name),
           risk_level = paste0(risk_level, ' risk')) %>% 
    ggplot(aes(x = broad_age, group = imd_quintile, col = imd_quintile)) + 
    geom_errorbar(aes(ymin = l, ymax = u), width = 0.4,
                  position = position_dodge(width = 0.4), alpha=1) +
    geom_point(aes(y = median), #shape = 1, 
               position = position_dodge(width = 0.4)) +
    geom_point(aes(y = true_value), shape = 4, size = 3, 
               position = position_dodge(width = 0.4), stroke = 0.8) +
    theme_lg() +
    scale_color_manual(values = imd_quintile_colors) + 
    facet_wrap(risk_level ~ name, scales = 'free', ncol = 2) + 
    # theme(legend.position = 'none') +
    scale_y_continuous(labels = scales::percent, limits = c(0,NA)) +
    labs(x = '', y = 'Healthcare attendance upon infection',
         color = 'IMD quintile',
         title = paste0(subtype_seasons$subtype[i], ', ', subtype_seasons$season[i]))
  
}

hc_plots <- map(.x = 1:nrow(subtype_seasons),
                .f = plot_healthcare_rates)

patchwork::wrap_plots(hc_plots) + plot_layout(guides = 'collect')
ggsave(gsub('data/posterior_epidemics.rds',figure_filename('healthcare_attendance'),.args[length(.args)]),
       width = 20, height = 10)

## add in start_date, reporting rates
fitted_epidemics[, start_of_epidemic := as.Date(paste0('01-09-', year), format = '%d-%m-%Y')]
fitted_epidemics[, date := start_of_epidemic + t] ## add date
fitted_epidemics[, c('t','start_of_epidemic') := NULL]

# aggregate week
fitted_epidemics[, date := last_monday(date)]
fitted_epidemics_agg <- fitted_epidemics[, .(infections = ceiling(sum(infections))),
                                         by = .(date, epidemic, chain_it_id, age_grp, imd_quintile, risk_level,
                                                subtype, vaccinated)]

fitted_epidemics_surv <- copy(fitted_epidemics_agg) ## copy for surveillance data

## sum over subtypes
fitted_epidemics_agg[, c('vaccinated', 'epidemic') := NULL]
fitted_epidemics_agg <- fitted_epidemics[, .(infections = sum(infections)),
                                         by = .(chain_it_id, date, age_grp, imd_quintile, risk_level)]
fitted_epidemics_agg[, chain_it_id := NULL]

#### SAVE EPIDEMICS ####
write_rds(fitted_epidemics_agg, .args[length(.args)])

#### PLOT EPIDEMICS ####

fitted_epidemics_agg_m <- rbind(
  fitted_epidemics_agg[, lapply(.SD, median), by = c('date', 'age_grp', 'imd_quintile', 'risk_level')][, measure := 'median'],
  fitted_epidemics_agg[, lapply(.SD, max), by = c('date', 'age_grp', 'imd_quintile', 'risk_level')][, measure := 'u'],
  fitted_epidemics_agg[, lapply(.SD, min), by = c('date', 'age_grp', 'imd_quintile', 'risk_level')][, measure := 'l'])

fitted_epidemics_agg_l <- dcast(fitted_epidemics_agg_m, 
                                date + age_grp + imd_quintile + risk_level ~ measure,
                                value.var = 'infections')

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
  geom_line(aes(date, median, group=interaction(imd_quintile,risk_level), col=imd_quintile, lty = risk_level)) +
  geom_point(aes(date, infections, group=interaction(imd_quintile,risk_level), col=imd_quintile)) +
  facet_wrap(age_grp~., scales='free') + theme_bw() +
  scale_linetype_manual(values = c(2,1)) + 
  scale_color_manual(values = imd_quintile_colors) +
  scale_fill_manual(values = imd_quintile_colors)
ggsave(gsub('data/posterior_epidemics.rds',figure_filename('posterior_epidemics'),.args[length(.args)]),
       width = 30, height = 20)

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

## CRUDE ATTACK RATE ESTIMATES
## (doing sum of median instead of median of sum, for ease of calculation (preliminary))
fitted_and_obs %>% group_by(age_grp, imd_quintile) %>% 
  summarise(med = sum(median)/length(years), 
            l = sum(l)/length(years),
            u = sum(u)/length(years),
            inf = sum(infections/length(years))) %>% 
  left_join(imd_age_pop %>% mutate(imd_quintile=as.character(imd_quintile)), 
            by = c('age_grp', 'imd_quintile')) %>% 
  ggplot() + 
  geom_bar(aes(x=age_grp, y=med/pop, fill=imd_quintile),
           stat='identity',position='dodge') +
  geom_errorbar(aes(x=age_grp, ymin=l/pop, ymax=u/pop, group=imd_quintile),
                position = position_dodge(width = 0.9), width=0.4, alpha=0.7) +
  theme_bw() + scale_fill_manual(values = imd_quintile_colors) +
  labs(x='', fill='IMD quintile', y='Attack rate')

#### SUBTYPE PROPORTIONS ####

epids_subtype <- copy(fitted_epidemics_surv)

## sum over subtypes
epids_subtype[, epidemic := NULL]
epids_subtype_agg <- epids_subtype[, .(infections = sum(infections)),
                                         by = .(chain_it_id, date, subtype)]

epids_subtype_agg_w <- dcast(epids_subtype_agg, 
                             chain_it_id + date ~ subtype,
                             value.var = 'infections')

# replace NA with 0
epids_subtype_agg_w[is.na(AH1N1), AH1N1 := 0]
epids_subtype_agg_w[is.na(AH3N2), AH3N2 := 0]
epids_subtype_agg_w[is.na(B), B := 0]

epids_subtype_agg_w[, total := AH1N1 + AH3N2 + B]

epids_subtype_agg_l <- melt(epids_subtype_agg_w,
                            id.vars = c('chain_it_id','date','total'))

epids_subtype_agg_l[, modelled_proportion := value/total]

setnames(epids_subtype_agg_l, 'variable', 'subtype')
epids_subtype_agg_l[, c('total','value','chain_it_id') := NULL]

epids_subtype_m <- rbind(
  epids_subtype_agg_l[, lapply(.SD, median), by = c('date','subtype')][, measure := 'median'],
  epids_subtype_agg_l[, lapply(.SD, max), by = c('date','subtype')][, measure := 'u'],
  epids_subtype_agg_l[, lapply(.SD, min), by = c('date','subtype')][, measure := 'l'])

epids_subtype_w <- dcast(epids_subtype_m,
                         date + subtype ~ measure, 
                         value.var = 'modelled_proportion')

## true proportions
true_proportions <- ukhsa_subtype_data %>% rename(date = date_formatted) %>% 
  select(date, subtype, proportion, season, total_tests) %>% 
  group_by(date, season) %>% 
  complete(subtype = unique(ukhsa_subtype_data$subtype), fill = list(proportion = NA,
                                                                     total_tests = NA))

## combine
fitted_and_obs_proportions <- epids_subtype_w %>% 
  left_join(true_proportions, by = c('date','subtype')) 

fitted_and_obs_proportions %>% 
  filter(!is.na(season)) %>% 
  ggplot() + 
  geom_ribbon(aes(date, ymin = l, ymax = u, fill = subtype), alpha = 0.3) +
  geom_line(aes(date, median, col = subtype, group = season), lwd = 0.8) +
  geom_line(aes(date, proportion, col = subtype, group = season), lwd = 0.4, lty = 2) +
  geom_point(aes(date, proportion, col = subtype, size = total_tests)) +
  theme_lg() +
  facet_grid(subtype ~ .) +
  scale_size_continuous(range = c(0.2, 4)) + 
  scale_x_date(breaks = "1 month", labels=date_format("%b\n%y")) +
  scale_color_manual(values = subtype_colors, labels = subtype_names) +
  scale_fill_manual(values = subtype_colors, labels = subtype_names) +
  labs(y = 'Proportion of infections', col = '', fill = '', size = 'Total weekly\nUKHSA tests', x='')

ggsave(gsub('data/posterior_epidemics.rds',figure_filename('subtype_proportions'),.args[length(.args)]),
       width = 16, height = 8)


#### FITTED SURVEILLANCE DATA ####

fitted_surv_dat <- fitted_epidemics_surv %>% 
  left_join(broad_ages, by = 'age_grp') %>% 
  left_join(mcmc_surveillance_rates_w,
            by = c('broad_age','imd_quintile','risk_level','chain_it_id','epidemic','subtype')) %>% 
  left_join(known_pars$proportion_observed %>% mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('age_grp','imd_quintile','risk_level')) %>% 
  left_join(vaccinated_data %>% 
              select(age_grp, imd_quintile, risk_level, start_of_season, 
                     subtype, VE_INF, VE_HOSP) %>% 
              mutate(epidemic = start_of_season - year(fitted_epidemics$date[1]) + 1,
                     imd_quintile = as.factor(imd_quintile)),
            by = c('age_grp','imd_quintile','risk_level','epidemic','subtype')) %>% 
  # scale down IHR for ineffectively vaccinated infections
  mutate(
    ve_hosp_multiplier = case_when(
      !vaccinated ~ 1, vaccinated ~ (1 - VE_HOSP)/(1 - VE_INF)),
    # if VE_HOSP < VE_INF, assume no protection against hospitalisation
    ve_hosp_multiplier = case_when(
      ve_hosp_multiplier > 1 ~ 1, T ~ ve_hosp_multiplier),
    secondary_rate = case_when(
      !vaccinated ~ secondary_rate, vaccinated ~ ve_hosp_multiplier*secondary_rate)) %>% 
  select(!ve_hosp_multiplier) %>% 
  mutate(observed_primary = infections*primary_rate*OS_COVERAGE,
         observed_secondary = infections*secondary_rate*OS_COVERAGE) %>% 
  group_by(chain_it_id, date, age_grp, imd_quintile, risk_level) %>% 
  summarise(observed_primary = sum(observed_primary),
            observed_secondary = sum(observed_secondary)) %>% 
  left_join(surveillance_data %>% ungroup() %>% select(!index) %>% 
              rename(date = week_start) %>% 
              filter(date %in% unique(fitted_epidemics$date)) %>% 
              mutate(imd_quintile = as.factor(imd_quintile)), 
            by = c('date','age_grp','imd_quintile','risk_level')) 

fitted_surv_dat <- data.table(fitted_surv_dat)

fitted_surv_dat[, chain_it_id := NULL]

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
ggsave(gsub('data/posterior_epidemics.rds',figure_filename('primary_data'),.args[length(.args)]),
       width = 16, height = 12)
fitted_secondary_plot
ggsave(gsub('data/posterior_epidemics.rds',figure_filename('secondary_data'),.args[length(.args)]),
       width = 16, height = 12)

