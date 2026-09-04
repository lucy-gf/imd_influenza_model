## CREATE DUMMY DATA (INFECTIONS) ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(readr))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "contact_matrix.rds"),
  file.path("data", "dummy_data", "known_parameters.rds"),
  file.path("data", "dummy_data", "unknown_parameters.rds"),
  file.path("data", "dummy_data", "dummy_infections.rds")
) else commandArgs(trailingOnly = TRUE)
  
source(file.path('scripts','setup', 'base_functions.R'))
source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','seir_model.R'))

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
contact_matrix <- readRDS(.args[2])

## check the right number of rows (no extra defining variables)
if(nrow(contact_matrix) != nrow(contact_matrix %>% select(p_imd_q, c_imd_q, p_age_group, c_age_group) %>% unique())){
  warning('Number of rows in contact matrix wrong')
}

contact_matrix$p_age_group <- factor(contact_matrix$p_age_group,
                                     levels = age_labels)
contact_matrix$c_age_group <- factor(contact_matrix$c_age_group,
                                     levels = age_labels)

## ENSURE ALL POP ORDERS ARE IMD THEN AGE
contact_matrix <- contact_matrix %>% 
  arrange(p_imd_q, c_imd_q,
          p_age_group, c_age_group)

## into matrix
cm <- contact_matrix %>% ungroup() %>% 
  mutate(p_var = paste0(p_imd_q, '_', p_age_group),
         c_var = paste0(c_imd_q, '_', c_age_group)) %>% 
  select(p_var,c_var,n) %>% 
  pivot_wider(names_from = c_var, values_from = n) %>% 
  select(!p_var) %>% as.matrix()

## KNOWN PARAMETERS
known_pars <- readRDS(.args[3])
years <- known_pars$years

vaccinated_data <- known_pars$vaccinated_data

## UNKNOWN PARAMETERS
unknown_pars <- readRDS(.args[4])
subtype_vec <- unique(unknown_pars$epid_parameters$subtype)
subtype_seasons <- unknown_pars$epid_parameters %>% 
  select(year, subtype) %>% unique()

## check R0
cat('\n')
R0_vec <- c(); Reff_vec <- c()
for(i in 1:nrow(subtype_seasons)){
  
    subtype_season_i <- subtype_seasons[i,]
    
    # cat(subtype_season_i$year, ', ', subtype_season_i$subtype, ': ', sep = '')
    
    v_p <- vaccinated_data %>%
      filter(start_of_season == subtype_season_i$year) %>%
      group_by(age_grp, imd_quintile) %>% 
      summarise(effectively_vaccinated_population = sum(effectively_vaccinated_population), 
                pop = sum(pop)) %>% ungroup() %>% 
      mutate(eff_v_p = effectively_vaccinated_population/pop) %>% 
      arrange(imd_quintile, age_grp) 
    
    pars <- unknown_pars$epid_parameters %>% filter(year == subtype_season_i$year,
                                                    subtype == subtype_season_i$subtype)
    
    absolute_susceptibility <- fcn_assign_ages(
      pars$susceptibility*pars$rel_susc_children,
      pars$susceptibility,
      pars$susceptibility*pars$rel_susc_older_adults,
      age_labels
    )
                                 
    if(sum(absolute_susceptibility > 1) > 0){stop('Some susceptibility over 100%')}
    
    R0 <- R0_func(absolute_susceptibility,
                  pars$infectious_period,
                  pars$transmissibility,
                  cm)
    Reff <- R0_func((1 - v_p$eff_v_p)*rep(absolute_susceptibility, 5),
                    pars$infectious_period,
                    pars$transmissibility,
                    cm)
    # cat('R0 = ', round(R0,3), ', ',
    #     'Reff = ', round(Reff,3),
    #     '\n', sep = '')
    
    R0_vec <- c(R0_vec, R0); Reff_vec <- c(Reff_vec, Reff)
  
}

R_dat <- data.table(subtype_seasons %>% mutate(R0 = R0_vec, Reff = Reff_vec))
print(R_dat)

#### EXPAND CONTACT MATRIX ####

## make matrix per capita
per_cap_matrix_45 <- t(t(cm)/imd_age_pop$pop)

## check calculation assumption is correct
# test_m<-matrix(c(1,2,3, 11,12,12), nrow = 2, ncol = 3, byrow = TRUE)
# should_be_111_1164 <- t(t(test_m)/c(1,2,3))

## scale up for both risk groups 
## (assuming per capita contacts are independent of risk level)
pc_cm <- expand_contact_matrix(per_cap_matrix_45)

ng <- nrow(per_cap_matrix_45)
ndim <- nrow(pc_cm)
if(!all.equal(2*ng, ndim)){warning('dimensions not adding up')}

#### RUN EACH SEASON ####

seasonal_seir_outputs <- list()

for(i in 1:nrow(subtype_seasons)){
  
  subtype_season_i <- subtype_seasons[i,]
  
  vaccinated_pop_seasonal <- vaccinated_data %>% 
    filter(start_of_season == subtype_season_i$year,
           subtype == subtype_season_i$subtype) %>% 
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
  
  pars <- unknown_pars$epid_parameters %>% 
    filter(year == subtype_season_i$year,
           subtype == subtype_season_i$subtype)
  
  init_infected_num <- pars$init_infected
  init_infected_vec <- (pop_stratified - pop_vaccinated)*init_infected_num/(tot_pop-sum(pop_vaccinated))
  if(!all.equal(sum(init_infected_vec), init_infected_num)){warning('init infected not adding up')}
  
  susceptibility_vec <- fcn_assign_ages(
    pars$susceptibility*pars$rel_susc_children,
    pars$susceptibility,
    pars$susceptibility*pars$rel_susc_older_adults,
    age_labels
  )
  
  time_series <- run_model(
    pop = pop_stratified,
    I0 = init_infected_vec,
    vacc_cov = pop_vaccinated,
    ve_inf = VE_INF,
    cm = pc_cm,
    trans = pars$transmissibility,
    susc = susceptibility_vec,
    lat_per = pars$latent_period,
    inf_per = pars$infectious_period
  )
  
  # check not rising at end of time
  n_tot_subpops <- nrow(time_series[t == max(t), ])
  n_subpops <- sum(time_series[t == max(t) - 1, ]$infections < time_series[t == max(t), ]$infections)
  if(n_subpops > 0){
    warning(paste0('Epidemic ', paste0(subtype_season_i, collapse = ' '),' still rising in ', n_subpops, '/', n_tot_subpops, ' subpopulations'))
  }
  
  # add population data
  
  time_series <- time_series[vaccinated_pop_seasonal %>% mutate(imd_quintile=as.character(imd_quintile)) %>% 
                               select(age_grp, imd_quintile, risk_level, pop),
                             on = c('age_grp','imd_quintile','risk_level')]
  
  if(F) {
    time_series %>% group_by(t) %>% summarise(inf = sum(infections), pop = sum(pop)) %>% ggplot() + geom_line(aes(t,inf/pop))
  }
  
  n_subpops_min <- sum(time_series[t == max(t), ]$infections> 1)
  if(n_subpops_min > 0){
    warning(paste0('Epidemic ', paste0(subtype_season_i, collapse = ' '),' still above 1 in ', n_subpops_min, '/', n_tot_subpops, ' subpopulations'))
  }
  
  time_series[, start_date := pars$start_date]
  time_series[, subtype := subtype_season_i$subtype]
  
  time_series <- time_series[t %in% 0:365] ## take data from the start of each day
  
  seasonal_seir_outputs[[i]] <- time_series
  
}

seasonal_seir_outputs_names <- paste0(subtype_seasons$year, ' (', subtype_seasons$subtype, ')')
names(seasonal_seir_outputs) <- seasonal_seir_outputs_names

#### PLOT EXAMPLES ####
  
seasonal_seir_outputs[[1]] %>% 
  ggplot() +
  geom_line(aes(t, infections/pop, col = imd_quintile, lty = vaccinated)) +
  scale_color_manual(values = imd_quintile_colors) +
  facet_grid(age_grp ~ risk_level, scales = 'free')

seasonal_seir_outputs[[1]] %>% 
  mutate(age_grp = case_when(
    grepl('18|26|35|50', age_grp) ~ '18-69',
    T ~ age_grp
  )) %>% 
  mutate(age_grp = factor(age_grp, levels = c('0-4','5-11','12-17','18-69','70-79','80+'))) %>% 
  group_by(t, age_grp, imd_quintile) %>% 
  summarise(inf = sum(infections), pop = sum(pop)) %>% 
  ggplot() +
  geom_line(aes(t, 100000*inf/pop, col = imd_quintile), lwd = 0.8) +
  scale_color_manual(values = imd_quintile_colors) +
  facet_wrap(age_grp ~ ., scales = 'free') + theme_bw() +
  labs(x = 'Day of epidemic', y = 'Infections per 100,000', col = 'IMD')

#### FINAL SIZE ####

plot_final_size <- function(k){
  
  seasonal_seir_outputs[[k]] %>% 
    group_by(age_grp, imd_quintile, risk_level, pop) %>% 
    summarise(infections = sum(infections)) %>% 
    ggplot() + 
    geom_bar(aes(x = age_grp, y = 100*infections/pop, 
                 fill = imd_quintile),
             stat = 'identity', position = 'dodge') + 
    theme_bw() + 
    scale_fill_manual(values = imd_quintile_colors) +
    facet_grid(.~risk_level) +
    labs(x = 'Age group', y = 'Final size (%)', fill = 'IMD quintile') +
    ggtitle(names(seasonal_seir_outputs)[k])

}

final_size_plots <- map(.x = 1:length(seasonal_seir_outputs), .f = plot_final_size)

patchwork::wrap_plots(final_size_plots, nrow = 3)

R_dat <- R_dat[, final_size := (rbindlist(seasonal_seir_outputs, idcol = "id") %>% 
                 group_by(id, age_grp, imd_quintile, risk_level, pop) %>% 
                 summarise(infections = sum(infections)) %>% 
                 group_by(id) %>% summarise(fs = sum(infections)/sum(pop)))$fs]

#### FLU WEEKLY PROPORTIONS ####

plot_weekly_props <- function(year_i){
  
  dat_vec <- R_dat[year==year_i]
  
  rbindlist(seasonal_seir_outputs, idcol = "id") %>% 
    filter(substr(id, 1, 4) == as.character(year_i)) %>% 
    group_by(t, age_grp, subtype) %>% 
    summarise(infections = sum(infections), 
              pop = sum(pop)) %>% 
    ggplot() + 
    geom_bar(aes(x = t, y = infections/pop, 
                 fill = subtype),
             stat = 'identity', position = 'stack', width = 1) + 
    theme_bw() + 
    scale_fill_manual(values = subtype_colors, labels = subtype_names) +
    facet_wrap(.~age_grp) +
    labs(x = 'Day of epidemic', y = 'Infected proportion', fill = '',
         title = year_i,
         subtitle = paste0(
           dat_vec$subtype, ': R0 ',
           round(dat_vec$R0, 2), ', Reff ',
           round(dat_vec$Reff, 2), ', AR ',
           100*round(dat_vec$final_size, 4), '%', collapse = '\n'))

 }

weekly_prop_plots <- map(.x = years, .f = plot_weekly_props)

patchwork::wrap_plots(weekly_prop_plots, nrow = 3)

ggsave(file.path("output", "figures", "dummy_infections", "dummy_infections.png"),
       height = 13, width = 9)

#### SAVE DUMMY DATA ####

saveRDS(seasonal_seir_outputs, .args[5])


