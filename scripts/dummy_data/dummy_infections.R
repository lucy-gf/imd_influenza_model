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

## check R0
cat('\n')
R0_vec <- c(); Reff_vec <- c()
for(i in 1:length(years)){
  for(j in c("A", "B")){
    cat(years[i], ', ', j, ': ', sep = '')
    
    v_p <- vaccinated_data %>%
      filter(start_of_season == years[i]) %>%
      group_by(age_grp, imd_quintile) %>% 
      summarise(effectively_vaccinated_population = sum(effectively_vaccinated_population), 
                pop = sum(pop)) %>% ungroup() %>% 
      mutate(eff_v_p = effectively_vaccinated_population/pop) %>% 
      arrange(imd_quintile, age_grp) 
    
    pars <- unknown_pars[[paste0('epid_parameters_s', i, '_', j)]]
    R0 <- R0_func(pars$susceptibility,
                  pars$infectious_period,
                  pars$transmissibility,
                  cm)
    Reff <- R0_func((1 - v_p$eff_v_p)*rep(pars$susceptibility, 5),
                    pars$infectious_period,
                    pars$transmissibility,
                    cm)
    cat('R0 = ', round(R0,3), ', ',
        'Reff = ', round(Reff,3),
        '\n', sep = '')
    
    R0_vec <- c(R0_vec, R0); Reff_vec <- c(Reff_vec, Reff)
  }
}

R_dat <- data.table(
  CJ(year = years, strain = c("A", "B")),
  R0 = R0_vec, Reff = Reff_vec
)

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

for(i in 1:length(years)){
  
  for(j in unique(vaccinated_data$strain)){
    
    list_number <- 2*i + (j == 'B') - 1
    
    vaccinated_pop_seasonal <- vaccinated_data %>% 
      filter(start_of_season == years[i],
             strain == j) %>% 
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
    
    pars <- unknown_pars[[paste0('epid_parameters_s', i, '_', j)]]
    
    init_infected_num <- pars$init_infected
    init_infected_vec <- (pop_stratified - pop_vaccinated)*init_infected_num/(tot_pop-sum(pop_vaccinated))
    if(!all.equal(sum(init_infected_vec), init_infected_num)){warning('init infected not adding up')}
    
    time_series <- run_model(
      pop = pop_stratified,
      I0 = init_infected_vec,
      vacc_cov = pop_vaccinated,
      ve_inf = VE_INF,
      cm = pc_cm,
      trans = pars$transmissibility,
      susc = pars$susceptibility,
      lat_per = pars$latent_period,
      inf_per = pars$infectious_period
    ) 
    
    time_series <- time_series[vaccinated_pop_seasonal %>% mutate(imd_quintile=as.character(imd_quintile)) %>% 
                                 select(age_grp, imd_quintile, risk_level, pop),
                               on = c('age_grp','imd_quintile','risk_level')]
    
    time_series[, start_date := pars$start_date]
    time_series[, strain := j]
    
    time_series <- time_series[t %in% 0:365] ## take data from the start of each day
    
    seasonal_seir_outputs[[list_number]] <- time_series
    
  }
  
}

seasonal_seir_outputs_names <- c()
for(k in 1:length(seasonal_seir_outputs)){
  seasonal_seir_outputs_names <- c(seasonal_seir_outputs_names,
                                   paste0(year(seasonal_seir_outputs[[k]]$start_date[1]),
                                          ' (', seasonal_seir_outputs[[k]]$strain[1], ')'))
}
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

#### FLU WEEKLY PROPORTIONS ####

plot_weekly_props <- function(year_i){
  
  rbindlist(seasonal_seir_outputs, idcol = "id") %>% 
    filter(substr(id, 1, 4) == as.character(year_i)) %>% 
    group_by(t, age_grp, strain) %>% 
    summarise(infections = sum(infections), 
              pop = sum(pop)) %>% 
    ggplot() + 
    geom_bar(aes(x = t, y = infections/pop, 
                 fill = strain),
             stat = 'identity', position = 'stack', width = 1) + 
    theme_bw() + 
    scale_fill_manual(values = strain_colors, labels = strain_names) +
    facet_wrap(.~age_grp) +
    labs(x = 'Day of epidemic', y = 'Infected proportion', fill = '',
         title = year_i,
         subtitle = paste0(
           'R0 A = ', round(R_dat[year==year_i]$R0[1], 2),
           ', R0 B = ', round(R_dat[year==year_i]$R0[2], 2),
           ', Reff A = ', round(R_dat[year==year_i]$Reff[1], 2),
           ', Reff B = ', round(R_dat[year==year_i]$Reff[2], 2)))

 }

weekly_prop_plots <- map(.x = years, .f = plot_weekly_props)

patchwork::wrap_plots(weekly_prop_plots, nrow = 3)

ggsave(file.path("output", "figures", "dummy_infections", "dummy_infections.png"),
       height = 12, width = 9)

plot_weekly_props_FILT <- function(year_i){
  
  rbindlist(seasonal_seir_outputs, idcol = "id") %>% 
    filter(substr(id, 1, 4) == as.character(year_i)) %>% 
    group_by(t, age_grp, strain) %>% 
    summarise(infections = sum(infections), 
              pop = sum(pop)) %>% 
    ggplot() + 
    geom_bar(aes(x = t, y = infections/pop, 
                 fill = strain),
             stat = 'identity', position = 'stack', width = 1) + 
    theme_bw() + 
    scale_fill_manual(values = strain_colors, labels = strain_names) +
    facet_wrap(.~age_grp) +
    labs(x = 'Day of epidemic', y = 'Infected proportion', fill = '',
         # title = year_i,
         subtitle = paste0(
           # 'R0 A = ', round(R_dat[year==year_i]$R0[1], 2),
           # ', R0 B = ', round(R_dat[year==year_i]$R0[2], 2),
           'Reff A = ', round(R_dat[year==year_i]$Reff[1], 2),
           ', Reff B = ', round(R_dat[year==year_i]$Reff[2], 2)))
  
}

weekly_prop_plots <- map(.x = c(2023, 2025), .f = plot_weekly_props_FILT)

patchwork::wrap_plots(weekly_prop_plots, nrow = 2)

ggsave(file.path("output", "figures", "dummy_infections", "dummy_infections_FILT.png"),
       height = 12, width = 9)

#### SAVE DUMMY DATA ####

saveRDS(seasonal_seir_outputs, .args[5])


