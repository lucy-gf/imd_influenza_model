## UNKNOWN PARAMETERS ##

#### SETUP ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(patchwork))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "imd_age_pop.rds"),
  file.path("data", "inputs", "subtype_years.rds"),
  file.path("data", "dummy_data", "unknown_parameters.rds")
) else commandArgs(trailingOnly = TRUE)

set.seed(60)

source(file.path("scripts","setup","colors.R"))
source(file.path('scripts','setup','base_functions.R'))

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

imd_age_pop <- imd_age_pop %>% 
  arrange(imd_quintile, age_grp)

broad_ages_susc <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages(
    'children', 
    'adults',
    'older_adults',
    age_labels
    )
)
broad_ages_care_rates <- data.table(
  age_grp = age_labels, 
  broad_age = fcn_assign_ages(
    'children', 
    'adults',
    'older_adults',
    age_labels
  )
)


## number of years of data
years <- 2023:2025 # 2023-24 to 2025-26
subtype_vec <- c('AH1N1','AH3N2','B')
subtype_years <- read_rds(.args[2])
subtype_years <- subtype_years %>% 
  mutate(subtype = convert_to_subtype(subtype),
         year = as.numeric(substr(season, 1, 4)))

#### EPIDEMIOLOGICAL PARAMETERS ####

epid_periods <- c(2, 3) # latent and infectious periods

## adults' susceptibility, transmissibility
epid_pars <- subtype_years %>% 
  select(year, subtype) %>% unique() %>% 
  mutate(susceptibility = rnorm(n = nrow(subtype_years), mean = 0.4, sd = 0.005),
         transmissibility = rnorm(n = nrow(subtype_years), mean = 0.12, sd = 0.001)) %>% 
  mutate(susceptibility = case_when(subtype == 'B' ~ 0.75*susceptibility, T ~ susceptibility))

## relative susceptibility (relative to adults' susceptibility)
rel_susceptibility <- data.table(cross_join(
  subtype_years %>% 
    select(year, subtype) %>% unique(),
  data.table(broad_age = c('children','older_adults'),
             rel_susceptibility = c(1.5, 1.1))
))

rel_susceptibility[year == 2023 & broad_age == 'children', rel_susceptibility := rel_susceptibility*(1.01)]
rel_susceptibility[year == 2023 & broad_age == 'older_adults', rel_susceptibility := rel_susceptibility*(0.96)]
rel_susceptibility[year == 2024 & broad_age == 'older_adults', rel_susceptibility := rel_susceptibility*(1.02)]
rel_susceptibility[year == 2025 & broad_age == 'children', rel_susceptibility := rel_susceptibility*(0.98)]
rel_susceptibility[subtype == 'B' & broad_age == 'children', rel_susceptibility := rel_susceptibility*(1.3)]
rel_susceptibility[subtype == 'B' & broad_age == 'older_adults', rel_susceptibility := rel_susceptibility*(0.5)]
rel_susceptibility[subtype == 'AH3N2' & broad_age == 'children', rel_susceptibility := rel_susceptibility*(1.01)]
rel_susceptibility[subtype == 'AH3N2' & broad_age == 'older_adults', rel_susceptibility := rel_susceptibility*(0.99)]
if(nrow(rel_susceptibility[rel_susceptibility<0])>0){stop('Negative rel_susceptibility')}

setorder(rel_susceptibility, year, subtype)

rel_susceptibility %>% 
  ggplot() + 
  geom_bar(aes(x = broad_age, y = rel_susceptibility, fill = subtype), 
             position = 'dodge', stat = 'identity', width = 1) +
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_fill_manual(values = subtype_colors) +
  theme_bw() + labs(x = 'Age group', 
                    y = 'Relative susceptibility') +
  facet_grid(. ~ year) +
  theme(text = element_text(size = 14))

#### MAKE INTO DATA TABLE ####

epid_parameters <- subtype_years %>% select(year, subtype) %>% 
  left_join(epid_pars, by = c('year','subtype')) %>% 
  mutate(latent_period = epid_periods[1],
         infectious_period = epid_periods[2],
         start_date = as.Date(paste0('01-09-', year), "%d-%m-%Y"),
         init_infected = floor(rnorm(n = nrow(subtype_years), mean = 300, sd = 10))) %>% 
  mutate(init_infected = case_when(subtype == 'B' ~ 0.5*init_infected, T ~ init_infected)) %>% 
  left_join(rel_susceptibility %>% mutate(broad_age = paste0('rel_susc_', broad_age)) %>% 
              pivot_wider(names_from = broad_age, values_from = rel_susceptibility),
            by = c('year','subtype'))

#### REPORTING RATES ####

# TODO for now these are the same in each season and for each subtype

gp_rate <- c(1, 0.1, 2, 5, 3, 6)/20
hosp_rate <- c(2, 0.1, 2, 8, 3, 15)/100
## LOW RISK CHILDREN, ADULTS, OLDER ADULTS, 
## HIGH RISK CHILDREN, ADULTS, OLDER ADULTS

care_rate_df <- data.frame(
  broad_age = rep(unique(broad_ages_care_rates$broad_age), 2),
  risk_level = rep(c('low','high'), each = 3),
  gp_rate = gp_rate,
  hosp_rate = hosp_rate
)

care_rate_age_df <- data.table(cross_join(
  care_rate_df,
  broad_ages_care_rates))[broad_age.x == broad_age.y,]
care_rate_age_df[, c('broad_age.x','broad_age.y') := NULL]

imd_spline_pars <- data.table(
  primary = c(-0.1, 0.12),
  secondary = c(-0.07, 0.05)
)
rel_imd_rep_rates <- data.frame(imd_quintile = 1:5,
                                rel_primary_rates = imd_spline(imd_spline_pars$primary),
                                rel_secondary_rates = imd_spline(imd_spline_pars$secondary))

care_rate_imd_df <- cross_join(
  care_rate_age_df,
  rel_imd_rep_rates
  ) %>% 
  mutate(gp_rate = gp_rate*rel_primary_rates,
         hosp_rate = hosp_rate*rel_secondary_rates)

care_rate_imd_df$age_grp <- factor(care_rate_imd_df$age_grp, levels = age_labels) 

ratep1 <- care_rate_imd_df %>% 
  pivot_longer(c(gp_rate, hosp_rate)) %>% 
  ggplot() + 
  geom_line(aes(age_grp, value, group = imd_quintile, 
                col = as.factor(imd_quintile)), lwd = 0.8) +
  geom_point(aes(x = age_grp, group = imd_quintile, 
                 y = value), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = imd_quintile, 
                 col = as.factor(imd_quintile), y = value), stroke=1.5, size = 3, shape = 1) +
  theme_bw() + ylim(c(0,NA)) + facet_grid(name ~ risk_level, scales = 'free') + 
  scale_color_manual(values = imd_quintile_colors) +
  labs(y='Rate', col = 'IMD quintile'); ratep1

ratep2 <- rel_imd_rep_rates %>% 
  pivot_longer(!imd_quintile) %>%
  ggplot() + 
  geom_line(aes(imd_quintile, value, group = name, lty = name), lwd = 0.8) +
  geom_point(aes(imd_quintile, value, group = name), 
             color = 'white', size = 3) +
  geom_point(aes(imd_quintile, value, group = name), 
             shape = 1, stroke = 2, size = 3) +
  theme_bw() + ylim(c(0.8,1.2)) +
  scale_linetype_manual(labels = c("Primary care", "Secondary care"),
                        values = c(1,2)) +
  theme(text = element_text(size = 12)) +
  labs(y='Relative rate (baseline = IMD 3)', lty = 'Care setting',
       x = 'IMD quintile'); ratep2

ratep1 + ratep2 + plot_layout(nrow = 1, widths = c(2,1))
ggsave(file.path('output','figures','dummy_infections','reporting_rates.png'), width = 16, height = 7)

#### MAKE INTO LIST ####

unknown_pars <- list(
  epid_parameters = epid_parameters,
  care_rates = care_rate_imd_df,
  imd_spline_pars = imd_spline_pars,
  primary_care_rates = gp_rate,
  secondary_care_rates = hosp_rate
)

## plot all epi pars

plot_epi_par <- function(par_name){
  
  epid_parameters %>% 
    ggplot() + 
    geom_bar(aes(x = year, y = !!sym(par_name), fill = subtype),
               width = 0.8, position = 'dodge', stat = 'identity') +
    scale_fill_manual(values = subtype_colors, labels = subtype_names) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c('2023/24', '2024/25', '2025/26')) +
    theme_bw()
  
}

epi_plots <- map(.x = c('transmissibility', 'susceptibility', 
                        'rel_susc_children', 'rel_susc_older_adults',
                        'init_infected'),
    .f = plot_epi_par)

patchwork::wrap_plots(epi_plots, guides = 'collect')
ggsave(file.path("output","figures","dummy_infections","epid_parameters.png"),
       width = 8, height = 7)

#### SAVE UNKNOWN PARAMETERS ####

saveRDS(unknown_pars, .args[3])

