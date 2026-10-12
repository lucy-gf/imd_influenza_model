## PLOTTING UKHSA PRIORS ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(patchwork))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(viridis))
suppressMessages(require(readr))
suppressMessages(require(readODS))
options(dplyr.summarise.inform = FALSE) 
options(scipen = 9999)

.args <- if (interactive()) c(
  file.path("data", "dummy_data", "dummy_infections.rds"),
  file.path("output", "figures", "exploration", "supervisor_slides", "demo_epids.png")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))

## read in epidemic data
dummy_epids <- readRDS(.args[1])

epidA <- dummy_epids[[1]] %>% 
  filter(age_grp == '0-4', imd_quintile == 1, 
         risk_level == 'low', vaccinated == FALSE)

epid_combined <- epidA %>% left_join(epidA %>% 
                      mutate(t = t + 65,
                             B_infections = infections/6) %>% 
                      select(t, B_infections), 
                    by = 't') %>% 
  mutate(B_infections = case_when(
    is.na(B_infections) ~ 0, T ~ B_infections
  )) 

epid_combined %>% 
  ggplot() + 
  geom_ribbon(aes(t, ymin = 0, ymax = infections), lwd = 1, alpha = 0.4, fill = '#F15156') +
  geom_line(aes(t, infections), lwd = 0.8, col = '#F15156') +
  theme_bw() + labs(x = 'time', y = 'infections') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub('epids','epid',.args[2]), width = 3.5, height = 3)

epid_combined %>% 
  ggplot() + 
  geom_ribbon(aes(t, ymin = 0, ymax = infections), lwd = 1, alpha = 0.4, fill = flu_strain_colors['flu_a']) +
  geom_line(aes(t, infections), lwd = 0.8, col = flu_strain_colors['flu_a']) +
  geom_ribbon(aes(t, ymin = 0, ymax = B_infections), lwd = 1, alpha = 0.4, fill = flu_strain_colors['flu_b']) +
  geom_line(aes(t, B_infections), lwd = 0.8, col = flu_strain_colors['flu_b']) +
  theme_bw() + labs(x = 'time', y = 'infections') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(.args[2], width = 3.5, height = 3)

epid_combined %>% 
  ggplot() + 
  geom_bar(data = epid_combined %>% select(t, infections, B_infections) %>% 
             pivot_longer(!t) %>% mutate(name = case_when(name == 'infections' ~ 'flu_a', T ~ 'flu_b')),
           aes(t, value, fill = name), width=1, alpha = 0.4, 
           stat = 'identity', position = 'stack') +
  geom_line(aes(t, infections + B_infections), lwd = 0.8) +
  theme_bw() + labs(x = 'time', y = 'infections') + 
  scale_fill_manual(values = flu_strain_colors) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + guides(fill="none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub(".png", "_summed.png", .args[2]), width = 3.5, height = 3)

epid_combined %>% 
  ggplot() + 
  geom_bar(data = epid_combined %>% select(t, infections, B_infections) %>% 
             pivot_longer(!t) %>% mutate(name = case_when(name == 'infections' ~ 'flu_a', T ~ 'flu_b')),
           aes(t, value, fill = name), width=1, alpha = 0.4, 
           stat = 'identity', position = 'stack') +
  geom_line(aes(t, infections + B_infections), lwd = 0.8) +
  geom_line(aes(t*1.5 - 48, infections*1.03), lwd = 0.8, lty = 2) +
  theme_bw() + labs(x = 'time', y = 'infections') + 
  scale_fill_manual(values = flu_strain_colors) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + guides(fill="none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub(".png", "_badly_fitted.png", .args[2]), width = 3.5, height = 3)

#### MULTIPLE YEARS ####

epidH1N1_23 <- epidA %>% mutate(subtype = 'AH1N1') %>% 
  select(t, infections, start_date, subtype)
epidH3N2_23 <- epidH1N1_23 %>% mutate(subtype = 'AH3N2',
                                      infections = 0.6*infections, t = t-10) 
epidH1N1_24 <- epidH1N1_23 %>% mutate(start_date = start_date + 365)
epidB_24 <- epidH1N1_23 %>% 
  mutate(t = t + 65,
         start_date = start_date + 365,
         infections = infections/6,
         subtype = 'B')
epidH3N2_25 <- epidH1N1_23 %>% mutate(start_date = start_date + 2*365,
                                      subtype = 'AH3N2')

mult_epid <- rbind(epidH1N1_23,
                   epidH3N2_23,
                   epidH1N1_24,
                   epidB_24,
                   epidH3N2_25) %>% 
  mutate(date = start_date + t)

mult_epid %>% 
  group_by(date) %>% mutate(total = sum(infections),
                            min = min(infections)) %>% 
  group_by(date, subtype) %>% mutate(line = case_when(
    subtype == 'AH1N1' ~ total, T ~ infections
  )) %>% 
  ggplot() + 
  geom_bar(aes(date, infections, fill = subtype),
           position = 'stack', stat = 'identity', width = 1) + 
  geom_line(aes(date, line, col = subtype), lwd = 0.5) +
  theme_bw() + labs(x = 'time', y = 'infections', fill = '', col = '') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  scale_fill_manual(values = subtype_colors, labels = subtype_names) + 
  scale_color_manual(values = subtype_colors, labels = subtype_names) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub(".png", "_multiyear.png", .args[2]), width = 6, height = 3)

mult_epid %>% 
  group_by(date) %>% mutate(total = sum(infections),
                            min = min(infections)) %>% 
  group_by(date, subtype) %>% mutate(line = case_when(
    subtype == 'AH1N1' ~ total, T ~ infections
  )) %>% 
  ggplot() + 
  geom_bar(aes(date, infections, fill = subtype),
           position = 'stack', stat = 'identity', width = 1) + 
  geom_line(aes(date, line, col = subtype), lwd = 0.5) +
  geom_vline(aes(xintercept = as.Date('15-07-2024', format = '%d-%m-%Y')), lty = 2) + 
  geom_vline(aes(xintercept = as.Date('15-07-2025', format = '%d-%m-%Y')), lty = 2) + 
  theme_bw() + labs(x = 'time', y = 'infections', fill = '', col = '') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  scale_fill_manual(values = subtype_colors, labels = subtype_names) + 
  scale_color_manual(values = subtype_colors, labels = subtype_names) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub(".png", "_multiyear_lines.png", .args[2]), width = 6, height = 3)

mult_epid %>% 
  group_by(subtype) %>% 
  mutate(prim = lag(infections, default = 0),
         sec = 0.3*lag(infections, n=2, default = 0)) %>% 
  ungroup() %>% 
  pivot_longer(!c(t, infections, start_date, subtype, date)) %>% 
  ggplot() + 
  geom_line(aes(date, value, lty = name, col = subtype), lwd = 1) +
  theme_bw() + labs(x = 'time', y = 'reported cases', col = '', lty = '') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 2000)) +
  scale_linetype_manual(values = c(2,3), labels = c('Primary','Secondary')) + 
  scale_color_manual(values = subtype_colors, labels = subtype_names) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub(".png", "_surveillance.png", .args[2]), width = 6, height = 3)

mult_epid %>% 
  group_by(subtype) %>% 
  mutate(prim = lag(infections, default = 0),
         sec = 0.3*lag(infections, n=2, default = 0)) %>% 
  ungroup() %>% 
  pivot_longer(!c(t, infections, start_date, subtype, date)) %>% 
  group_by(date, name) %>% summarise(value = sum(value)) %>% 
  ggplot() + 
  geom_line(aes(date, value, lty = name), lwd = 1) +
  theme_bw() + labs(x = 'time', y = 'reported cases', lty = '') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 3300)) +
  scale_linetype_manual(values = c(2,3), labels = c('Primary','Secondary')) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(gsub(".png", "_surveillance_agnostic.png", .args[2]), width = 6, height = 3)

p <- mult_epid %>% 
  group_by(subtype) %>% 
  mutate(infections = case_when(subtype == 'AH1N1' ~ 0.6*lag(infections, default = 0), T ~ 2*infections)) %>%
  group_by(date) %>% mutate(total = sum(infections),
                            min = min(infections)) %>% 
  group_by(date, subtype) %>% mutate(line = case_when(
    subtype == 'AH1N1' ~ total, T ~ infections
  )) %>% 
  filter(date > as.Date('01-07-2024',format = '%d-%m-%Y'),
         date < as.Date('01-07-2025',format = '%d-%m-%Y')) %>% 
  ggplot() + 
  geom_bar(aes(date, infections, fill = subtype),
           position = 'stack', stat = 'identity', width = 1) +
  geom_line(aes(date, line, col = subtype), lwd = 0.5) +
  theme_bw() + labs(x = 'time', y = 'infections', fill = '', col = '') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1200)) +
  scale_fill_manual(values = subtype_colors, labels = subtype_names) + 
  scale_color_manual(values = subtype_colors, labels = subtype_names) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()); p

ggsave(gsub("epids.png", "_MCMC_epids.png", .args[2]), width = 3, height = 2)

p1 <- mult_epid %>% 
  group_by(subtype) %>% 
  mutate(infections = case_when(subtype == 'AH1N1' ~ 0.6*lag(infections, default = 0), T ~ 2*infections)) %>%
  filter(date > as.Date('01-07-2024',format = '%d-%m-%Y'),
         date < as.Date('01-07-2025',format = '%d-%m-%Y')) %>% 
  ggplot() + 
  geom_bar(aes(date, infections, fill = subtype),
           position = 'fill', stat = 'identity', width = 1) +
  theme_bw() + labs(x = 'time', y = 'proportion of\ninfections', fill = '') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = subtype_colors, labels = subtype_names) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()); p1

p1b <- mult_epid %>% 
  group_by(subtype) %>% 
  mutate(sim_infections = case_when(subtype == 'AH1N1' ~ 0.6*lag(infections, default = 0), T ~ 2*infections)) %>%
  filter(date > as.Date('01-07-2024',format = '%d-%m-%Y'),
         date < as.Date('01-07-2025',format = '%d-%m-%Y')) %>% 
  group_by(date) %>% mutate(sim_prop = sim_infections/sum(sim_infections),
                            prop = infections/sum(infections)) %>%
  filter(subtype == 'AH1N1') %>% select(!c(infections, sim_infections)) %>% 
  pivot_longer(!c(t,start_date,subtype,date)) %>% mutate(sim = grepl('sim',name)) %>% 
  drop_na() %>% 
  ggplot() + 
  geom_line(aes(date, value, col = sim), lwd = 1) +
  theme_bw() + labs(x = 'time', y = 'proportion of\ninfections H1N1', col = '') + 
  scale_color_manual(values = c('#6247AA', '#FCA311'), labels = c('True data','Simulated')) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()); p1b

p2 <- mult_epid %>% 
  group_by(subtype) %>% 
  mutate(infections = case_when(subtype == 'AH1N1' ~ 0.6*lag(infections, default = 0), T ~ 2*infections)) %>%
  filter(date > as.Date('01-07-2024',format = '%d-%m-%Y'),
         date < as.Date('01-07-2025',format = '%d-%m-%Y')) %>% 
  group_by(subtype) %>% 
  mutate(prim = lag(infections, default = 0),
         sec = 0.3*lag(infections, n=2, default = 0)) %>% 
  ungroup() %>% 
  pivot_longer(!c(t, infections, start_date, subtype, date)) %>% 
  group_by(date, name) %>% summarise(value = sum(value)) %>% 
  ggplot() + 
  geom_line(aes(date, value, lty = name), lwd = 1) +
  theme_bw() + labs(x = 'time', y = '\nexpected cases', lty = '') + 
  scale_linetype_manual(values = c(2,3), labels = c('Primary','Secondary')) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1800)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()); p2

p2b <- mult_epid %>% 
  group_by(subtype) %>% 
  mutate(lag_infections = case_when(subtype == 'AH1N1' ~ 0.6*lag(infections, default = 0), T ~ 2*infections)) %>%
  filter(date > as.Date('01-07-2024',format = '%d-%m-%Y'),
         date < as.Date('01-07-2025',format = '%d-%m-%Y')) %>% 
  group_by(subtype) %>% 
  mutate(prim = lag(infections, default = 0),
         sec = 0.3*lag(infections, n=2, default = 0),
         lag_prim = lag(lag_infections, default = 0),
         lag_sec = 0.3*lag(lag_infections, n=2, default = 0)) %>% 
  ungroup() %>% 
  pivot_longer(!c(t, infections, lag_infections, start_date, subtype, date)) %>% 
  group_by(date, name) %>% summarise(value = sum(value)) %>% 
  mutate(lag = grepl('lag', name),
         name = gsub('lag_', '', name)) %>% 
  ggplot() + 
  geom_line(aes(date, value, lty = name, col = lag), lwd = 1) +
  theme_bw() + labs(x = 'time', y = '\ncases', lty = '', col = '') + 
  scale_linetype_manual(values = c(2,3), labels = c('Primary','Secondary')) +
  scale_color_manual(values = c('#6247AA', '#FCA311'), labels = c('True data','Simulated')) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1800)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()); p2b

p1 + p2 + plot_layout(nrow = 2, heights = c(1,2))
ggsave(gsub("epids.png", "_MCMC_surveillance.png", .args[2]), width = 5, height = 3.5)
p1b + p2b + plot_layout(nrow = 2, heights = c(1,2))
ggsave(gsub("epids.png", "_MCMC_surveillance_2.png", .args[2]), width = 4.8, height = 3.5)


#### PRIOR DISTRIBUTIONS ####

beta_pars <- function(mean, concentration) {
  c(a = mean * concentration, b = (1 - mean) * concentration)
}

prim_beta  <- beta_pars(0.02, 200) 
sec_beta   <- beta_pars(0.005, 200)
shape1_susc <- 2
shape2_susc <- 2
reff_mean <- 2
reff_sd   <- 0.4
imd_mean <- 0
imd_sd <- 0.5

min_trans <- 0; max_trans <- 1
min_susc <- 0; max_susc <- 1
min_rel_susc <- 1/10; max_rel_susc <- 5
min_reff <- 1; max_reff <- 4
min_log_init_inf <- 0; max_log_init_inf <- 4
min_reporting <- 0; max_reporting <- 1
min_spline <- -log(5); max_spline <- log(5) # equivalent to IMD ratios at most

## data
nsamps <- 10000

prior_values <- data.table(
  reff_values = seq(min_reff, max_reff, length.out = nsamps),
  reff_samples = rnorm(nsamps, mean = reff_mean, sd = reff_sd),
  reff_distr = dnorm(seq(min_reff, max_reff, length.out = nsamps), mean = reff_mean, sd = reff_sd),
  susc_values = seq(min_susc, max_susc, length.out = nsamps),
  susc_distr = dbeta(seq(min_susc, max_susc, length.out = nsamps), shape1_susc, shape2_susc),
  c_susc_samples = rbeta(nsamps, shape1_susc, shape2_susc),
  a_susc_samples = rbeta(nsamps, shape1_susc, shape2_susc),
  oa_susc_samples = rbeta(nsamps, shape1_susc, shape2_susc),
  imd_spline_values = seq(min_spline, max_spline, length.out = nsamps),
  imd_spline_distr = dnorm(seq(min_spline, max_spline, length.out = nsamps), mean = imd_mean, sd = imd_sd),
  imd_spline_samples = rnorm(nsamps, mean = imd_mean, sd = imd_sd),
  init_inf_values = seq(min_log_init_inf, max_log_init_inf, length.out = nsamps), 
  init_inf_distr = dunif(seq(min_log_init_inf, max_log_init_inf, length.out = nsamps), min_log_init_inf, max_log_init_inf),
  init_inf_samples = runif(nsamps, min_log_init_inf, max_log_init_inf),
  reporting_values = seq(min_reporting, max_reporting, length.out = nsamps),
  prim_rates = dbeta(seq(min_reporting, max_reporting, length.out = nsamps), shape1 = prim_beta['a'], shape2 = prim_beta['b']),
  sec_rates = dbeta(seq(min_reporting, max_reporting, length.out = nsamps), shape1 = sec_beta['a'], shape2 = sec_beta['b']),
  transmissibility = 0,
  imd_1_value = 0
)

pb <- txtProgressBar(min = 1, max = nrow(prior_values), style = 3)

for(i in 1:nrow(prior_values)){

  prior_values$transmissibility[i] <- R0_func(
    susceptibility    = fcn_assign_ages(prior_values$c_susc_samples[i],
                                        prior_values$a_susc_samples[i], 
                                        prior_values$oa_susc_samples[i],
                                        age_labels),
    inf_period        = epid_periods[2],
    beta_in           = 1,
    cm_in             = cm_input,
    per_capita        = TRUE,
    population_vector = demography_input$population,
    R0assumed         = prior_values$reff_samples[i],
    return_beta       = TRUE
  )
  
  prior_values$imd_1_value[i] <- exp(prior_values$imd_spline_samples[i])
  
  setTxtProgressBar(pb, i)
  
}
close(pb)

prior_values <- prior_values %>% 
  filter(transmissibility >= min_trans & transmissibility <= max_trans)

reff_plot <- prior_values %>% 
  ggplot() +
  geom_ribbon(aes(x = reff_values, ymax = reff_distr, ymin = 0), fill = '#F75C03',
              alpha = 0.4) +
  geom_line(aes(x = reff_values, y = reff_distr), col = '#F75C03',
               lwd = 0.8) +
  # geom_vline(xintercept = min_reff, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = max_reff, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = 'reff') + 
  theme_lg() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); reff_plot

susc_plot <- prior_values %>%
  ggplot() +
  geom_ribbon(aes(x = susc_values, ymax = susc_distr, ymin = 0), fill = '#820263',
              alpha = 0.4) +
  geom_line(aes(x = susc_values, y = susc_distr), col = '#820263',
            lwd = 0.8) +
  # geom_vline(xintercept = min_trans, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = max_trans, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = 'susceptibility') + 
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); susc_plot

trans_plot <- prior_values %>% 
  ggplot() +
  geom_density(aes(transmissibility), col = '#D90368', fill = '#D90368',
               alpha = 0.4, lwd = 0.8) +
  # geom_vline(xintercept = min_trans, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = max_trans, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = 'transmissibility') + xlim(c(0,1)) +
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); trans_plot

init_inf_plot <- prior_values %>% 
  ggplot() +
  geom_ribbon(aes(x = init_inf_values, ymax = init_inf_distr, ymin = 0), fill = '#D17B0F',
              alpha = 0.4) +
  geom_line(aes(x = init_inf_values, y = init_inf_distr), col = '#D17B0F',
            lwd = 0.8) +
  # geom_vline(xintercept = min_log_init_inf, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = max_log_init_inf, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = "raw initial infected parameter") + 
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); init_inf_plot

init_inf_plot_10 <- prior_values %>% 
  ggplot() +
  geom_density(aes(10^init_inf_samples/1000), fill = '#D17B0F', col = '#D17B0F',
               alpha = 0.4, lwd = 0.8) +
  # geom_vline(xintercept = 10^min_log_init_inf/1000, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = 10^max_log_init_inf/1000, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = "initial infected (1000s)") + 
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); init_inf_plot_10

imd_spline_plot <- prior_values %>% 
  ggplot() +
  geom_ribbon(aes(x = imd_spline_values, ymax = imd_spline_distr, ymin = 0), fill = '#3C153B',
              alpha = 0.4) +
  geom_line(aes(x = imd_spline_values, y = imd_spline_distr), col = '#3C153B',
            lwd = 0.8) +
  # geom_vline(xintercept = min_spline, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = max_spline, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = "raw IMD spline parameter") + 
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); imd_spline_plot

imd_1_plot <- prior_values %>% 
  ggplot() +
  geom_ribbon(aes(x = exp(imd_spline_values), ymax = imd_spline_distr, ymin = 0), fill = '#3C153B',
              alpha = 0.4) +
  geom_line(aes(x = exp(imd_spline_values), y = imd_spline_distr), col = '#3C153B',
            lwd = 0.8) +
  # geom_vline(xintercept = min_spline, lty = 2, alpha = 0.5, lwd = 0.8) + 
  # geom_vline(xintercept = max_spline, lty = 2, alpha = 0.5, lwd = 0.8) + 
  labs(y = '', x = "relative IMD 1 IHR\n(compared to IMD 3)") + 
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); imd_1_plot

primary_plot <- prior_values %>% 
  filter(reporting_values <= 0.2) %>% 
  ggplot() +
  geom_ribbon(aes(x = reporting_values, ymax = prim_rates, ymin = 0), fill = '#8B1E3F',
              alpha = 0.4) +
  geom_line(aes(x = reporting_values, y = prim_rates), col = '#8B1E3F',
            lwd = 0.8) +
  labs(y = '', x = 'infection-gp ratio') + 
  theme_lg() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); primary_plot

secondary_plot <- prior_values %>% 
  filter(reporting_values <= 0.05) %>% 
  ggplot() +
  geom_ribbon(aes(x = reporting_values, ymax = sec_rates, ymin = 0), fill = '#DB4C40',
              alpha = 0.4) +
  geom_line(aes(x = reporting_values, y = sec_rates), col = '#DB4C40',
            lwd = 0.8) +
  labs(y = '', x = 'infection-hospitalisation ratio') + 
  theme_lg() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()); secondary_plot

reff_plot + susc_plot + trans_plot + init_inf_plot + init_inf_plot_10 +
  primary_plot + secondary_plot + imd_spline_plot + imd_1_plot
ggsave(gsub("epids.png", "prior_pars.png", .args[2]), width = 12, height = 10, dpi = 600)
 

### subtype -Inf binomial example

if(F){
  weekly_props %>% 
    ggplot() +
    geom_line(aes(date_formatted, modelled_proportion), lwd = 0.8) + 
    geom_label(x=as.Date('25-12-2023', format = '%d-%m-%Y'),y=0.8, 
               label = 'modelled proportion') +  
    geom_label(x=as.Date('25-12-2023', format = '%d-%m-%Y'),y=0.55, 
               label = 'ukhsa test\nproportion', col = 'red') + 
    geom_line(aes(date_formatted, proportion), col = 'red', lwd = 0.8, alpha = 0.5) + 
    geom_point(aes(date_formatted, proportion, size = total_flu), col = 'red') + 
    theme_lg() + 
    labs(x = '', y = 'proportion H1N1', size = 'total positive\nflu tests')
}








