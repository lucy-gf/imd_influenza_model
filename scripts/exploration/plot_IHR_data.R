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
  file.path("data", "ukhsa", "parameter_batch_AH1N1.rds"),
  file.path("output", "figures", "exploration", "ukhsa", "ukhsa_IHR_strain.png")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','age_grp_assignment.R'))

## read in UKHSA data
ah1n1 <- readRDS(.args[1])
ah3n2 <- readRDS(gsub("H1N1", "H3N2", .args[1]))
b <- readRDS(gsub("AH1N1", "B", .args[1]))

priors <- readRDS(gsub("parameter_batch_AH1N1", "hospital_rate_prior", .args[1]))

## plot priors

prior_samples <- priors$samples %>% 
  mutate(age_group = case_when(
    name %like% "-0" ~ '0-4',
    name %like% "-5" ~ '5-14',
    name %like% "-15" ~ '15-44',
    name %like% "-45" ~ '45-64',
    name %like% "-65" ~ '65+'
  ),
  interval = case_when(
    name %like% "-0" ~ '[0, 5)',
    name %like% "-5" ~ '[5, 15)',
    name %like% "-15" ~ '[15, 45)',
    name %like% "-45" ~ '[45, 65)',
    name %like% "-65" ~ '[65, Inf)'
  ),
  risk_group = case_when(
    name %like% 'low_risk' ~ 'low_risk',
    name %like% 'high_risk' ~ 'high_risk'
  )) %>% 
  left_join(priors$log_normal, by = c('interval','risk_group'))

prior_samples$age_group <- factor(prior_samples$age_group,
                                  levels = c('0-4','5-14','15-44','45-64','65+'))

prior_samples %>% 
  ggplot() + 
  geom_violin(aes(age_group, value, fill = risk_group, color = risk_group,
                  group = interaction(risk_group, age_group)),
              alpha = 0.4, lwd = 0.8) +
  geom_point(aes(x = age_group, exp(meanlog), color = risk_group),
             position = position_dodge(width = 0.9)) + 
  labs(x = 'Age group', y = 'Infection-Hospitalisation Ratio (prior)', fill = '', color = '') + 
  scale_fill_manual(values = risk_colors) +
  scale_color_manual(values = risk_colors) +
  theme_bw()

## plot posteriors

posteriors <- rbind(ah1n1, ah3n2, b) %>% 
  filter(!grepl('flat',p_lbl)) %>% 
  mutate(val = as.numeric(gsub('hosp_','',par))) %>% 
  mutate(age_num = case_when(val < 15 ~ val-9, T~ val-14), 
         risk_group = case_when(val < 15 ~ 'low_risk', T~'high_risk'),
         age_group = case_when(
           age_num == 1 ~ '0-4',
           age_num == 2 ~ '5-14',
           age_num == 3 ~ '15-44',
           age_num == 4 ~ '45-64',
           age_num == 5 ~ '65+'
         )) 

posteriors$age_group <- factor(posteriors$age_group,
                               levels = c('0-4','5-14','15-44','45-64','65+'))

posteriors %>% 
  ggplot() + 
  geom_density(aes(x = value, fill = risk_group, col = risk_group, 
                   group = interaction(season, risk_group)), 
               alpha = 0.5) + 
  scale_fill_manual(values = risk_colors) +
  scale_color_manual(values = risk_colors) +
  theme_bw() + facet_grid(age_group~subtype, scales = 'free')

posteriors %>% 
  ggplot() + 
  geom_violin(aes(age_group, value, fill = risk_group, color = risk_group,
                  group = interaction(risk_group, age_group)),
              alpha = 0.4, lwd = 0.8) +
  facet_grid(season ~ subtype, scales = 'free') +
  labs(x = 'Age group', y = 'Infection-Hospitalisation Ratio (posterior)', fill = '', color = '') + 
  scale_fill_manual(values = risk_colors) +
  scale_color_manual(values = risk_colors) +
  theme_bw()

## by subtype

plot_pri_post <- function(strain){
  
  posteriors %>% 
    filter(subtype == strain) %>% 
    rename(posterior = value) %>% 
    left_join(prior_samples %>% select(risk_group, age_group, np, value, meanlog) %>% rename(prior = value), 
              by = c('risk_group','age_group','np')) %>% 
    group_by(age_group, risk_group, season) %>% 
    mutate(mean_posterior = mean(posterior)) %>% 
    ggplot() + 
    geom_density(aes(x = posterior, fill = risk_group, color = risk_group,
                     group = interaction(risk_group, age_group)),
                 alpha = 0.4, lwd = 0.8) +
    geom_density(aes(x = prior, color = risk_group,
                     group = interaction(risk_group, age_group)),
                 alpha = 0, lwd = 0.8, lty = 2) +
    geom_vline(aes(xintercept = exp(meanlog), color = risk_group), lty = 2) +
    geom_vline(aes(xintercept = mean_posterior, color = risk_group)) +
    facet_grid(age_group ~ season, scales = 'free') +
    scale_x_log10() +
    labs(x = 'Infection-Hospitalisation Ratio (log 10 scale)', fill = '', color = '') + 
    scale_fill_manual(values = risk_colors) +
    scale_color_manual(values = risk_colors) +
    theme_bw()
  
}

plot_pri_post('AH1N1')
plot_pri_post('AH3N2')
plot_pri_post('B')







