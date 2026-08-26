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
source(file.path('scripts','setup','age_grp_assignment.R'))

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





