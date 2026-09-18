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
        panel.background = element_blank()); p3

p1 + p2 + plot_layout(nrow = 2, heights = c(1,2))
ggsave(gsub("epids.png", "_MCMC_surveillance.png", .args[2]), width = 5, height = 3.5)
p1b + p2b + plot_layout(nrow = 2, heights = c(1,2))
ggsave(gsub("epids.png", "_MCMC_surveillance_2.png", .args[2]), width = 4.8, height = 3.5)



