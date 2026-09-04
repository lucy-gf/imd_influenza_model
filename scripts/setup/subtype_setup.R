## SETTING UP SUBTYPE-YEARS ##

## PLOTTING UKHSA WEEKLY STRAIN DATA ##

#### SETUP ####
suppressMessages(require(ggplot2))
suppressMessages(require(patchwork))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(viridis))
suppressMessages(require(readr))
suppressMessages(require(readODS))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "ukhsa", "annual_influenza_2025_2026.ods"),
  file.path("data", "inputs", "subtype_years.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','base_functions.R'))

## read in UKHSA data
ukhsa_dat <- read_ods(.args[1], sheet = 48, skip = 3)

# rename columns
colnames(ukhsa_dat) <- c('date', 'week', 'flu_a_unsubtyped', 'flu_a_h1n1pdm09',
                         'flu_a_h3n2', 'flu_b', 'positivity')

ukhsa_dat <- ukhsa_dat %>% mutate(positivity = positivity/100,
                                  date_formatted = as.Date(date, format = '%d %b %Y'),
                                  flu_a = flu_a_unsubtyped + flu_a_h1n1pdm09 + flu_a_h3n2,
                                  flu_tot = flu_a + flu_b) %>% 
  mutate(season = case_when(
    week <= 26 ~ paste0(year(date_formatted + 3) - 1, '/', substr(year(date_formatted + 3), 3, 4)),
    T ~ paste0(year(date_formatted + 3), '/', substr(year(date_formatted + 3) + 1, 3, 4))))

absolute_plot <- function(min_year = 2021,
                          all_a = F){
  
  if(!all_a){
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, flu_a, flu_b) %>% 
      pivot_longer(!date_formatted) %>% 
      group_by(date_formatted) %>% mutate(TOTPOS = sum(value)) %>% 
      ggplot() + 
      geom_line(aes(date_formatted, TOTPOS),
                lwd = 0.6, lty = 2) +
      geom_line(aes(date_formatted, value, col = name, group = name),
                lwd = 0.8) +
      scale_color_manual(values = flu_subtype_colors,
                         labels = flu_subtype_names) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(x = '', y = 'Positive tests', col = '') +
      theme_bw(); p
  }else{
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, starts_with('flu_')) %>% 
      select(!c(flu_a, flu_tot)) %>% 
      pivot_longer(!date_formatted) %>% 
      group_by(date_formatted) %>% mutate(TOTPOS = sum(value)) %>% 
      ggplot() + 
      geom_line(aes(date_formatted, TOTPOS),
                lwd = 0.6, lty = 2) +
      geom_line(aes(date_formatted, value, col = name, group = name),
                lwd = 0.8) +
      scale_color_manual(values = flu_subtype_colors,
                         labels = flu_subtype_names) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(x = '', y = 'Positive tests', col = '') +
      theme_bw(); p
  }
  
  p
  
}

abs_timeseries <- absolute_plot(min_year = 2021, all_a = T); abs_timeseries

min_season <- 2023

ukhsa_dat %>% 
  mutate(season = case_when(
    week <= 26 ~ paste0(year(date_formatted) - 1, '/', substr(year(date_formatted), 3, 4)),
    T ~ paste0(year(date_formatted), '/', substr(year(date_formatted) + 1, 3, 4)))) %>% 
  filter(as.numeric(substr(season, 1, 4)) >= min_season) %>% 
  select(season, date_formatted, starts_with('flu_')) %>% 
  select(!c(flu_a, flu_tot)) %>% 
  pivot_longer(!c(season, date_formatted)) %>% 
  rename(subtype = name) %>% 
  group_by(season, subtype) %>% 
  summarise(mean = mean(value)) %>% 
  filter(!subtype %like% 'subtype') %>% 
  ggplot() + geom_bar(aes(x = season, fill = subtype, y = mean),
                      position = 'dodge', stat = 'identity') + 
  theme_bw() +
  scale_fill_manual(values = flu_subtype_colors)
  
ukhsa_dat %>% 
  filter(as.numeric(substr(season, 1, 4)) >= min_season) %>% 
  select(season, date_formatted, starts_with('flu_')) %>% 
  select(!c(flu_a, flu_tot)) %>% 
  pivot_longer(!c(season, date_formatted)) %>% 
  rename(subtype = name) %>% 
  filter(!subtype %like% 'subtype') %>% 
  group_by(season, subtype) %>% 
  mutate(week = 1:n(),
         peak = which.max(value)) %>% 
  filter(week == peak)

ukhsa_dat %>% 
  filter(as.numeric(substr(season, 1, 4)) >= min_season) %>% 
  select(season, date_formatted, starts_with('flu_')) %>% 
  select(!c(flu_a, flu_tot)) %>% 
  pivot_longer(!c(season, date_formatted)) %>% 
  rename(subtype = name) %>% 
  filter(!subtype %like% 'subtype') %>% 
  group_by(season, subtype) %>% 
  mutate(week = 1:n(),
         peak = which.max(value)) %>% 
  ggplot() + 
  geom_hline(aes(yintercept = 51), lty = 3, alpha = 0.5) + 
  geom_vline(aes(xintercept = date_formatted, 
                 col = subtype,
                 alpha = (week == peak)),
             lty = 2) +
  geom_line(aes(date_formatted, value, col = subtype, group = subtype),
            lwd = 0.8) +
  scale_alpha_manual(values = c(0, 1)) + 
  scale_color_manual(values = flu_subtype_colors,
                     labels = flu_subtype_names) +
  scale_y_continuous(expand = expansion(c(0,0.1))) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = '', y = 'Positive tests', col = '') +
  theme_bw()

subtype_info <- ukhsa_dat %>% 
  filter(as.numeric(substr(season, 1, 4)) >= min_season) %>% 
  select(season, date_formatted, starts_with('flu_')) %>% 
  select(!c(flu_a, flu_tot)) %>% 
  pivot_longer(!c(season, date_formatted)) %>% 
  rename(subtype = name) %>% 
  filter(!subtype %like% 'subtype') %>% 
  group_by(season, subtype) %>% 
  mutate(mean = mean(value),
         week = 1:n(),
         peak = which.max(value)) %>% 
  select(season, subtype, peak, mean) %>% unique() %>% ungroup() %>% 
  mutate(peak_qual = case_when(peak > 30 ~ 'Late', T ~ 'Normal')) %>% 
  filter(mean > 50)

write_rds(subtype_info, .args[2])

