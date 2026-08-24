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
  file.path("output", "figures", "exploration", "ukhsa_weekly_strain.png")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','colors.R'))
source(file.path('scripts','setup','age_grp_assignment.R'))

## read in UKHSA data
ukhsa_dat <- read_ods(.args[1], sheet = 48, skip = 3)

# rename columns
colnames(ukhsa_dat) <- c('date', 'week', 'flu_a_unsubtyped', 'flu_a_h1n1pdm09',
                         'flu_a_h3n2', 'flu_b', 'positivity')

ukhsa_dat <- ukhsa_dat %>% mutate(positivity = positivity/100,
                                  date_formatted = as.Date(date, format = '%d %b %Y'),
                                  flu_a = flu_a_unsubtyped + flu_a_h1n1pdm09 + flu_a_h3n2,
                                  flu_tot = flu_a + flu_b)

p_positivity <- ukhsa_dat %>% 
  ggplot() + 
  geom_line(aes(date_formatted, positivity)) +
  theme_bw() + labs(x = '', y = 'Overall test positivity')

p_tests <- ukhsa_dat %>% 
  mutate(tests = flu_tot/positivity) %>% 
  mutate(cumulative_na = cumsum(positivity == 0)) %>% 
  filter(positivity > 0, flu_tot > 3) %>% 
  ggplot() + 
  geom_line(aes(date_formatted, tests, col = positivity, group = cumulative_na)) +
  scale_y_continuous(breaks = 2000*0:10, limits = c(0,NA)) +
  scale_color_viridis() + 
  theme_bw() + labs(x = '', y = 'Overall tests carried out', col = 'Positivity')

stacked_plot <- function(min_year = 2021,
                         all_a = F){
  if(!all_a){
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, flu_a, flu_b) %>% 
      pivot_longer(!date_formatted) %>% 
      ggplot() + 
      geom_bar(aes(date_formatted, value, fill = name),
               stat = 'identity', position = 'stack', width = 7) +
      scale_fill_manual(values = flu_strain_colors,
                        labels = flu_strain_names) +
      labs(x = '', y = 'Positive tests', fill = '') +
      theme_bw()
  }else{
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, starts_with('flu_')) %>% 
      select(!c(flu_a, flu_tot)) %>% 
      pivot_longer(!date_formatted) %>% 
      ggplot() + 
      geom_bar(aes(date_formatted, value, fill = name),
               stat = 'identity', position = 'stack', width = 7) +
      scale_fill_manual(values = flu_subtype_colors,
                        labels = flu_subtype_names) +
      labs(x = '', y = 'Positive tests', fill = '') +
      theme_bw()
  }
  
  p
  
}

absolute_plot <- function(min_year = 2021,
                          all_a = F){
  
  if(!all_a){
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, flu_a, flu_b) %>% 
      pivot_longer(!date_formatted) %>% 
      ggplot() + 
      geom_line(aes(date_formatted, value, col = name, group = name),
                lwd = 0.8) +
      scale_color_manual(values = flu_strain_colors,
                         labels = flu_strain_names) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(x = '', y = 'Positive tests', col = '') +
      theme_bw()
  }else{
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, starts_with('flu_')) %>% 
      select(!c(flu_a, flu_tot)) %>% 
      pivot_longer(!date_formatted) %>% 
      ggplot() + 
      geom_line(aes(date_formatted, value, col = name, group = name),
                lwd = 0.8) +
      scale_color_manual(values = flu_subtype_colors,
                         labels = flu_subtype_names) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(x = '', y = 'Positive tests', col = '') +
      theme_bw()
  }
  
  p
  
}

bw_plot <- function(min_year = 2021,
                    min_tests = 50){
  
  p_tests_bw <- ukhsa_dat %>% 
    mutate(cumulative_na = cumsum(flu_tot < min_tests)) %>% 
    filter(positivity > 0, flu_tot >= min_tests, year(date_formatted) >= min_year) %>% 
    ggplot() + 
    geom_line(aes(date_formatted, flu_tot, group = cumulative_na)) +
    scale_y_continuous(limits = c(0,NA)) +
    theme_bw() + labs(x = '', y = 'Pos. tests')
  
  p_tests_bw
  
}

proportion_plot <- function(min_year = 2021,
                            min_tests = 50, 
                            all_a = F){
 
  if(!all_a){
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year,
             flu_tot >= min_tests) %>% 
      select(date_formatted, flu_tot, flu_a, flu_b) %>% 
      pivot_longer(! c(date_formatted, flu_tot)) %>% 
      ggplot() + 
      geom_bar(aes(date_formatted, value/flu_tot, fill = name, group = name),
               stat = 'identity', position = 'stack', width = 7) +
      scale_fill_manual(values = flu_strain_colors,
                        labels = flu_strain_names) +
      labs(x = '', y = 'Proportion of positive tests', fill = '',
           title = paste0('Proportion of positive tests by strain, in weeks with at least ', min_tests, ' tests')) +
      theme_bw()
  }else{
    p <- ukhsa_dat %>% 
      filter(year(date_formatted) >= min_year,
             flu_tot >= min_tests) %>% 
      select(date_formatted, starts_with('flu_')) %>% 
      select(!c(flu_a)) %>% 
      pivot_longer(! c(date_formatted, flu_tot)) %>% 
      ggplot() + 
      geom_bar(aes(date_formatted, value/flu_tot, fill = name, group = name),
               stat = 'identity', position = 'stack', width = 7) +
      scale_fill_manual(values = flu_subtype_colors,
                        labels = flu_subtype_names) +
      labs(x = '', y = 'Proportion of positive tests', fill = '',
           title = paste0('Proportion of positive tests by subtype, in weeks with at least ', min_tests, ' tests')) +
      theme_bw()
  }
  
  p 
  
}

### make patchwork plots

height_val <- 8; width_val <- 9

p_positivity + p_tests + plot_layout(nrow = 2)
ggsave(gsub('strain.png', 'positivity.png', .args[2]), width = width_val, height = height_val)

stacked_plot(min_year = 2021) + stacked_plot(min_year = 2021, all_a = T) +
  plot_layout(nrow = 2)
ggsave(gsub('.png', '_stacked.png', .args[2]), width = width_val, height = height_val)

absolute_plot(min_year = 2021) + absolute_plot(min_year = 2021, all_a = T) +
  plot_layout(nrow = 2)
ggsave(.args[2], width = width_val, height = height_val)

bw_plot() + proportion_plot() + proportion_plot(all_a = T) +
  plot_layout(nrow = 3, heights = c(1,9,9))
ggsave(gsub('.png', '_proportion.png', .args[2]), width = width_val + 3, height = height_val + 2)

bw_plot() + theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) + 
  proportion_plot() + theme(title=element_blank()) + plot_layout(nrow = 2, heights = c(1,9))
ggsave(gsub('.png', '_proportion_A_B.png', .args[2]), width = width_val + 1, height = height_val - 2)

bw_plot() + theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) + 
  proportion_plot(all_a = T) + theme(title=element_blank()) + plot_layout(nrow = 2, heights = c(1,9))
ggsave(gsub('.png', '_proportion_A_ALL_B.png', .args[2]), width = width_val + 2, height = height_val - 2)

#### AGE GROUPED #############

## read in UKHSA data by age group
ukhsa_dat_age_2526 <- read_ods(.args[1], sheet = 23, skip = 3)
ukhsa_dat_age_2425 <- read_ods(gsub('25_2026','24_2025',.args[1]), sheet = 22, skip = 3)
ukhsa_dat_age_2324 <- read_ods(gsub('25_2026','23_2024',.args[1]), sheet = 28, skip = 3)

# rename columns
colnames(ukhsa_dat_age_2526) <- c('date', 'week', 'age_group', 'flu_a_unsubtyped', 'flu_a_h1n1pdm09',
                             'flu_a_h3n2', 'flu_b')
colnames(ukhsa_dat_age_2425) <- c('date', 'age_group', 'flu_a', 'flu_b')
colnames(ukhsa_dat_age_2324) <- c('age_group', 'flu_a_h1n1pdm09', 'flu_a_h3n2', 'flu_a_unsubtyped', 'flu_b')

ukhsa_dat_age_2526 <- ukhsa_dat_age_2526 %>% 
  mutate(age_group = gsub('Less than ', '0-',
                          gsub(' and above| years and over| years-over', '+',
                               gsub('Between | years of age| year of age','',
                                    gsub(' and ', '-', age_group)))))
ukhsa_dat_age_2425 <- ukhsa_dat_age_2425 %>% 
  mutate(age_group = gsub('Less than ', '0-',
                          gsub(' and above| years and over| years-over', '+',
                               gsub('Between | years of age| year of age','',
                                    gsub(' and ', '-', age_group)))))
ukhsa_dat_age_2324 <- ukhsa_dat_age_2324 %>% 
  mutate(age_group = gsub('Less than ', '0-',
                          gsub(' and above|-above| years and over| years-over', '+',
                               gsub('Between | years of age| year of age','',
                                    gsub(' and ', '-', age_group)))))

ukhsa_dat_age_2526$age_group <- factor(ukhsa_dat_age_2526$age_group, 
                                       levels = unique(ukhsa_dat_age_2526$age_group))
ukhsa_dat_age_2425$age_group <- factor(ukhsa_dat_age_2425$age_group, 
                                       levels = unique(ukhsa_dat_age_2425$age_group))
ukhsa_dat_age_2324$age_group <- factor(ukhsa_dat_age_2324$age_group, 
                                       levels = unique(ukhsa_dat_age_2324$age_group))

ukhsa_dat_age_2526 <- ukhsa_dat_age_2526 %>% 
  mutate(date_formatted = as.Date(date, format = '%d %b %Y'),
         flu_a = flu_a_unsubtyped + flu_a_h1n1pdm09 + flu_a_h3n2,
         flu_tot = flu_a + flu_b)

ukhsa_dat_age_2425 <- ukhsa_dat_age_2425 %>% 
  mutate(date_formatted = as.Date(date, format = '%d %b %Y'),
         flu_tot = flu_a + flu_b)
  
proportion_plot_age <- function(min_year = 2021,
                                all_a = F){
  
  if(!all_a){
    ukhsa_dat_age_2526 %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, age_group, flu_tot, flu_a, flu_b) %>% 
      pivot_longer(! c(date_formatted, age_group, flu_tot)) %>% 
      rbind(ukhsa_dat_age_2425 %>% 
              filter(year(date_formatted) >= min_year) %>% 
              select(date_formatted, age_group, flu_tot, flu_a, flu_b) %>% 
              pivot_longer(! c(date_formatted, age_group, flu_tot))) %>% 
      filter(flu_tot > 0) %>% 
      ggplot() + 
      geom_bar(aes(date_formatted, value/flu_tot, fill = name, group = name),
               stat = 'identity', position = 'stack', width = 7) +
      scale_fill_manual(values = flu_strain_colors,
                        labels = flu_strain_names) +
      facet_wrap(age_group ~ .) + 
      labs(x = '', y = 'Proportion of positive tests', fill = '',
           title = paste0('Proportion of positive tests by strain (2024-2026)')) +
      theme_bw()
  }else{
    ukhsa_dat_age_2526 %>% 
      filter(year(date_formatted) >= min_year) %>% 
      select(date_formatted, age_group, starts_with('flu_')) %>% 
      select(!c(flu_a)) %>% 
      pivot_longer(! c(date_formatted, age_group, flu_tot)) %>% 
      filter(flu_tot > 0) %>% 
      ggplot() + 
      geom_bar(aes(date_formatted, value/flu_tot, fill = name, group = name),
               stat = 'identity', position = 'stack', width = 7) +
      scale_fill_manual(values = flu_subtype_colors,
                        labels = flu_subtype_names) +
      facet_wrap(age_group ~ .) + 
      labs(x = '', y = 'Proportion of positive tests', fill = '',
           title = paste0('Proportion of positive tests by subtype (2025-2026 only)')) +
      theme_bw()
  }
  
}

proportion_plot_age() + proportion_plot_age(all_a = T) + plot_layout(nrow = 2)
ggsave(gsub('.png', '_age_weekly.png', .args[2]), width = width_val + 3, height = height_val + 2)

### over whole season

p2324 <- ukhsa_dat_age_2324 %>% 
  pivot_longer(!age_group) %>% 
  ggplot() + 
  geom_bar(aes(x = age_group, y = value, fill = name),
           position = 'fill', stat = 'identity') +
  scale_fill_manual(values = flu_subtype_colors,
                    labels = flu_subtype_names) +
  labs(x = '', y = 'Proportion of positive tests', fill = '',
       title = '2023 - 2024') +
  theme_bw()

p2425 <- ukhsa_dat_age_2425 %>% 
  select(!c(date, flu_tot)) %>% 
  pivot_longer(!c(age_group, date_formatted)) %>% 
  group_by(age_group, name) %>% 
  summarise(value = sum(value)) %>% 
  ggplot() + 
  geom_bar(aes(x = age_group, y = value, fill = name),
           position = 'fill', stat = 'identity') +
  scale_fill_manual(values = flu_strain_colors,
                    labels = flu_strain_names) +
  labs(x = '', y = 'Proportion of positive tests', fill = '',
       title = '2024 - 2025') +
  theme_bw()

p2526 <- ukhsa_dat_age_2526 %>% 
  select(!c(week, date, flu_tot, flu_a)) %>% 
  pivot_longer(!c(age_group, date_formatted)) %>% 
  group_by(age_group, name) %>% 
  summarise(value = sum(value)) %>% 
  ggplot() + 
  geom_bar(aes(x = age_group, y = value, fill = name),
           position = 'fill', stat = 'identity') +
  scale_fill_manual(values = flu_subtype_colors,
                    labels = flu_subtype_names) +
  labs(x = '', y = 'Proportion of positive tests', fill = '',
       title = '2025 - 2026') +
  theme_bw()

p2324 + (p2425 + p2526 + plot_layout(nrow = 1)) +
  plot_layout(nrow = 2, guides = 'collect')
ggsave(gsub('.png', '_age.png', .args[2]), width = width_val + 3, height = height_val + 2)





