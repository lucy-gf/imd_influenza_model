## holiday mixing changes, by broad age group ##

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(readxl))
suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(viridis))
suppressPackageStartupMessages(require(socialmixr))
suppressPackageStartupMessages(require(contactsurveys))
options(dplyr.summarise.inform = FALSE) 

.args <- if (interactive()) c(
  file.path("data", "inputs", "contact_matrix.rds"),
  file.path('data','inputs','imd_age_pop.rds'),
  file.path("data", "contact_matrix", "reconnect_part_holidays.rds"),
  file.path("data", "inputs", "holiday_contact_matrix.rds")
) else commandArgs(trailingOnly = TRUE)

source(file.path('scripts','setup','base_functions.R'))
source(file.path('scripts','setup','colors.R'))

broad_age_limits <- c(0, 18, 65)
broad_age_vals <- broad_age_limits[broad_age_limits != 0]
broad_age_labels <- c(paste0(c(0, broad_age_vals[1:length(broad_age_vals)-1]), '-', c(broad_age_vals-1)), paste0(broad_age_vals[length(broad_age_vals)], '+'))

## read in main contact matrix
cm <- readRDS(.args[1])

#### POPULATION DATA ####

imd_age_pop <- readRDS(.args[2])
broad_ages <- data.table(
  age_grp = unique(imd_age_pop$age_grp),
  age_group = fcn_assign_ages('[0,18)','[18,65)','[65,Inf)', unique(imd_age_pop$age_grp)))
broad_age_pop <- imd_age_pop %>% 
  left_join(broad_ages, by = 'age_grp') %>%
  group_by(age_group) %>% 
  summarise(population = sum(pop))

#### RECONNECT DATA ####

c_survey <- contactsurveys::download_survey("https://doi.org/10.5281/zenodo.16845074",
                                            overwrite = T)
reconnect <- load_survey(c_survey, participant_key = c("part_id", "cont_id"))

# change age grouping of participants
reconnect$participants <- reconnect$participants %>% 
  mutate(part_age_group = cut(part_age_exact, 
                              breaks = c(-Inf, broad_age_vals, Inf),
                              labels = broad_age_labels,
                              right = F))

# check old issue resolved
if(sum(table(reconnect$participants$part_id) != 1) != 0){warning('Some participant IDs repeated')}

#### ADD WEIGHTS ####
## load population representativeness data
{
  # from census 2021 (https://www.ons.gov.uk/peoplepopulationandcommunity/culturalidentity
  # /ethnicity/datasets/ethnicgroupbyageandsexinenglandandwales)
  eth_age_sex <- suppressWarnings(data.table(read_xlsx(file.path("data","population","ethnicgroupagesex11.xlsx"), sheet = 6, skip = 3)))
  eth_age_sex <- eth_age_sex[get(colnames(eth_age_sex)[1])=='K04000001',]
  eth_age_sex[grepl('100 or over', Age), Age := 100]
  eth_age_sex <- eth_age_sex[, lapply(.SD, c_to_0)]
  eth_age_sex[, 3:ncol(eth_age_sex)] <- lapply(eth_age_sex[, 3:ncol(eth_age_sex)], as.numeric)
  for(ethn in c('Asian','Black','Mixed','White','Other')){
    for(sex in c('Female','Male')){
      vec <- (substr(colnames(eth_age_sex),1,5) == ethn) & grepl(sex,colnames(eth_age_sex))
      eth_age_sex$next_col <- rowSums(eth_age_sex[, ..vec])
      colnames(eth_age_sex)[length(colnames(eth_age_sex))] <- paste0(ethn,'_',sex)
    }
  }
  eth_age_sex <- eth_age_sex %>% select('Age',contains('_')) %>% mutate('p_age_group' = cut(Age,
                                                                                            breaks = c(-Inf, broad_age_vals, Inf),
                                                                                            labels = broad_age_labels,
                                                                                            right = F)) %>% 
    mutate('p_adult_child' = cut(Age, breaks = c(-Inf, 18, Inf), labels = c('Child','Adult'), right = F)) %>% 
    select(!Age) %>% pivot_longer(!c(p_age_group, p_adult_child)) %>% 
    separate_wider_delim(name, delim = "_", names = c("p_ethnicity", "p_gender")) %>% 
    summarise(value = sum(value), .by = c(p_age_group, p_adult_child, p_ethnicity, p_gender)) %>% ungroup() %>% 
    mutate(proportion = value/sum(value)) %>% complete(p_adult_child, p_age_group, p_ethnicity, p_gender,
                                                       fill = list(value = 0, proportion = 0))
}

## join weights
part_data <- reconnect$participants %>% 
  mutate(part_gender = case_when(part_gender == 'F' ~ 'female', part_gender == 'M' ~ 'male'))
colnames(part_data) <- gsub('part_', 'p_', colnames(part_data))
weights <- weight_participants(part_data, 
                               weighting = c('p_age_group','p_gender','p_ethnicity','day_week'),
                               group_vars = 'p_age_group')
if(F){ ## can plot
  ggplot(weights) + 
  geom_hline(yintercept = 1, lty = 2, alpha = 0.5) +
  geom_jitter(aes(p_ethnicity, y = post_strat_weight, col = p_ethnicity, shape = day_week), height = 0.01, alpha = 0.7) + 
  theme_lg() + facet_grid(p_gender~p_age_group, scales = 'fixed')
  }

# merge
reconnect$participants <- reconnect$participants %>% left_join(weights %>% 
                                       select(p_id, post_strat_weight) %>% rename(part_id = p_id),
                                     by = 'part_id')

reconnect <- reconnect %>% weigh("post_strat_weight")

# not important about assigning large group contacts' ages since the
# age groups are the same as the estimated bands used in the survey
RC_dist <- socialmixr::contact_age_distribution(reconnect)

#### ADDITIONAL PARTICIPANT COLUMNS ####

# add country and holiday factor
rc_holidays <- readRDS(.args[3]) %>% rename(part_id = p_id) %>% 
  mutate(IS_XMAS = case_when(p_termtime == 'xmas_holiday' ~ T, T ~ F))

reconnect$participants <- reconnect$participants %>% 
  left_join(rc_holidays, by = 'part_id')

reconnect$participants %>% filter(p_country == 'England') %>% 
  group_by(part_age_group, c_contact_date, p_broad_hols, p_termtime) %>% 
  count() %>% 
  drop_na() %>% 
  ggplot() +
  geom_bar(aes(x = c_contact_date, y = n, fill = p_termtime),
           stat = 'identity', position = 'stack', width = 1) +
  facet_grid(part_age_group ~ ., scales = 'free') +
  theme_lg() + labs(x = '', fill = '')

## filter to England only
reconnect$participants <- reconnect$participants %>% 
  filter(p_country == 'England')
reconnect$contacts <- reconnect$contacts %>% 
  filter(part_id %in% reconnect$participants$part_id)

#### MEAN CONTACTS ####
## by broad age and holiday

hols_analysis <- copy(reconnect$participants) %>% 
  select(part_id, part_age_group, IS_XMAS) %>% 
  left_join(reconnect$contacts %>% 
              group_by(part_id) %>% count(), 
            by = 'part_id') %>% filter(!is.na(IS_XMAS))
hols_analysis$n <- na_to_0(hols_analysis$n)
hols_analysis <- hols_analysis %>% select(!part_id)

mean_contacts_hols <- hols_analysis[, lapply(.SD, neg_bin_fcn), by = c('part_age_group','IS_XMAS')]

mean_contacts_hols %>% 
  ggplot() + geom_point(aes(x = part_age_group, y = n, col = IS_XMAS),
                        size = 3) + 
  theme_lg() +
  ylim(c(0,NA))

#### CONTACT MATRICES ####

# filtered to IS_XMAS
reconnect_XMAS <- copy(reconnect)
reconnect_XMAS$participants <- reconnect_XMAS$participants %>% 
  filter(IS_XMAS)
reconnect_XMAS$contacts <- reconnect_XMAS$contacts %>% 
  filter(part_id %in% reconnect_XMAS$participants$part_id)

# filtered to !IS_XMAS
reconnect_NOT_XMAS <- copy(reconnect)
reconnect_NOT_XMAS$participants <- reconnect_NOT_XMAS$participants %>% 
  filter(!IS_XMAS)
reconnect_NOT_XMAS$contacts <- reconnect_NOT_XMAS$contacts %>% 
  filter(part_id %in% reconnect_NOT_XMAS$participants$part_id)

# calculated
xmas_matrix <- calc_matrix(reconnect_XMAS, 
                           reciprocal = T,
                           age_pop = broad_age_pop,
                           age_limits_in = broad_age_limits)
not_xmas_matrix <- calc_matrix(reconnect_NOT_XMAS, 
                               reciprocal = T,
                               age_pop = broad_age_pop,
                               age_limits_in = broad_age_limits)

# combine
holiday_cm <- not_xmas_matrix %>% 
  left_join(xmas_matrix, by = c('part_age','cont_age'), suffix = c('_nx','_x')) %>% 
  mutate(diff = value_x - value_nx,
         propdiff = diff/value_nx)

#### PLOT ####
max_value <- max(c(xmas_matrix$value,
                   not_xmas_matrix$value))

px <- xmas_matrix %>% 
  ggplot() + 
  geom_tile(aes(x = part_age, y = cont_age, 
                fill = value)) +
  geom_label(aes(x = part_age, y = cont_age,
                label = round(value, 2)), col = 'black', fill = 'white', alpha = 0.5) + 
  labs(x = 'Participant age group', y = 'Contact age group', fill = 'Mean contacts',
       title = 'Winter holidays') +
  theme_lg() + scale_fill_viridis(limits = c(0,max_value)) +
  coord_fixed(); px

pnx <- not_xmas_matrix %>% 
  ggplot() + 
  geom_tile(aes(x = part_age, y = cont_age, 
                fill = value)) +
  geom_label(aes(x = part_age, y = cont_age,
                 label = round(value, 2)), col = 'black', fill = 'white', alpha = 0.5) + 
  labs(x = 'Participant age group', y = 'Contact age group', fill = 'Mean contacts',
       title = 'Not winter holidays') +
  theme_lg() + scale_fill_viridis(limits = c(0,max_value)) +
  coord_fixed(); pnx

pdiff <- holiday_cm %>% 
  ggplot() + 
  geom_tile(aes(x = part_age, y = cont_age, 
                fill = diff)) +
  geom_label(aes(x = part_age, y = cont_age,
                 label = round(value_x - value_nx, 2)), col = 'black', fill = 'white', alpha = 0.5) + 
  labs(x = 'Participant age group', y = 'Contact age group', fill = 'Difference',
       title = 'Difference in mean contacts\n(positive = more in winter holidays)') +
  theme_lg() + scale_fill_gradient2(low = muted("blue"), mid = "white", 
                                    high = muted("red")) +
  coord_fixed(); pdiff

pdiffprop <- holiday_cm %>% 
  ggplot() + 
  geom_tile(aes(x = part_age, y = cont_age, 
                fill = propdiff)) +
  geom_label(aes(x = part_age, y = cont_age,
                 label = paste0(100*round((value_x - value_nx)/value_nx, 2),'%')), 
             col = 'black', fill = 'white', alpha = 0.5) + 
  labs(x = 'Participant age group', y = 'Contact age group', fill = 'Difference',
       title = 'Percentage difference in mean contacts\n(positive = more in winter holidays)') +
  theme_lg() + scale_fill_gradient2(low = muted("green"), mid = "white", 
                                    high = muted("purple")) +
  coord_fixed(); pdiffprop

pnx + px + 
  pdiff + pdiffprop + plot_layout(nrow = 2)
ggsave(file.path('output','figures','exploration','holiday_contacts.png'),
       width = 16, height = 13)

#### SAVE HOLIDAY MATRIX ####
write_rds(holiday_cm, .args[4])
