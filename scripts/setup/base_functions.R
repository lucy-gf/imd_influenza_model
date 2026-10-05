#### BASE FUNCTIONS FOR ANALYSIS ####

na_to_0 <- function(x){
  x[is.na(x)] <- 0
  x
}

c_to_0 <- function(v){if(length(v[v=='c']) > 0){v[v=='c'] <- 0; v}else{v}}

susc_vector <- function(x,
                        rep1 = 3,
                        rep2 = 4,
                        rep3 = 2){
  
  x <- unname(unlist(x))
  
  if(length(x) %notin% 2:3){stop('Length not 2 or 3')}
  if(length(x) == 2){x <- c(x[1], 1, x[2])}
  
  c(rep(x[1],rep1), rep(x[2],rep2), rep(x[3],rep3))
  
}

imd_spline <- function(imd_pars){
  
  if(length(imd_pars) != 2){stop('Input length not 2')}
  
  # IMD 1 is exp(imd_pars[1]) times the reporting rate of IMD 3,
  # IMD 5 is exp(imd_pars[1]) times the reporting rate of IMD 3
  
  x <- exp(imd_pars)
  
  c(x[1], 1 + (x[1] - 1)/2, 1, 1 + (x[2] - 1)/2, x[2])
  
}

## GGPLOT THEME ##

theme_lg <- function(gridline_x = FALSE, gridline_y = TRUE) {
  
  gridline <- element_line(
    linewidth = 0.15,
    color = "#999999"
  )
  
  gridline_x <- if (isTRUE(gridline_x)) {
    gridline
  } else {
    element_blank()
  }
  
  gridline_y <- if (isTRUE(gridline_y)) {
    gridline
  } else {
    element_blank()
  }
  
  theme_bw(
  ) +
    theme(
      text = element_text(family = "Helvetica"),
      plot.title = element_text(
        size = 18,
        face = "bold",
        color = "black",
        margin = margin(b = 10, l = 15)
      ),
      plot.subtitle = element_text(
        size = 14,
        color = "#333333",
        margin = margin(b = 10)
      ),
      plot.caption = element_text(
        size = 13,
        color = "#777777",
        margin = margin(t = 15),
        hjust = 0
      ),
      axis.text = element_text(
        size = 12,
        color = "#1F1F1F"
      ),
      axis.title = element_text(
        size = 14,
        color = "#1F1F1F"
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = gridline_x,
      panel.grid.major.y = gridline_y,
      axis.ticks.x = element_line(
        linetype = "solid",
        linewidth = 0.25,
        color = "#333333"
      ),
      axis.ticks.y = element_blank(),
      axis.ticks.length.x = unit(4, units = "pt"),
      strip.background = element_rect(
        fill = "#F0F0F0",
        color = "#777777"
      ),
      strip.text = element_text(
        size = 12,
        # face = "semibold",
        color = "black"
      )
    ) 
  
}

# library(datasauRus)
# 
# datasaurus_dozen %>%
#   filter(dataset %in% c('dino', 'star', 'bullseye', 'slant_down')) %>%
#   ggplot() +
#   geom_point(aes(x, y), size = 2) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(c(0.00075, 0.025))) +
#   theme_lg() + ggtitle('Datasaurus') + labs(x = '', y = '') +
#   facet_wrap(dataset ~ .)

## FUNCTION TO CONVERT BETWEEN SUBTYPE SYNTAX ##

convert_to_subtype <- function(string_vec){
  
  for(k in 1:length(string_vec)){
    
    text <- string_vec[k]
    
    updated_text <- gsub('flu|pdm09|_', '', text)
    
    string_vec[k] <- toupper(updated_text)
    
  }
  
  string_vec
  
}


## FUNCTION TO ASSIGN AGE GROUPS TO BROAD AGE GROUPS ##

# PURPOSE
# 1. risk: Risk groups (if not using UKHSA data)
# 2. vacc_cov: Vaccination coverage
# 3. vacc_eff: Vaccine efficacy
# 4. susc: Population susceptibility
# 5. reporting rates: Vaccine efficacy

fcn_assign_ages <- function(
    child_value,
    adult_value,
    older_adult_value,
    age_labels_in,
    purpose = NULL
    ){
  
  n_c <- 3
  n_a <- 3
  n_o_a <- length(age_labels_in) - n_c - n_a
  
  vector <- c(rep(child_value, n_c),
              rep(adult_value, n_a),
              rep(older_adult_value, n_o_a))

  return(vector)
  
  }


## function to weight participants

weight_participants <- function(part,
                                eth_age_sex_structure = eth_age_sex,
                                weighting = NULL,
                                group_vars = NULL,
                                truncation_percentile = c(0.05, 0.95) # c(0,1)
) {
  # align p_gender case
  part <- part %>% mutate(p_gender = tolower(p_gender))
  eth_age_sex_structure <- eth_age_sex_structure %>%
    mutate(p_gender = tolower(p_gender)) %>%
    select(!value)
  
  # add dataframe for day_week
  if ("day_week" %in% weighting) {
    day_week_structure <- data.table(
      day_week = c("weekday", "weekend"),
      proportion = c(5 / 7, 2 / 7)
    )
    if ("c_day_of_week" %in% group_vars) {
      stop("Can't weight by day_week and group by c_day_of_week")
    }
  }
  
  # specify which of age/gender/ethnicity are being weighted for
  age_gender_vec <- weighting[weighting %like% "age|gender|ethn"]
  
  # age_gender_vec plus p_adult_child if this is in the group_vars vector and p_age_group is in weighting vec
  a_g_e_grouping_vars <- if ("p_adult_child" %in% group_vars & "p_age_group" %in% weighting) {
    c(age_gender_vec, "p_adult_child")
  } else {
    age_gender_vec
  }
  
  a_g_e_df <- if (length(age_gender_vec) > 0) {
    eth_age_sex_structure %>%
      group_by(!!!syms(a_g_e_grouping_vars)) %>%
      summarise(proportion = sum(proportion))
  } else {
    data.frame()
  }
  
  # weight by age and/or gender within specified groups
  part_age <- part %>%
    group_by(!!!syms(unique(c(group_vars, age_gender_vec)))) %>%
    summarise(sample_totals = n()) %>%
    ungroup() %>%
    group_by(!!!syms(group_vars)) %>%
    mutate(group_total = sum(sample_totals))
  
  if (nrow(a_g_e_df) > 0) {
    part_age <- part_age %>% left_join(a_g_e_df, by = (a_g_e_grouping_vars))
    
    not_in_prop_groups <- part_age %>%
      group_by(!!!syms(group_vars)) %>%
      summarise(sum_in = sum(proportion, na.rm = T)) %>%
      ungroup()
    
    true_props <- a_g_e_df %>%
      group_by(!!!syms(group_vars)) %>%
      summarise(true_prop = sum(proportion, na.rm = T)) %>%
      ungroup()
    
    not_in_prop_groups <- not_in_prop_groups %>% 
      left_join(true_props, by = group_vars) %>% 
      mutate(proportion_of_group = sum_in/true_prop) %>% 
      select(!c(sum_in, true_prop))
    
    for (i in 1:nrow(part_age)) {
      if (is.na(part_age$proportion[i]) & "p_gender" %in% age_gender_vec) {
        if ((is.na(part_age$p_gender[i]) | part_age$p_gender[i] == "other")) {
          if ("p_ethnicity" %in% age_gender_vec) {
            if (part_age$p_ethnicity[i] %like% "Prefer") {
              if ("p_age_group" %in% a_g_e_grouping_vars) {
                part_age[i, ] <- suppressMessages(part_age[i, ] %>% rows_update(
                  a_g_e_df %>% group_by(p_age_group) %>%
                    summarise(proportion = sum(proportion)) %>%
                    left_join(part_age %>% group_by(!!!syms(a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "ethn|gender"])) %>%
                                summarise(sample_totals = sum(sample_totals))),
                  by = a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "ethn|gender"],
                  unmatched = "ignore"
                ))
              }
            }
          } else {
            part_age[i, ] <- suppressMessages(part_age[i, ] %>% rows_update(
              a_g_e_df %>% group_by(!!!syms(a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "gender"])) %>%
                summarise(proportion = sum(proportion)) %>%
                left_join(part_age %>% group_by(!!!syms(a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "gender"])) %>%
                            summarise(sample_totals = sum(sample_totals))),
              by = a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "gender"],
              unmatched = "ignore"
            ))
          }
        }
      } else {
        if (is.na(part_age$proportion[i]) & "p_ethnicity" %in% age_gender_vec) {
          if (part_age$p_ethnicity[i] %like% "Prefer" & length(a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "ethn"]) > 0) {
            part_age[i, ] <- suppressMessages(part_age[i, ] %>% rows_update(
              a_g_e_df %>% group_by(!!!syms(a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "ethn"])) %>%
                summarise(proportion = sum(proportion)) %>%
                left_join(part_age %>% group_by(!!!syms(a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "ethn"])) %>%
                            summarise(sample_totals = sum(sample_totals))),
              by = a_g_e_grouping_vars[!a_g_e_grouping_vars %like% "ethn"],
              unmatched = "ignore"
            ))
          }
        }
      }
    }
    
    # scale up proportions to account for missing categories
    if (length(group_vars) > 0) {
      part_age <- part_age %>%
        left_join(not_in_prop_groups, by = group_vars) %>%
        mutate(proportion = proportion / proportion_of_group) %>%
        select(!proportion_of_group)
    } else {
      part_age <- part_age %>%
        mutate(proportion_of_group = not_in_prop_groups$proportion_of_group[1]) %>%
        mutate(proportion = proportion / proportion_of_group) %>%
        select(!proportion_of_group)
    }
  } else {
    part_age <- part_age %>% mutate(proportion = NA)
  }
  
  part_age <- part_age %>% mutate(
    sample_prop = sample_totals / nrow(part),
    post_strat_weight = case_when(
      is.na(proportion) ~ 1,
      T ~ proportion / sample_prop
    )
  )
  
  weighted_data <- part %>%
    select(p_id, !!!(group_vars), !!!syms(weighting))
  if (length(age_gender_vec) == 0) {
    weighted_data <- weighted_data %>% mutate(post_strat_weight = 1)
  } else {
    weighted_data <- weighted_data %>%
      left_join(unique(part_age %>% select(!!!(group_vars), !!!syms(age_gender_vec), post_strat_weight)), by = unique(c(group_vars, age_gender_vec)))
  }
  
  if ("day_week" %in% weighting) {
    weighted_data <- weighted_data %>%
      left_join(
        day_week_structure %>% left_join(
          part %>% group_by(day_week, !!!syms(group_vars)) %>%
            summarise(sample_tot = n()) %>% ungroup() %>%
            group_by(!!!syms(group_vars)) %>% mutate(pop = sum(sample_tot)) %>%
            mutate(sample_proportion = sample_tot / pop) %>% ungroup(),
          by = "day_week"
        ) %>%
          mutate(weight_mult = proportion / sample_proportion) %>% select(day_week, !!!syms(group_vars), weight_mult),
        by = c("day_week", group_vars)
      ) %>%
      mutate(post_strat_weight = case_when(
        is.na(day_week) ~ post_strat_weight,
        T ~ post_strat_weight * weight_mult
      )) %>%
      select(!weight_mult)
  }
  
  # truncate the weights that fall above the 99th percentile or below the 1% percentile
  percs <- quantile(weighted_data$post_strat_weight, truncation_percentile)
  if (n_distinct(weighted_data$post_strat_weight) >= 20) {
    weighted_data <- weighted_data %>%
      mutate(post_strat_weight = case_when(
        post_strat_weight < percs[1] ~ percs[1],
        post_strat_weight > percs[2] ~ percs[2],
        T ~ post_strat_weight
      ))
  }
  
  weighted_data
}

neg_bin_fcn <- function(vec){
  
  # calc. mean and variance, used for initial values
  m <- mean(vec)
  v <- var(vec)
  # cat('Mean = ', m, ', Var = ', v, ' Vec = ', vec, ' -- ', sep = '')
  
  if(length(vec) > 1){
    if(sum(vec != 0)){
      outs = optim(c(mu = m, k = (v - m)/m^2), lower = c(mu = 1e-5, k = 1e-5), nb_loglik, x = vec, method = "L-BFGS-B")
      ret <- as.numeric(outs$par)
    } else{
      ret <- c(0,0)
    }
  }else{
    if(length(vec) == 1){ret <- c(vec, 0)}
  }
  
  # paste(ret, collapse = '_')
  ret[1]
  
}

# Log-likelihood function for negative binomial
nb_loglik <- function(x, par) {
  k <- par[["k"]]
  mean <- par[["mu"]]
  ll <- rep(NA_real_, length(x))
  ll <- dnbinom(x, mu = mean, size = 1/k, log = TRUE)
  return(-sum(ll))
}

## function to calculate matrix from survey data

calc_matrix <- function(data,
                        reciprocal = F,
                        age_pop,
                        age_limits_in = seq(0, 80, 5)){
  
  data_grouped <- assign_age_groups(data, estimated_contact_age = RC_dist,
                                    age_limits = age_limits_in)
  
  matrix <- compute_matrix(data_grouped)
  
  if(reciprocal){
    
    for(k in 1:nrow(age_pop)){
      
      if(age_pop$age_group[k] != colnames(matrix$matrix)[k]){
        stop(paste0('Age group names not aligned (', k, ')'))
      }
      
    }
    
    m_ij_p_i <- age_pop$population * matrix$matrix
    m_ji_p_j <- t(m_ij_p_i)
    denominator <- 2*age_pop$population
    
    reciprocal_matrix <- (m_ij_p_i + m_ji_p_j)/denominator
    
    dt <- reciprocal_matrix 
    
  }else{
    
    dt <- matrix$matrix
    
  }
  
  gg_matrix <- data.table(dt) %>% 
    mutate(part_age = colnames(data.table(matrix$matrix))) %>% 
    pivot_longer(!part_age) %>% 
    rename(cont_age = name) 
  
  gg_matrix$part_age <- factor(gg_matrix$part_age,
                               levels = unique(gg_matrix$part_age))
  gg_matrix$cont_age <- factor(gg_matrix$cont_age,
                               levels = unique(gg_matrix$cont_age))
  
  gg_matrix 
  
}

## function to plot from survey data

ggplot_matrix <- function(data,
                          age_limits_in = seq(0, 80, 5)){
  
  calc_matrix(data,
              age_limits_in = age_limits_in) %>% 
    ggplot() + 
    geom_tile(aes(x = part_age, y = cont_age, 
                  fill = value)) + 
    labs(x = 'Participant age group', y = 'Contact age group', fill = 'Mean contacts') +
    theme_bw() + scale_fill_viridis(limits = c(0,NA)) +
    coord_fixed()
  
}

# function to plot from survey data

ggplot_reciprocal_matrix <- function(data,
                                     age_pop = age_population,
                                     age_limits_in = seq(0, 80, 5)){
  
  calc_matrix(data, reciprocal = T,
              age_limits_in = age_limits_in) %>% 
    ggplot() + 
    geom_tile(aes(x = part_age, y = cont_age, 
                  fill = value)) + 
    labs(x = 'Participant age group', y = 'Contact age group', fill = 'Mean contacts') +
    theme_bw() + scale_fill_viridis(limits = c(0,NA)) +
    coord_fixed()
  
}




