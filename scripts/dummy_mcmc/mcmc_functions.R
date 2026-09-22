## MCMC FUNCTIONS ##

TEXT_SAVE_DIR <- file.path('mcmc_output')
if(!dir.exists(TEXT_SAVE_DIR)){dir.create(TEXT_SAVE_DIR)}

## MCMC FITTING, WITH UNKNOWN REPORTING RATES ##
run_mcmc_inference <- function(
    demography_input, 
    vaccinated_input,
    subtype_season,
    cm_input, 
    epidemic_to_fit, 
    epid_periods,
    coverage_rates,
    care_delays,
    initial_parameters,
    n_samples, 
    nburn, 
    thinning,
    n_chains,
    txt_output = NULL
) {
  
  ## SET UP DATA FRAMES ETC. ##
  broad_ages <- data.table(
    age_grp = age_labels, 
    broad_age = fcn_assign_ages('children','adults','older_adults',age_labels)
  )
  epidemic_dt <- as.data.table(epidemic_to_fit)
  coverage_rates$imd_quintile <- factor(coverage_rates$imd_quintile)
  n_subtypes <- n_distinct(subtype_season$subtype)
  n_pars_per_subtype <- (length(initial_parameters) - 4)/2
  
  ll_call_count <- 0
  ll_total_calls <- (nburn + n_samples * thinning) * n_chains
  
  txt_out <- file.path('mcmc_output',paste0('index_',txt_output,'.txt'))
  
  # Define the log likelihood function
  llikelihood <- function(pars) {
    
    # Progress tracking
    if(ll_call_count == 0){ll_start_time <<- Sys.time()}
    mod_val <- if(ll_total_calls < 20){1}else{if(ll_total_calls < 1000){50}else{200}}
    ll_call_count <<- ll_call_count + 1
    if(ll_call_count %% mod_val == 0) {
      elapsed    <- as.numeric(difftime(Sys.time(), ll_start_time, units = 'mins'))
      rate       <- ll_call_count / elapsed
      remaining  <- (ll_total_calls - ll_call_count) / rate
      writeLines(sprintf(
        "INDEX %d: Iteration %d / %d (%.1f%%) | Elapsed: %.1f min | Est. remaining: %.1f min\n",
        txt_output, ll_call_count, ll_total_calls, 
        100 * ll_call_count / ll_total_calls,
        elapsed, remaining
      ), txt_out)
      cat()
    }
    
    epidemic_1 <- unique(subtype_season$subtype)[1]
    epidemic_2 <- ifelse(n_subtypes == 0, NA, unique(subtype_season$subtype)[2])
    
    ##----------------------------##
    #### PARAMETERS: EPIDEMIC 1 ####
    ##----------------------------##
    
    transmissibility_1 <- pars[1]

    adult_susceptibility_1 <- pars[2]
    rel_susceptibility_1 <- susc_vector(pars[3:4]) 
    susceptibility_vec_1 <- fcn_assign_ages(
      adult_susceptibility_1*rel_susceptibility_1[1],
      adult_susceptibility_1,
      adult_susceptibility_1*rel_susceptibility_1[2],
      age_labels
    )
 
    init_infected_1 <- 10^(pars[5])
    care_rate_df_1 <- data.frame(
      broad_age = rep(unique(broad_ages$broad_age), 2),
      risk_level = rep(c('low','high'), each = 3),
      primary_care = pars[6:11],
      secondary_care = pars[12:17]
    )
    care_rate_age_df_1 <- data.table(cross_join(
      care_rate_df_1,
      broad_ages))[broad_age.x == broad_age.y,]
    care_rate_age_df_1[, c('broad_age.x','broad_age.y') := NULL]
    
    imd_spline_pars <- data.table(
      primary = pars[1:2 + 2*n_pars_per_subtype],
      secondary = pars[3:4 + 2*n_pars_per_subtype]
    )
    rel_imd_rep_rates <- data.frame(imd_quintile = 1:5,
                                    rel_primary_rates = imd_spline(imd_spline_pars$primary),
                                    rel_secondary_rates = imd_spline(imd_spline_pars$secondary))
    
    reporting_rates_1 <- cross_join(
      care_rate_age_df_1,
      rel_imd_rep_rates
    ) %>% 
      mutate(primary_care = primary_care*rel_primary_rates,
             secondary_care = secondary_care*rel_secondary_rates) %>% 
      select(!c(rel_primary_rates, rel_secondary_rates))
    
    reporting_dt_1 <- as.data.table(reporting_rates_1)
    reporting_long_1 <- melt(reporting_dt_1,
                           measure.vars = c('primary_care', 'secondary_care'),
                           variable.name = 'setting',
                           value.name = 'rate')
    
    # any out-of-bounds proposals slipping past the prior
    if(any(is.na(pars)) || 
       transmissibility_1 < min_trans || 
       sum(susceptibility_vec_1 < min_susc) > 0  || sum(susceptibility_vec_1 > max_susc) > 0 ||
       init_infected_1 < 0 || init_infected_1 > min(demography_input$population) ) {
      return(-Inf)
    }
    
    if(n_subtypes == 2){
      
      ##----------------------------##
      #### PARAMETERS: EPIDEMIC 2 ####
      ##----------------------------##
      
      transmissibility_2 <- pars[1 + n_pars_per_subtype]
      
      adult_susceptibility_2 <- pars[2 + n_pars_per_subtype]
      rel_susceptibility_2 <- susc_vector(pars[3:4 + n_pars_per_subtype]) 
      susceptibility_vec_2 <- fcn_assign_ages(
        adult_susceptibility_2*rel_susceptibility_2[1],
        adult_susceptibility_2,
        adult_susceptibility_2*rel_susceptibility_2[2],
        age_labels
      )
      
      init_infected_2 <- 10^(pars[5 + n_pars_per_subtype])
      care_rate_df_2 <- data.frame(
        broad_age = rep(unique(broad_ages$broad_age), 2),
        risk_level = rep(c('low','high'), each = 3),
        primary_care = pars[6:11 + n_pars_per_subtype],
        secondary_care = pars[12:17 + n_pars_per_subtype]
      )
      care_rate_age_df_2 <- data.table(cross_join(
        care_rate_df_2,
        broad_ages))[broad_age.x == broad_age.y,]
      care_rate_age_df_2[, c('broad_age.x','broad_age.y') := NULL]
      
      reporting_rates_2 <- cross_join(
        care_rate_age_df_2,
        rel_imd_rep_rates
      ) %>% 
        mutate(primary_care = primary_care*rel_primary_rates,
               secondary_care = secondary_care*rel_secondary_rates) %>% 
        select(!c(rel_primary_rates, rel_secondary_rates))
      
      reporting_dt_2 <- as.data.table(reporting_rates_2)
      reporting_long_2 <- melt(reporting_dt_2,
                               measure.vars = c('primary_care', 'secondary_care'),
                               variable.name = 'setting',
                               value.name = 'rate')
      
      # any out-of-bounds proposals slipping past the prior
      if(any(is.na(pars)) || 
         transmissibility_2 < min_trans || 
         sum(susceptibility_vec_2 < min_susc) > 0  || sum(susceptibility_vec_2 > max_susc) > 0 ||
         init_infected_2 < 0 || init_infected_2 > min(demography_input$population) ) {
        return(-Inf)
      }
      
    }
    
    ##-----------------------##
    #### EPIDEMIC FUNCTION ####
    ##-----------------------##
    
    produce_epi_data <- function(k){
      
      epid_subtype <- get(paste0('epidemic_', k))
      
      vaccinated_subtype <- vaccinated_input[vaccinated_input$subtype == epid_subtype, ]
      pop_vaccinated <- vaccinated_subtype$vaccinated_population
      VE_INF <- vaccinated_subtype$VE_INF
      
      init_infected_vec <- (demography_input$population - pop_vaccinated)*get(paste0('init_infected_', k))/
        (sum(demography_input$population)-sum(pop_vaccinated))
      
      time_series <- run_model(
        pop = demography_input$population,
        I0 = init_infected_vec,
        vacc_cov = pop_vaccinated,
        ve_inf = VE_INF,
        cm = cm_input,
        trans = get(paste0('transmissibility_', k)),
        susc = get(paste0('susceptibility_vec_', k)),
        lat_per = epid_periods[1],
        inf_per = epid_periods[2]
      )
      
      ## add start date (using first of september throughout)
      start_of_epidemic <- as.Date(paste0('01-09-',year(epidemic_to_fit$week_start[1])), format = '%d-%m-%Y')
      time_series <- time_series[, date := start_of_epidemic + t] ## add date
      time_series[, t := NULL]
      
      if(any(is.na(time_series))) return(-Inf)
      
      weekly_infections <- time_series[, .(infections = sum(infections)), by = .(date)]$infections
      
      # Calculate value of OBSERVED infections 
      time_series <- time_series[coverage_rates, on = c('age_grp','imd_quintile','risk_level')]
      time_series_fit <- time_series[, .(infections = sum(OS_COVERAGE*infections)),
                                     by = .(date, age_grp, imd_quintile, risk_level, vaccinated)]
      time_series_fit[, imd_quintile := as.numeric(imd_quintile)]
      
      if(any(is.na(time_series_fit))) return(-Inf)
      
      setorder(time_series_fit, age_grp, imd_quintile, risk_level, vaccinated)
      
      # Aggregate to weekly and join with surveillance data
      time_series_fit[, week_start := last_monday(date)]
      time_series_weekly <- time_series_fit[, .(infections = sum(infections)), 
                                            by = .(week_start, age_grp, imd_quintile, risk_level, vaccinated)]
      
      # Join reporting rates (make primary and secondary care versions)
      reporting_data <- get(paste0('reporting_long_', k))
      time_series_long <- merge(rbind(cbind(time_series_weekly, setting = unique(reporting_data$setting)[1]),
                                      cbind(time_series_weekly, setting = unique(reporting_data$setting)[2])), 
                                reporting_data,
                                by = colnames(reporting_data)[colnames(reporting_data) != 'rate'],
                                all.x = TRUE)
      
      # Scale down hospitalisation rate for those ineffectively vaccinated
      # (VE_HOSP must be geq VI_INF)
      time_series_long <- merge(time_series_long, 
                                vaccinated_subtype[start_of_season == year(start_of_epidemic), 
                                                 c('age_grp','imd_quintile','risk_level','VE_INF','VE_HOSP')],
                                by = c('age_grp','imd_quintile','risk_level'),
                                all.x = TRUE)
      time_series_long[
        VE_HOSP >= VE_INF & vaccinated & setting == 'secondary_care', rate := rate*((1 - VE_HOSP)/(1 - VE_INF))
      ]
      
      time_series_long[, c('VE_HOSP','VE_INF') := NULL]
      
      # Validate reporting rates
      if(any(is.na(time_series_long$rate))) return(-Inf)
      if(any(time_series_long$rate <= 0 | time_series_long$rate >= 1)) return(-Inf)
      
      # Aggregate over vaccination status into a mean number of cases per {group x setting}
      # (as the rate varies between vaccinated/unvaccinated, but observations are not vaccination-specific)
      time_series_agg <- time_series_long[, .(expected_cases = sum(infections*rate)),
                                          by = list(age_grp, imd_quintile, risk_level, setting, week_start)
      ]
      
      return(list(weekly_infections, 
                  time_series_agg))
      
    }
    
    if(n_subtypes == 2){
      
      ##--------------------------------##
      #### RUNNING MULTIPLE EPIDEMICS ####
      ##--------------------------------##
    
      epi_data_1 <- produce_epi_data(1)
      epi_data_2 <- produce_epi_data(2)
      
      weekly_proportions <- data.table(
        date_formatted = as.Date(paste0('01-09-',year(epidemic_to_fit$week_start[1])), format = '%d-%m-%Y') +
          1:length(epi_data_1[[1]]),
        modelled_infections = (epi_data_1[[1]] + epi_data_2[[1]]),
        modelled_proportion = epi_data_1[[1]]/(epi_data_1[[1]] + epi_data_2[[1]])
        ) %>% drop_na()
      
      weekly_true_proportions <- subtype_season %>% filter(subtype == epidemic_1) %>% select(date_formatted, proportion)
      
      weekly_props <- weekly_true_proportions %>% 
        left_join(weekly_proportions, by = 'date_formatted') %>% drop_na()
      
      # Merge subtype-specific data and aggregate
      MODELLED_EPI <- rbind(epi_data_1[[2]], epi_data_2[[2]])[, .(expected_cases = sum(expected_cases)),
                                                                   by = list(age_grp, imd_quintile, risk_level, setting, week_start)]
      
    }else{
      
      ##-----------------------------##
      #### RUNNING SINGLE EPIDEMIC ####
      ##-----------------------------##
      
      epi_data_1 <- produce_epi_data(1)
      
      MODELLED_EPI <- epi_data_1[[2]]
      
    }
    
    # Pivot observed data longer
    observed_data <- melt(epidemic_dt, 
                          measure.vars = c('primary_care', 'secondary_care'),
                          variable.name = 'setting', 
                          value.name = 'observations')
    
    ## Shift observations by the care delays
    observed_data <- rbind(observed_data[setting=='primary_care'][,
                                                                  observations := shift(observations, n = -delays['primary'], fill = 0),
                                                                  by = .(age_grp, imd_quintile, risk_level, index)],
                           observed_data[setting=='secondary_care'][,
                                                                    observations := shift(observations, n = -delays['secondary'], fill = 0),
                                                                    by = .(age_grp, imd_quintile, risk_level, index)])
    
    # Merge with observed data
    time_series_joint <- merge(MODELLED_EPI, 
                               observed_data,
                               by = c('age_grp', 'imd_quintile', 'week_start', 'risk_level', 'setting'), 
                               all.x = F)
    
    if(F){
      time_series_joint %>% ## plot {observations} against {infections x reporting rates}
        ggplot() + theme_bw() +
        geom_line(aes(week_start,expected_cases, lty=risk_level, col=setting, group=interaction(setting, risk_level))) +
        geom_point(aes(week_start,observations, shape=risk_level, col=setting, group=interaction(setting, risk_level)), alpha = 0.4) +
        facet_grid(age_grp ~ imd_quintile, scales = 'free')
      }
    
    # OLD (FOR BINOMIAL LLIKELIHOOD):
    # Make log likelihood use max(infections, observations) to avoid -Inf where poss
    # I.e. if(infections[i] < observations[i]){infections[i] <- observations[i]}
    # time_series_maxed <- pmax(time_series_shifted$infections,
    #                           time_series_shifted$observations)
    
    # Vectorised log likelihood
    total_ll <- sum(dpois(
      x    = time_series_joint$observations,
      lambda = time_series_joint$expected_cases,
      log  = TRUE
    ), na.rm = TRUE)
    
    if(n_subtypes == 2){
      
      total_ll <- total_ll + sum(dbinom(
        x    = round(weekly_props$modelled_infections*weekly_props$modelled_proportion),
        size = round(weekly_props$modelled_infections),
        p    = weekly_props$proportion,
        log  = TRUE
      ), na.rm = TRUE)
      
    }
    
    if(is.nan(total_ll) | is.infinite(total_ll)) return(-Inf)
    
    return(total_ll)
  
  }
  
  llprior <- function(pars) {
    
    # pars[1] is transmissibility, pars[2:4] are susceptibility (children, adults, older adults), 
    # pars[5] is log of initial infected on 1st september
    
    susceptibility_vec_1 <- fcn_assign_ages(
      pars[2]*pars[3],
      pars[2],
      pars[2]*pars[4],
      age_labels
    )
    susceptibility_vec_2 <- fcn_assign_ages(
      pars[2 + n_pars_per_subtype]*pars[3 + n_pars_per_subtype],
      pars[2 + n_pars_per_subtype],
      pars[2 + n_pars_per_subtype]*pars[4 + n_pars_per_subtype],
      age_labels
    )
    
    reporting_vec <- c(6:n_pars_per_subtype, n_pars_per_subtype + 6:n_pars_per_subtype)
    reporting_prim_vec <- c(6:11, n_pars_per_subtype + 6:11)
    reporting_sec_vec <- reporting_vec[reporting_vec %notin% reporting_prim_vec]
    
    if(
      
      ## TRANSMISSIBILITY
      pars[1] < min_trans || pars[1] > max_trans ||
      pars[1 + n_pars_per_subtype] < min_trans || pars[1 + n_pars_per_subtype] > max_trans ||
      
      ## SUSCEPTIBILITY
      sum(c(susceptibility_vec_1, susceptibility_vec_2) < min_susc) > 0 || sum(c(susceptibility_vec_1, susceptibility_vec_2) > max_susc) > 0 ||
      
      ## INITIAL INFECTED
      pars[5] < min_log_init_inf || pars[5] > max_log_init_inf ||
      pars[5 + n_pars_per_subtype] < min_log_init_inf || pars[5 + n_pars_per_subtype] > max_log_init_inf ||
      
      ## REPORTING
      any(pars[reporting_vec] < min_reporting) || any(pars[reporting_vec] > max_reporting) ||
      
      ## IMD SPLINES
      any(pars[1:4 + 2*n_pars_per_subtype] < min_spline) || any(pars[1:4 + 2*n_pars_per_subtype] > max_spline)
      
    ) {return(-Inf)}
    
    lprob <- 0
    
    R0_1 <- R0_func(susceptibility = susceptibility_vec_1,
                    inf_period = epid_periods[2],
                    beta_in = pars[1],
                    cm_in = cm_input,
                    per_capita = T,
                    population_vector = demography_input$population)
    
    R0_2 <- R0_func(susceptibility = susceptibility_vec_2,
                    inf_period = epid_periods[2],
                    beta_in = pars[1 + n_pars_per_subtype],
                    cm_in = cm_input,
                    per_capita = T,
                    population_vector = demography_input$population)
    
    if(sum(c(R0_1, R0_2) < min_r0) > 0 || sum(c(R0_1, R0_2) > max_r0) > 0) { return(-Inf) }
    
    # Scaled normal prior on R0 (pnorm(min_r0) is the log normalising constant)
    lprob <- lprob + dnorm(R0_1, mean = r0_mean, sd = r0_sd, log = TRUE) -
      pnorm(min_r0, mean = r0_mean, sd = r0_sd, lower.tail = FALSE, log.p = TRUE)
    # Uniform prior on transmissibility 
    lprob <- lprob + dunif(unname(pars[1]), min = min_trans, max = max_trans, log = TRUE)
    # Uniform prior on susceptibility, x2
    lprob <- lprob + sum(dunif(unname(pars[2])*c(1, unname(pars[3:4])), min = min_susc, max = max_susc, log = TRUE))
    # Uniform prior on initial infected
    lprob <- lprob + dunif(unname(pars[5]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
    # Beta prior on primary care reporting rates (pars 5:10), centred at 2%
    lprob <- lprob + sum(dbeta(pars[6:11], 
                               shape1 = prim_beta['a'], 
                               shape2 = prim_beta['b'], log = TRUE))
    # Beta prior on secondary care reporting rates (pars 11:16), centred at 0.5%
    lprob <- lprob + sum(dbeta(pars[12:17], 
                               shape1 = sec_beta['a'], 
                               shape2 = sec_beta['b'], log = TRUE))
    
    if(n_subtypes == 2){
      
      # priors on second set of epi parameters
      
      lprob <- lprob + dnorm(R0_2, mean = r0_mean, sd = r0_sd, log = TRUE) -
        pnorm(min_r0, mean = r0_mean, sd = r0_sd, lower.tail = FALSE, log.p = TRUE)
      
      lprob <- lprob + dunif(unname(pars[1 + n_pars_per_subtype]), min = min_trans, max = max_trans, log = TRUE)
      
      lprob <- lprob + sum(dunif(unname(pars[2 + n_pars_per_subtype])*c(1, unname(pars[3:4 + n_pars_per_subtype])), min = min_susc, max = max_susc, log = TRUE))
      
      lprob <- lprob + dunif(unname(pars[5 + n_pars_per_subtype]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
      
      lprob <- lprob + sum(dbeta(pars[6:11 + n_pars_per_subtype], 
                                 shape1 = prim_beta['a'], 
                                 shape2 = prim_beta['b'], log = TRUE))
      
      lprob <- lprob + sum(dbeta(pars[12:17 + n_pars_per_subtype], 
                                 shape1 = sec_beta['a'], 
                                 shape2 = sec_beta['b'], log = TRUE))
      
    }
    
    # Uniform prior on spline parameters
    lprob <- lprob + sum(dunif(unname(pars[1:4 + 2*n_pars_per_subtype]), min = min_spline, max = max_spline, log = TRUE))
    
    return(lprob)
  }
  
  # get beta parameters from mean and concentration
  beta_pars <- function(mean, concentration) {
    c(a = mean * concentration, b = (1 - mean) * concentration)
  }
  # TODO Should these be lognormal instead?
  
  # Primary care: centred at 0.02, secondary: centred at 0.005
  prim_beta  <- beta_pars(0.02, 200) 
  sec_beta   <- beta_pars(0.005, 200)
  # R0 centred at 2 
  r0_mean <- 2
  r0_sd   <- 0.4
  
  ## set bounds
  min_trans <- 0; max_trans <- 1
  min_susc <- 0; max_susc <- 1
  min_r0 <- 1; max_r0 <- 4
  min_log_init_inf <- 0; max_log_init_inf <- log10(min(demography_input$population))
  min_reporting <- 0; max_reporting <- 1
  min_spline <- -log(5); max_spline <- log(5) # equivalent to IMD ratios at most 
  # 5x higher/lower than each other in IMD 1 vs 3 or IMD 5 vs 3,
  # which would lead to huge ratios between IMD 1 and IMD 5 at the extremes
  
  ## set up sampler
  lower_vals <- c(min_trans, rep(min_susc, 3), min_log_init_inf,
                  rep(min_reporting, 12),
                  min_trans, rep(min_susc, 3), min_log_init_inf,
                  rep(min_reporting, 12),
                  rep(min_spline, 4))
  upper_vals <- c(max_trans, rep(max_susc, 3), max_log_init_inf,
                  rep(max_reporting, 12),
                  max_trans, rep(max_susc, 3), max_log_init_inf,
                  rep(max_reporting, 12),
                  rep(max_spline, 4))
  
  sampler <- function(n = 1){
    out <- matrix(NA, nrow = n, ncol = length(initial_parameters))
    for(j in 1:n){
      valid <- FALSE
      while(!valid){
        repeat {
          R0_1 <- rnorm(1, mean = r0_mean, sd = r0_sd)
          R0_2 <- rnorm(1, mean = r0_mean, sd = r0_sd)
          if (R0_1 >= min_r0 && R0_1 <= max_r0 && R0_2 >= min_r0 && R0_2 <= max_r0) break
        }
        
        susc_11 <- runif(1, min_susc, max_susc)
        susc_21 <- runif(1, min_susc, max_susc)
        susc_31 <- runif(1, min_susc, max_susc)
        
        adult_susc_1 <- susc_21
        child_rel_susc_1 <- susc_11/susc_21
        older_adult_rel_susc_1 <- susc_31/susc_21
        
        susc_12 <- runif(1, min_susc, max_susc)
        susc_22 <- runif(1, min_susc, max_susc)
        susc_32 <- runif(1, min_susc, max_susc)
        
        adult_susc_2 <- susc_22
        child_rel_susc_2 <- susc_12/susc_22
        older_adult_rel_susc_2 <- susc_32/susc_22
        
        trans_1 <- R0_func(
          susceptibility    = fcn_assign_ages(susc_11, susc_21, susc_31, age_labels),
          inf_period        = epid_periods[2],
          beta_in           = 1,
          cm_in             = cm_input,
          per_capita        = TRUE,
          population_vector = demography_input$population,
          R0assumed         = R0_1,
          return_beta       = TRUE
        )
        
        trans_2 <- R0_func(
          susceptibility    = fcn_assign_ages(susc_12, susc_22, susc_32, age_labels),
          inf_period        = epid_periods[2],
          beta_in           = 1,
          cm_in             = cm_input,
          per_capita        = TRUE,
          population_vector = demography_input$population,
          R0assumed         = R0_2,
          return_beta       = TRUE
        )
        
        if(trans_1 >= min_trans & trans_1 <= max_trans & 
           trans_2 >= min_trans & trans_2 <= max_trans) valid <- TRUE
      }
      
      # Check trans is within bounds
      if(trans_1 < min_trans | trans_1 > max_trans |
         trans_2 < min_trans | trans_2 > max_trans) next
      
      # Draw reporting rates from Beta priors
      prim_rates_1 <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
      sec_rates_1  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
      prim_rates_2 <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
      sec_rates_2  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
      
      out[j, ] <- c(
        trans_1, adult_susc_1, child_rel_susc_1, older_adult_rel_susc_1,
        runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
        prim_rates_1,                                     # primary care rates
        sec_rates_1,                                      # secondary care rates
        trans_2, adult_susc_2, child_rel_susc_2, older_adult_rel_susc_2,
        runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
        prim_rates_2,                                     # primary care rates
        sec_rates_2,                                      # secondary care rates
        runif(4, min_spline, max_spline)                # IMD spline, uniform
      )
    }
    return(out)
  }
  
  ## prior
  prior <- createPrior(density = llprior, 
                       sampler = sampler,
                       lower = lower_vals,
                       upper = upper_vals)
  
  bayesianSetup <- createBayesianSetup(
    likelihood = llikelihood, 
    prior = prior
  )
  
  settings <- list(
    iterations = nburn + n_samples*thinning, ## setup to save all (pre-thinning etc.)
    burnin = 0,
    thin = 1,
    message = T, nrChains=n_chains, parallel = F
  ) 
  
  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = 'DEzs', settings = settings)
  
  return(out)
}


## FUNCTIONS TO PLOT MCMC SAMPLES ##

plot_density <- function(var, filtered = T){
  
  data <- if(filtered){mcmc_samples_filtered}else{mcmc_samples}
  
  if(var != 'likelihood'){
    data <- data %>% select(!likelihood)
  }
  
  epid_pars_joining <- epid_pars %>% 
    mutate(true_value = value) %>% 
    filter(name == var) %>% 
    select(epidemic, epidemic_of_season, season, subtype, name, true_value)
  
  if(!grepl('imd_spline', var)){
    
    data %>%
      select(!starts_with('imd_')) %>% 
      pivot_longer(!c(iteration,epidemic,chain)) %>%
      mutate(epidemic_of_season = as.numeric(substr(name, nchar(name), nchar(name))),
             name = substr(name, 1, nchar(name) - 7)) %>% 
      filter(name == var) %>%
      left_join(epid_pars_joining, by = c('epidemic','epidemic_of_season','name')) %>% 
      filter(!is.na(season)) %>% # remove non-real epidemics
      mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                       T ~ F)) %>% 
      mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
      mutate(epidemic = paste0("Season ", epidemic),
             subtype_season = paste0(season, ": ", subtype)) %>% 
      ggplot() +
        geom_density(aes(x = value, fill = as.factor(chain), group = chain), alpha = 0.4) +
        geom_vline(aes(xintercept = true_value, alpha = true_value_NA), lty=2) +
        scale_alpha_manual(values = c(1,0)) +
        theme_bw() + labs(y = ifelse(var=='R0', 'R0 (calculated after)', var)) +
        facet_grid(.~subtype_season, scales = 'free') +
        theme(legend.position = 'none')
    
  }else{
    
    epid_pars_joining <- epid_pars_joining %>% 
      select(!c(subtype, epidemic_of_season)) %>% unique()
    
    data %>%
      pivot_longer(!c(iteration,epidemic,chain)) %>%
      filter(name == var) %>%
      left_join(epid_pars_joining, by = c('epidemic','name')) %>% 
      filter(!is.na(season)) %>% # remove non-real epidemics
      mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                       T ~ F)) %>% 
      mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
      mutate(epidemic = paste0("Season ", epidemic)) %>% 
      ggplot() +
      geom_density(aes(x = value, fill = as.factor(chain), group = chain), alpha = 0.4) +
      geom_vline(aes(xintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      theme_bw() + labs(y = ifelse(var=='R0', 'R0 (calculated after)', var)) +
      facet_grid(.~season, scales = 'free') +
      theme(legend.position = 'none')
    
  }
    
}

plot_trace <- function(var, filtered = F){
  
  data <- if(filtered){mcmc_samples_filtered}else{mcmc_samples} 
  
  if(var != 'likelihood'){
    data <- data %>% select(!likelihood)
  }
  
  Y_LAB <- if(var=='R0'){'R0 (calculated after)'}else{
    if(var=='init_infected'){'Initial infected (log10)'}else{var}
  }
  
  epid_pars_joining <- epid_pars %>% 
    mutate(true_value = value) %>% 
    filter(name == var) %>% 
    select(epidemic, epidemic_of_season, season, subtype, name, true_value)
  
  if(!grepl('imd_spline', var)){
    
    p <- data %>%
      select(!starts_with('imd_')) %>% 
      pivot_longer(!c(iteration,epidemic,chain)) %>%
      mutate(epidemic_of_season = as.numeric(substr(name, nchar(name), nchar(name))),
             name = substr(name, 1, nchar(name) - 7)) %>% 
      filter(name == var) %>%
      left_join(epid_pars_joining, by = c('epidemic','epidemic_of_season','name')) %>% 
      filter(!is.na(season)) %>% # remove non-real epidemics
      mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                       T ~ F)) %>% 
      mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
      mutate(epidemic = paste0("Season ", epidemic),
             subtype_season = paste0(season, ": ", subtype)) %>% 
      ggplot() +
      geom_line(aes(x = iteration/1000, y = value, col = as.factor(chain), group = chain)) +
      geom_hline(aes(yintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      # geom_vline(xintercept = burn_in, lty=3, alpha = 0.5) +
      facet_grid(.~subtype_season) +
      theme_bw() + labs(y = Y_LAB) +
      # scale_color_manual(values = var_cols) +
      theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
    
  }else{
    
    epid_pars_joining <- epid_pars_joining %>% 
      select(!c(subtype, epidemic_of_season)) %>% unique()
    
    p <- data %>%
      pivot_longer(!c(iteration,epidemic,chain)) %>%
      filter(name == var) %>%
      left_join(epid_pars_joining, by = c('epidemic','name')) %>% 
      filter(!is.na(season)) %>% # remove non-real epidemics
      mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                       T ~ F)) %>% 
      mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
      mutate(epidemic = paste0("Season ", epidemic)) %>% 
      ggplot() +
      geom_line(aes(x = iteration/1000, y = value, col = as.factor(chain), group = chain)) +
      geom_hline(aes(yintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      # geom_vline(xintercept = burn_in, lty=3, alpha = 0.5) +
      facet_grid(.~season) +
      theme_bw() + labs(y = Y_LAB) +
      # scale_color_manual(values = var_cols) +
      theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
    
  }
  
  if(var=='init_infected'){
    p <- p + scale_y_log10()
  }
  
  p
}


## OLD FITTING FUNCTION, WHERE REPORTING RATES WERE KNOWN ##
#
# run_mcmc_inference_old <- function(
#     demography_input, 
#     vaccinated_input,
#     cm_input, 
#     epidemic_to_fit, 
#     epid_periods,
#     reporting_rates,
#     coverage_rates,
#     care_delays,
#     initial_parameters,
#     n_samples, 
#     nburn, 
#     thinning,
#     n_chains
# ) {
#   
#   ll_call_count <- 0
#   ll_total_calls <- (nburn + n_samples * thinning) * n_chains
#   
#   txt_out <- file.path('mcmc_output',paste0('index_',txt_output,'.txt'))
#   
#   # Define the log likelihood function
#   llikelihood <- function(pars) {
#     
#     # Progress tracking
#     if(ll_call_count == 0){ll_start_time <<- Sys.time()}
#     mod_val <- if(ll_total_calls < 20){1}else{if(ll_total_calls < 1000){50}else{200}}
#     ll_call_count <<- ll_call_count + 1
#     if(ll_call_count %% mod_val == 0) {
#       elapsed    <- as.numeric(difftime(Sys.time(), ll_start_time, units = 'mins'))
#       rate       <- ll_call_count / elapsed
#       remaining  <- (ll_total_calls - ll_call_count) / rate
#       writeLines(sprintf(
#         "INDEX %d: Iteration %d / %d (%.1f%%) | Elapsed: %.1f min | Est. remaining: %.1f min\n",
#         txt_output, ll_call_count, ll_total_calls, 
#         100 * ll_call_count / ll_total_calls,
#         elapsed, remaining
#       ), txt_out)
#       cat()
#     }
#     
#     transmissibility <- pars[1]
#     susceptibility <- susc_vector(pars[2:3]) 
#     init_infected_num <- 10^(pars[4])
#     
#     # any out-of-bounds proposals slipping past the prior
#     if(any(is.na(pars)) || 
#        transmissibility < min_trans || 
#        sum(susceptibility < min_susc) > 0  || sum(susceptibility > max_susc) > 0 ||
#        init_infected_num < 0 || init_infected_num > min(demography_input$population) ) {
#       return(-Inf)
#     }
#     
#     pop_vaccinated <- vaccinated_input$effectively_vaccinated_population
#     
#     init_infected_vec <- (demography_input$population - pop_vaccinated)*init_infected_num/
#       (sum(demography_input$population)-sum(pop_vaccinated))
#     
#     time_series <- run_model(
#       pop = demography_input$population,
#       I0 = init_infected_vec,
#       vacc = pop_vaccinated,
#       cm = cm_input,
#       trans = transmissibility,
#       susc = susceptibility,
#       lat_per = epid_periods[1],
#       inf_per = epid_periods[2]
#     )
#     
#     ## add start date (using first of september throughout)
#     start_of_epidemic <- as.Date(paste0('01-09-',year(epidemic_to_fit$week_start[1])), format = '%d-%m-%Y')
#     time_series <- time_series[, date := start_of_epidemic + t] ## add date
#     time_series[, t := NULL]
#     
#     if(any(is.na(time_series))) return(-Inf)
#     
#     # Take rounded value of OBSERVED infections (need an integer)
#     coverage_rates$imd_quintile <- factor(coverage_rates$imd_quintile)
#     time_series <- time_series[coverage_rates, on = c('age_grp','imd_quintile','risk_level')]
#     time_series_fit <- time_series[, .(infections = round(sum(OS_COVERAGE*infections))),
#                                         by = .(date, age_grp, imd_quintile, risk_level)]
#     time_series_fit[, imd_quintile := as.numeric(imd_quintile)]
#     
#     if(any(is.na(time_series_fit))) return(-Inf)
#     
#     # Check epidemic is growing at start, and not at the end
#     setorder(time_series_fit, age_grp, imd_quintile, risk_level)
#     if(time_series_fit$infections[1] > time_series_fit$infections[5]) return(-Inf)
#     if(time_series_fit$infections[nrow(time_series_fit)] > time_series_fit$infections[nrow(time_series_fit)-5]) return(-Inf)
#     
#     # Aggregate to weekly and join with surveillance data
#     time_series_fit[, week_start := last_monday(date)]
#     time_series_weekly <- time_series_fit[, .(infections = sum(infections)), 
#                                           by = .(week_start, age_grp, imd_quintile, risk_level)]
#     
#     # Convert epidemic_to_fit to data.table if not already
#     epidemic_dt <- as.data.table(epidemic_to_fit)
#     
#     # Join with observed data
#     time_series_joint <- merge(time_series_weekly, epidemic_dt, 
#                                by = c('age_grp', 'imd_quintile', 'week_start', 'risk_level'), 
#                                all.x = F)
#     
#     # Pivot primary_care and secondary_care to long format
#     time_series_long <- melt(time_series_joint, 
#                              measure.vars = c('primary_care', 'secondary_care'),
#                              variable.name = 'setting', 
#                              value.name = 'observations')
#     
#     # Join reporting rates
#     reporting_dt <- as.data.table(reporting_rates)
#     reporting_long <- melt(reporting_dt,
#                            measure.vars = c('primary_care', 'secondary_care'),
#                            variable.name = 'setting',
#                            value.name = 'rate')
#     
#     time_series_long <- merge(time_series_long, reporting_long,
#                               by = colnames(reporting_long)[colnames(reporting_long) != 'rate'],
#                               all.x = TRUE)
#     
#     # Validate reporting rates
#     if(any(is.na(time_series_long$rate))) return(-Inf)
#     if(any(time_series_long$rate <= 0 | time_series_long$rate >= 1)) return(-Inf)
#     
#     ## shift observations by the care delays
#     time_series_shifted <- rbind(
#       time_series_long[setting=='primary_care'][,
#       observations := shift(observations, n = -delays['primary'], fill = 0),
#       by = .(age_grp, imd_quintile, risk_level, index, rate)],
#       time_series_long[setting=='secondary_care'][,
#       observations := shift(observations, n = -delays['secondary'], fill = 0),
#       by = .(age_grp, imd_quintile, risk_level, index, rate)])
# 
#     # time_series_shifted %>% ## plot {observations} against {infections x reporting rates}
#     #   ggplot() + theme_bw() +
#     #   geom_line(aes(week_start,rate*infections, lty=risk_level, group=interaction(setting, risk_level))) +
#     #   geom_line(aes(week_start,observations, lty=risk_level, group=interaction(setting, risk_level)),col=2) +
#     #   facet_grid(age_grp ~ imd_quintile, scales = 'free')
# 
#     # Vectorised log likelihood
#     total_ll <- sum(dbinom(
#       x    = time_series_shifted$observations,
#       size = time_series_shifted$infections,
#       prob = time_series_shifted$rate,
#       log  = TRUE
#     ), na.rm = TRUE)
#     
#     if(is.nan(total_ll) | is.infinite(total_ll)) return(-Inf)
#     
#     return(total_ll)
#   }
#   
#   llprior <- function(pars) {
#     
#     # pars[1] is transmissibility, pars[2:3] are susceptibility (children, adults, older adults), 
#     # pars[4] is log of initial infected on 1st september
#     
#     if(
#       pars[1] < min_trans ||
#       pars[2] < min_susc || pars[2] > max_susc ||
#       pars[3] < min_susc || pars[3] > max_susc ||
#       pars[4] < min_log_init_inf || pars[4] > max_log_init_inf 
#     ) {return(-Inf)}
#     
#     lprob <- 0
#     
#     R0 <- R0_func(susceptibility = susc_vector(pars[2:3]),
#                   inf_period = epid_periods[2],
#                   beta_in = pars[1],
#                   cm_in = cm_input,
#                   per_capita = T,
#                   population_vector = demography_input$population)
#     
#     if(R0 < min_r0 || R0 > max_r0) { return(-Inf) }
#     
#     # Uniform prior on R0 between 1 and 3
#     lprob <- lprob + dunif(R0, min = min_r0, max = max_r0, log = TRUE)
#     # Uniform prior on transmissibility
#     lprob <- lprob + dunif(unname(pars[1]), min = min_trans, max = max_trans, log = TRUE)
#     # Uniform prior on susceptibility, x3
#     lprob <- lprob + dunif(unname(pars[2]), min = min_susc, max = max_susc, log = TRUE)
#     lprob <- lprob + dunif(unname(pars[3]), min = min_susc, max = max_susc, log = TRUE)
#     lprob <- lprob + dunif(unname(pars[4]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
# 
#     return(lprob)
#   }
#   
#   ## set bounds
#   min_trans <- 0; max_trans <- 0.4
#   min_susc <- 0; max_susc <- 5
#   min_r0 <- 1; max_r0 <- 3
#   min_log_init_inf <- 0; max_log_init_inf <- log10(min(demography_input$population))
#   
#   ## set up sampler
#   lower_vals <- c(min_trans, rep(min_susc, 2), min_log_init_inf)
#   upper_vals <- c(max_trans, rep(max_susc, 2), max_log_init_inf)
#   
#   sampler <- function(n = 1){
#     out <- matrix(NA, nrow = n, ncol = length(initial_parameters))
#     for(j in 1:n){
#       R0 <- 0
#       while(R0 < 1 | R0 > 3){
#         trans <- runif(1, lower_vals[1], upper_vals[1])
#         susc_1  <- runif(1, lower_vals[2], upper_vals[2])
#         susc_2  <- runif(1, lower_vals[3], upper_vals[3])
#         R0 <- R0_func(
#           susceptibility    = susc_vector(c(susc_1, susc_2)),
#           inf_period        = epid_periods[2],
#           beta_in           = trans,
#           cm_in             = cm_input,
#           per_capita        = TRUE,
#           population_vector = demography_input$population
#         )
#       }
#       out[j, ] <- c(trans, susc_1, susc_2, runif(1, lower_vals[4], upper_vals[4]))
#     }
#     return(out)
#   }
#   
#   ## prior
#   prior <- createPrior(density = llprior, 
#                        sampler = sampler,
#                        lower = lower_vals,
#                        upper = upper_vals)
#   
#   bayesianSetup <- createBayesianSetup(
#     likelihood = llikelihood, 
#     prior = prior
#   )
#   
#   settings <- list(
#     iterations = nburn + n_samples*thinning, 
#     burnin = 0,
#     thin = 1,
#     message = T, nrChains=n_chains, parallel = F
#   ) ## setup to save all (pre-thinning etc.)
#   
#   # cat(llprior(initial_parameters))
#   
#   out <- runMCMC(bayesianSetup = bayesianSetup, sampler = 'DEzs', settings = settings)
#   
#   #plot(out); summary(out)
#   
#   return(out)
# }



