## MCMC FUNCTIONS ##

TEXT_SAVE_DIR <- file.path('mcmc_output')
if(!dir.exists(TEXT_SAVE_DIR)){dir.create(TEXT_SAVE_DIR)}

## MCMC FITTING ##
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
    n_pop = 100,
    n_cores = 1,
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
  ll_total_calls <- (nburn + n_samples * thinning)/n_pop
  
  txt_out <- file.path('mcmc_output',paste0('index_',txt_output,'.txt'))
  
  LLcache <- new.env(hash = TRUE, parent = emptyenv())
  
  make_key <- function(par) paste(sprintf("%.15g", par), collapse = "_")
  
  # Define the log likelihood function
  single_ll <- function(pars) {
    
    epidemic_1 <- unique(subtype_season$subtype)[1]
    epidemic_2 <- ifelse(n_subtypes == 0, NA, unique(subtype_season$subtype)[2])
    
    ##----------------------------##
    #### PARAMETERS: EPIDEMIC 1 ####
    ##----------------------------##
    
    transmissibility_1 <- pars[1]
     
    susceptibility_vec_1 <- fcn_assign_ages(
      pars[2],
      pars[3],
      pars[4],
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
      
      susceptibility_vec_2 <- fcn_assign_ages(
        pars[2 + n_pars_per_subtype],
        pars[3 + n_pars_per_subtype],
        pars[4 + n_pars_per_subtype],
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
      start_of_epidemic <- as.Date(paste0('01-09-',year(epidemic_dt$week_start[1])), format = '%d-%m-%Y')
      time_series <- time_series[, date := start_of_epidemic + t] ## add date
      time_series[, t := NULL]
      
      if(any(is.na(time_series))) return(-Inf)
      
      weekly_infections <- time_series[, .(infections = sum(infections)), by = .(date)]$infections
      
      # Calculate value of OBSERVED infections 
      coverage_rates$imd_quintile <- factor(coverage_rates$imd_quintile,
                                            levels = unique(coverage_rates$imd_quintile))
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
      ) %>% drop_na() %>% 
        ## IF 0, CHANGE TO VERY LOW NUMBER TO AVOID -INF LOG-LIKELIHOOD
        ## IF 1, CHANGE TO VERY HIGH NUMBER TO AVOID -INF LOG-LIKELIHOOD
        mutate(
          min_modelled_proportion = min(modelled_proportion[modelled_proportion != 0]),
          max_modelled_proportion = max(modelled_proportion[modelled_proportion != 1]),
          modelled_proportion = case_when(
            modelled_proportion == 0 ~ min_modelled_proportion/10000,
            modelled_proportion == 1 ~ 1 - (1 - max_modelled_proportion)/10000,
            T ~ modelled_proportion
        )) %>% select(!c(min_modelled_proportion, max_modelled_proportion))
      
      weekly_true_proportions <- subtype_season %>% filter(subtype == epidemic_1) %>% 
        rename(epidemic_1_positive = value) %>% 
        select(date_formatted, proportion, total_flu, epidemic_1_positive)
      
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
    
    # Vectorised log likelihood
    ll_1 <- sum(dpois(
      x    = time_series_joint$observations,
      lambda = time_series_joint$expected_cases,
      log  = TRUE
    ), na.rm = TRUE)
    
    if(n_subtypes == 2){
      
      ll_2 <- sum(dbinom(
        x    = weekly_props$epidemic_1_positive,
        size = weekly_props$total_flu,
        p    = weekly_props$modelled_proportion,
        log  = TRUE
      ), na.rm = TRUE)
      
    }else{ ll_2 <- 0 }
    
    assign(make_key(pars), c(ll_1, ll_2), envir = LLcache)
    
    total_ll <- ll_1 + ll_2
    
    if(is.nan(total_ll) | is.infinite(total_ll)) return(-Inf)
    
    return(total_ll)
    
  }
  
  llikelihood <- function(pars) {
    
    if (is.matrix(pars)) {
      
      # progress tracking now happens ONCE PER GENERATION, in the parent process
      if (ll_call_count == 0) { ll_start_time <<- Sys.time() }
      ll_call_count <<- ll_call_count + 1
      mod_val <- if (ll_total_calls < 20) 1 else if (ll_total_calls < 1000) 10 else 50
      if (ll_call_count %% mod_val == 0) {
        elapsed   <- as.numeric(difftime(Sys.time(), ll_start_time, units = 'mins'))
        rate      <- ll_call_count / elapsed
        remaining <- (ll_total_calls - ll_call_count) / rate
        writeLines(sprintf(
          "INDEX %d: Generation %d / %d (%.1f%%) | Elapsed: %.1f min | Est. remaining: %.1f min\n",
          txt_output, ll_call_count, ll_total_calls,
          100 * ll_call_count / ll_total_calls, elapsed, remaining
        ), txt_out)
      }
      
      res <- parallel::mclapply(
        seq_len(nrow(pars)),
        function(j) single_ll(pars[j, ]),
        mc.cores = n_cores
      )
      return(unlist(res))
      
    } else {
      return(single_ll(pars))
    }
    
  }
  
  llprior <- function(pars) {
    
    # pars[1] is transmissibility, pars[2:4] are susceptibility (children, adults, older adults), 
    # pars[5] is log of initial infected on 1st september
    
    susceptibility_vec_1 <- fcn_assign_ages(
      pars[2],
      pars[3],
      pars[4],
      age_labels
    )
    susceptibility_vec_2 <- fcn_assign_ages(
      pars[2 + n_pars_per_subtype],
      pars[3 + n_pars_per_subtype],
      pars[4 + n_pars_per_subtype],
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
    # Beta prior on susceptibility
    lprob <- lprob + sum(dbeta(unname(pars[2:4]), shape1 = shape1_susc, shape2 = shape2_susc, log = TRUE))
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
      
      lprob <- lprob + sum(dbeta(unname(pars[2:4 + n_pars_per_subtype]), shape1 = shape1_susc, shape2 = shape2_susc, log = TRUE))
      
      lprob <- lprob + dunif(unname(pars[5 + n_pars_per_subtype]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
      
      lprob <- lprob + sum(dbeta(pars[6:11 + n_pars_per_subtype], 
                                 shape1 = prim_beta['a'], 
                                 shape2 = prim_beta['b'], log = TRUE))
      
      lprob <- lprob + sum(dbeta(pars[12:17 + n_pars_per_subtype], 
                                 shape1 = sec_beta['a'], 
                                 shape2 = sec_beta['b'], log = TRUE))
      
    }
    
    # Normal prior on spline parameters
    lprob <- lprob + sum(dnorm(unname(pars[1:4 + 2*n_pars_per_subtype]), mean = imd_mean, sd = imd_sd, log = TRUE))
    
    return(lprob)
  }
  
  # get beta parameters from mean and concentration
  beta_pars <- function(mean, concentration) {
    c(a = mean * concentration, b = (1 - mean) * concentration)
  }
  # TODO Should these be lognormal instead?
  
  # set up priors and lower/upper bounds
  {
  # Primary care: centred at 0.02, secondary: centred at 0.005
  prim_beta  <- beta_pars(0.02, 200) 
  sec_beta   <- beta_pars(0.005, 200)
  # Susceptibility beta distributed, mean at 0.5
  shape1_susc <- 2
  shape2_susc <- 2
  # R0 centred at 2 
  r0_mean <- 2
  r0_sd   <- 0.4
  # IMD splines normal around 0
  imd_mean <- 0
  imd_sd <- 0.5
  
  ## set bounds
  min_trans <- 0; max_trans <- 1
  min_susc <- 0; max_susc <- 1
  min_r0 <- 1; max_r0 <- 4
  min_log_init_inf <- 0; max_log_init_inf <- 4
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
  }
  
  sampler <- function(n = 1){
    out <- matrix(NA, nrow = n, ncol = length(initial_parameters))
    for(j in 1:n){
      repeat {
      valid <- FALSE
      while(!valid){
        repeat {
          R0_1 <- rnorm(1, mean = r0_mean, sd = r0_sd)
          R0_2 <- rnorm(1, mean = r0_mean, sd = r0_sd)
          if (R0_1 >= min_r0 && R0_1 <= max_r0 && R0_2 >= min_r0 && R0_2 <= max_r0) break
        }
        
        susc_1 <- rbeta(3, shape1_susc, shape2_susc)
        susc_2 <- rbeta(3, shape1_susc, shape2_susc)
        
        trans_1 <- R0_func(
          susceptibility    = fcn_assign_ages(susc_1[1],
                                              susc_1[2],
                                              susc_1[3], 
                                              age_labels),
          inf_period        = epid_periods[2],
          beta_in           = 1,
          cm_in             = cm_input,
          per_capita        = TRUE,
          population_vector = demography_input$population,
          R0assumed         = R0_1,
          return_beta       = TRUE
        )
        
        trans_2 <- R0_func(
          susceptibility    = fcn_assign_ages(susc_2[1],
                                              susc_2[2],
                                              susc_2[3],  
                                              age_labels),
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
         trans_2 < min_trans | trans_2 > max_trans ) next
      
      ## IMD splines normal around 0, SD = 0.5
      imd_spline_samples <- rnorm(4, mean = imd_mean, sd = imd_sd) 
      
      # Check IMD splines are within bounds
      if(sum(imd_spline_samples < min_spline) + 
         sum(imd_spline_samples > max_spline) > 0) next
      
      # Draw reporting rates from Beta priors
      prim_rates_1 <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
      sec_rates_1  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
      prim_rates_2 <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
      sec_rates_2  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
      
      out[j, ] <- c(
        trans_1, susc_1,
        runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
        prim_rates_1,                                     # primary care rates
        sec_rates_1,                                      # secondary care rates
        trans_2, susc_2,
        runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
        prim_rates_2,                                     # primary care rates
        sec_rates_2,                                      # secondary care rates
        imd_spline_samples       # IMD spline, normal around 0
      )
      break
      }
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
    prior = prior,
    parallel = "external"
  )
  
  settings <- list(
    iterations = nburn + n_samples*thinning, ## setup to save all (pre-thinning etc.)
    burnin = 0,
    thin = 1,
    message = T, 
    nrChains=n_chains,
    startValue = prior$sampler(n_pop) 
  ) 
  
  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = 'DEzs', settings = settings)
  
  nPar <- out$setup$numPars
  
  lookupComponents <- function(chainMat, cache, nPar) {
    parMat <- chainMat[, 1:nPar, drop = FALSE]
    comp <- t(apply(parMat, 1, function(par) {
      val <- mget(make_key(par), envir = cache, ifnotfound = list(c(NA_real_, NA_real_)))[[1]]
      val
    }))
    colnames(comp) <- c("LL1", "LL2")
    coda::mcmc(cbind(chainMat, comp))
  }
  
  LLcomponents <- if (coda::is.mcmc.list(out$chain)) {
    coda::as.mcmc.list(lapply(out$chain, lookupComponents, cache = LLcache, nPar = nPar))
  } else {
    lookupComponents(out$chain, LLcache, nPar)
  }
  
  out$LLcomponents <- LLcomponents
  
  return(out)
}


## FUNCTIONS TO PLOT MCMC SAMPLES ##

plot_density <- function(var, filtered = T){
  
  data <- if(filtered){mcmc_samples_filtered}else{mcmc_samples}
  
  if(var != 'likelihood'){
    data <- data %>% select(!likelihood)
  }
  
  var_label <- gsub('_rate_','_rate\n', var)
  var_label <- gsub('_spline_','_spline\n', var_label)
  
  var_label <- if(var=='R0'){'R0 (calculated after)'}else{
    if(var=='init_infected'){'Initial infected (log10)'}else{var_label}
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
        theme(legend.position = 'none') +
      labs(y = var_label)
    
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
      theme(legend.position = 'none') +
      labs(y = var_label)
    
  }
    
}

plot_trace <- function(var, filtered = F){
  
  data <- if(filtered){mcmc_samples_filtered}else{mcmc_samples} 
  
  if(var != 'likelihood'){
    data <- data %>% select(!likelihood)
  }
  
  var_label <- gsub('_rate_','_rate\n', var)
  var_label <- gsub('_spline_','_spline\n', var_label)
  
  var_label <- if(var=='R0'){'R0 (calculated after)'}else{
    if(var=='init_infected'){'Initial infected (log10)'}else{var_label}
  }
  
  use_colors <- n_distinct(data$chain) <= 10
  
  epid_pars_joining <- epid_pars %>% 
    mutate(true_value = value) %>% 
    filter(name == var) %>% 
    select(epidemic, epidemic_of_season, season, subtype, name, true_value)
  
  if(!grepl('imd_spline|likelihood', var)){
    
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
      geom_hline(aes(yintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      facet_grid(.~subtype_season) +
      theme_bw() + labs(y = var_label) +
      theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
    
  }else{
    
    if(var == 'likelihood'){
      
      epid_pars_joining <- epid_pars %>% 
        select(epidemic, season) %>% unique()
      
      p <- data %>%
        pivot_longer(!c(iteration,epidemic,chain)) %>%
        filter(name == var) %>%
        left_join(epid_pars_joining, by = c('epidemic')) %>% 
        filter(!is.na(season)) %>% # remove non-real epidemics
        mutate(epidemic = paste0("Season ", epidemic)) %>% 
        ggplot() +
        geom_line(aes(x = iteration/1000, y = value, col = as.factor(chain), group = chain)) +
        geom_hline(yintercept = 0, lty=2) +
        facet_grid(.~season) +
        theme_bw() + labs(y = var_label) +
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
        geom_hline(aes(yintercept = true_value, alpha = true_value_NA), lty=2) +
        scale_alpha_manual(values = c(1,0)) +
        facet_grid(.~season) +
        theme_bw() + labs(y = var_label) +
        theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
      
    }
    
  }
  
  if(var=='init_infected'){
    p <- p + scale_y_log10()
  }
  
  if(use_colors){
    
    p <- p +
      geom_line(aes(x = iteration/1000, y = value, col = as.factor(chain), group = chain))
      
  }else{
    
    p <- p +
      geom_line(aes(x = iteration/1000, y = value, group = chain), alpha = 0.2)
    
  }
  
  p
  
}












