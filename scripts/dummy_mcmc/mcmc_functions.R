## MCMC FUNCTIONS ##

## WITH UNKNOWN REPORTING RATES ##
run_mcmc_inference <- function(
    demography_input, 
    vaccinated_input,
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
    broad_age = c(rep('children', 3), rep('adults', 4), rep('older_adults', 2))
  )
  epidemic_dt <- as.data.table(epidemic_to_fit)
  coverage_rates$imd_quintile <- factor(coverage_rates$imd_quintile)
  
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
    
    transmissibility <- pars[1]
    susceptibility <- susc_vector(pars[2:3]) 
    init_infected_num <- 10^(pars[4])
    
    care_rate_df <- data.frame(
      broad_age = rep(unique(broad_ages$broad_age), 2),
      risk_level = rep(c('low','high'), each = 3),
      primary_care = pars[5:10],
      secondary_care = pars[11:16]
    )

    care_rate_age_df <- data.table(cross_join(
      care_rate_df,
      broad_ages))[broad_age.x == broad_age.y,]
    care_rate_age_df[, c('broad_age.x','broad_age.y') := NULL]
    
    imd_spline_pars <- data.table(
      primary = pars[17:18],
      secondary = pars[19:20]
    )
    rel_imd_rep_rates <- data.frame(imd_quintile = 1:5,
                                    rel_primary_rates = imd_spline(imd_spline_pars$primary),
                                    rel_secondary_rates = imd_spline(imd_spline_pars$secondary))
    
    reporting_rates <- cross_join(
      care_rate_age_df,
      rel_imd_rep_rates
    ) %>% 
      mutate(primary_care = primary_care*rel_primary_rates,
             secondary_care = secondary_care*rel_secondary_rates) %>% 
      select(!c(rel_primary_rates, rel_secondary_rates))
    
    reporting_dt <- as.data.table(reporting_rates)
    reporting_long <- melt(reporting_dt,
                           measure.vars = c('primary_care', 'secondary_care'),
                           variable.name = 'setting',
                           value.name = 'rate')
    
    # any out-of-bounds proposals slipping past the prior
    if(any(is.na(pars)) || 
       transmissibility < min_trans || 
       sum(susceptibility < min_susc) > 0  || sum(susceptibility > max_susc) > 0 ||
       init_infected_num < 0 || init_infected_num > min(demography_input$population) ) {
      return(-Inf)
    }
    
    pop_vaccinated <- vaccinated_input$effectively_vaccinated_population
    
    init_infected_vec <- (demography_input$population - pop_vaccinated)*init_infected_num/
      (sum(demography_input$population)-sum(pop_vaccinated))
    
    time_series <- run_model(
      pop = demography_input$population,
      I0 = init_infected_vec,
      vacc = pop_vaccinated,
      cm = cm_input,
      trans = transmissibility,
      susc = susceptibility,
      lat_per = epid_periods[1],
      inf_per = epid_periods[2]
    )
    
    ## add start date (using first of september throughout)
    start_of_epidemic <- as.Date(paste0('01-09-',year(epidemic_to_fit$week_start[1])), format = '%d-%m-%Y')
    time_series <- time_series[, date := start_of_epidemic + t] ## add date
    time_series[, t := NULL]
    
    if(any(is.na(time_series))) return(-Inf)
    
    # Take rounded value of OBSERVED infections (need an integer)
    time_series <- time_series[coverage_rates, on = c('age_grp','imd_quintile','risk_level')]
    time_series_fit <- time_series[, .(infections = floor(sum(OS_COVERAGE*infections))),
                                   by = .(date, age_grp, imd_quintile, risk_level)]
    time_series_fit[, imd_quintile := as.numeric(imd_quintile)]
    
    if(any(is.na(time_series_fit))) return(-Inf)
    
    # Check epidemic is growing at start, and not at the end
    setorder(time_series_fit, age_grp, imd_quintile, risk_level)
    if(time_series_fit$infections[1] > time_series_fit$infections[5]) return(-Inf)
    if(time_series_fit$infections[nrow(time_series_fit)] > time_series_fit$infections[nrow(time_series_fit)-5]) return(-Inf)
    
    # Aggregate to weekly and join with surveillance data
    time_series_fit[, week_start := last_monday(date)]
    time_series_weekly <- time_series_fit[, .(infections = sum(infections)), 
                                          by = .(week_start, age_grp, imd_quintile, risk_level)]
    
    # Join with observed data
    time_series_joint <- merge(time_series_weekly, epidemic_dt, 
                               by = c('age_grp', 'imd_quintile', 'week_start', 'risk_level'), 
                               all.x = F)
    
    # Pivot primary_care and secondary_care to long format
    time_series_long <- melt(time_series_joint, 
                             measure.vars = c('primary_care', 'secondary_care'),
                             variable.name = 'setting', 
                             value.name = 'observations')
    
    # Join reporting rates
    time_series_long <- merge(time_series_long, reporting_long,
                              by = colnames(reporting_long)[colnames(reporting_long) != 'rate'],
                              all.x = TRUE)
    
    # Validate reporting rates
    if(any(is.na(time_series_long$rate))) return(-Inf)
    if(any(time_series_long$rate <= 0 | time_series_long$rate >= 1)) return(-Inf)
    
    ## shift observations by the care delays
    time_series_shifted <- rbind(
      time_series_long[setting=='primary_care'][,
                                                observations := shift(observations, n = -delays['primary'], fill = 0),
                                                by = .(age_grp, imd_quintile, risk_level, index, rate)],
      time_series_long[setting=='secondary_care'][,
                                                  observations := shift(observations, n = -delays['secondary'], fill = 0),
                                                  by = .(age_grp, imd_quintile, risk_level, index, rate)])
    
    if(F){
      time_series_shifted %>% ## plot {observations} against {infections x reporting rates}
        ggplot() + theme_bw() +
        geom_line(aes(week_start,rate*infections, lty=risk_level, col=setting, group=interaction(setting, risk_level))) +
        geom_line(aes(week_start,observations, lty=risk_level, col=setting, group=interaction(setting, risk_level))) +
        facet_grid(age_grp ~ imd_quintile, scales = 'free')
      time_series_shifted %>%
        ggplot() + theme_bw() +
        geom_jitter(aes(week_start, infections >= observations, shape=risk_level, col=setting, 
                        group=interaction(setting, risk_level)), width = 1, height = 0.1) +
        facet_grid(age_grp ~ imd_quintile, scales = 'free') # should be TRUE everywhere
      }
    
    # Vectorised log likelihood
    total_ll <- sum(dbinom(
      x    = time_series_shifted$observations,
      size = time_series_shifted$infections,
      prob = time_series_shifted$rate,
      log  = TRUE
    ), na.rm = TRUE)
    
    if(is.nan(total_ll) | is.infinite(total_ll)) return(-Inf)
    
    return(total_ll)
  }
  
  llprior <- function(pars) {
    
    # pars[1] is transmissibility, pars[2:3] are susceptibility (children, adults, older adults), 
    # pars[4] is log of initial infected on 1st september
    
    if(
      pars[1] < min_trans || pars[1] > max_trans ||
      pars[2] < min_susc || pars[2] > max_susc ||
      pars[3] < min_susc || pars[3] > max_susc ||
      pars[4] < min_log_init_inf || pars[4] > max_log_init_inf ||
      any(pars[5:16] < min_reporting) || any(pars[5:16] > max_reporting) ||
      any(pars[17:20] < min_spline) || any(pars[17:20] > max_spline)
    ) {return(-Inf)}
    
    lprob <- 0
    
    R0 <- R0_func(susceptibility = susc_vector(pars[2:3]),
                  inf_period = epid_periods[2],
                  beta_in = pars[1],
                  cm_in = cm_input,
                  per_capita = T,
                  population_vector = demography_input$population)
    
    if(R0 < min_r0 || R0 > max_r0) { return(-Inf) }
    
    # Scaled Beta prior on R0
    r0_scaled <- (R0 - min_r0) / (max_r0 - min_r0)  # rescale to [0,1]
    lprob <- lprob + dbeta(r0_scaled,
                           shape1 = r0_beta['a'],
                           shape2 = r0_beta['b'], log = TRUE)
    # Uniform prior on transmissibility 
    lprob <- lprob + dunif(unname(pars[1]), min = min_trans, max = max_trans, log = TRUE)
    # Uniform prior on susceptibility, x2
    lprob <- lprob + sum(dunif(unname(pars[2:3]), min = min_susc, max = max_susc, log = TRUE))
    # Uniform prior on initial infected
    lprob <- lprob + dunif(unname(pars[4]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
    # Beta prior on primary care reporting rates (pars 5:10), centred at 2%
    lprob <- lprob + sum(dbeta(pars[5:10], 
                               shape1 = prim_beta['a'], 
                               shape2 = prim_beta['b'], log = TRUE))
    # Beta prior on secondary care reporting rates (pars 11:16), centred at 0.5%
    lprob <- lprob + sum(dbeta(pars[11:16], 
                               shape1 = sec_beta['a'], 
                               shape2 = sec_beta['b'], log = TRUE))
    # Uniform prior on spline parameters
    lprob <- lprob + sum(dunif(unname(pars[17:20]), min = min_spline, max = max_spline, log = TRUE))
    return(lprob)
  }
  
  # Helper to get Beta parameters from mean and concentration
  beta_pars <- function(mean, concentration) {
    c(a = mean * concentration, b = (1 - mean) * concentration)
  }
  
  # Primary care: centred at 0.02, secondary: centred at 0.005
  prim_beta  <- beta_pars(0.02, 200) 
  sec_beta   <- beta_pars(0.005, 200)
  # R0 centred at 2 (scaled mean = 0.5)
  r0_beta <- beta_pars(0.5, 10)  # a=5, b=5
  
  ## set bounds
  min_trans <- 0; max_trans <- 1
  min_susc <- 0; max_susc <- 5
  min_r0 <- 1; max_r0 <- 3
  min_log_init_inf <- 0; max_log_init_inf <- log10(min(demography_input$population))
  min_reporting <- 0; max_reporting <- 1
  min_spline <- -log(5); max_spline <- log(5) # equivalent to IMD ratios at most 
  # 5x higher/lower than each other in IMD 1 vs 3 or IMD 5 vs 3,
  # which would lead to huge ratios between IMD 1 and IMD 5 at the extremes
  
  ## set up sampler
  lower_vals <- c(min_trans, rep(min_susc, 2), min_log_init_inf,
                  rep(min_reporting, 12), rep(min_spline, 4))
  upper_vals <- c(max_trans, rep(max_susc, 2), max_log_init_inf,
                  rep(max_reporting, 12), rep(max_spline, 4))
  
  sampler <- function(n = 1){
    out <- matrix(NA, nrow = n, ncol = length(initial_parameters))
    for(j in 1:n){
      valid <- FALSE
      while(!valid){
        r0_scaled <- rbeta(1, shape1 = r0_beta['a'], shape2 = r0_beta['b'])
        R0 <- min_r0 + r0_scaled * (max_r0 - min_r0)
        
        susc_1 <- runif(1, min_susc, max_susc)
        susc_2 <- runif(1, min_susc, max_susc)
        
        trans <- R0_func(
          susceptibility    = susc_vector(c(susc_1, susc_2)),
          inf_period        = epid_periods[2],
          beta_in           = 1,
          cm_in             = cm_input,
          per_capita        = TRUE,
          population_vector = demography_input$population,
          R0assumed         = R0,
          return_beta       = TRUE
        )
        
        if(trans >= min_trans & trans <= max_trans) valid <- TRUE
      }
      
      # Check trans is within bounds
      if(trans < min_trans | trans > max_trans) next
      
      # Draw reporting rates from Beta priors
      prim_rates <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
      sec_rates  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
      
      out[j, ] <- c(
        trans, susc_1, susc_2,
        runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
        prim_rates,                                     # primary care rates
        sec_rates,                                      # secondary care rates
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
  data %>%
    pivot_longer(!c(iteration,epidemic,chain)) %>%
    filter(name == var) %>%
    mutate(epidemic = paste0("Epidemic ", epidemic)) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = as.factor(chain), group = chain), alpha = 0.4) +
    geom_vline(data = epid_pars %>% mutate(epidemic = paste0("Epidemic ", epidemic)) %>% filter(name==var),
               aes(xintercept = value), lty=2) +
    theme_bw() + labs(y = ifelse(var=='R0', 'R0 (calculated after)', var)) +
    facet_grid(.~epidemic) +
    # scale_fill_manual(values = var_cols) +
    theme(legend.position = 'none')
}

plot_trace <- function(var, filtered = F){
  data <- if(filtered){mcmc_samples_filtered}else{mcmc_samples}
  Y_LAB <- if(var=='R0'){'R0 (calculated after)'}else{
    if(var=='init_infected'){'Initial infected (log10)'}else{var}
  }
  p <- data %>%
    mutate(epidemic = paste0("Epidemic ", epidemic)) %>% 
    pivot_longer(!c(iteration,epidemic,chain)) %>%
    filter(name == var) %>%
    ggplot() +
    geom_line(aes(x = iteration, y = value, col = as.factor(chain), group = chain)) +
    geom_hline(data = epid_pars %>% mutate(epidemic = paste0("Epidemic ", epidemic)) %>% filter(name==var),
               aes(yintercept = value), lty=2) +
    # geom_vline(xintercept = burn_in, lty=3, alpha = 0.5) +
    facet_grid(.~epidemic) +
    theme_bw() + labs(y = Y_LAB) +
    # scale_color_manual(values = var_cols) +
    theme(legend.position = 'none')
  
  if(var=='init_infected'){
    p <- p + scale_y_log10()
  }
  
  p
}


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



