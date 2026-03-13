## MCMC FUNCTIONS ##

run_mcmc_inference <- function(
    demography_input, 
    vaccinated_input,
    cm_input, 
    epidemic_to_fit, 
    epid_periods,
    reporting_rates,
    coverage_rates,
    care_delays,
    initial_parameters,
    n_samples, 
    nburn, 
    thinning,
    n_chains
) {
  
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
    susceptibility <- susc_vector(pars[2:4]) 
    start_day <- floor(pars[5]) # start day must be an integer
    
    # any out-of-bounds proposals slipping past the prior
    if(any(is.na(pars)) || 
       transmissibility < 0.05 || transmissibility > 0.3 ||
       sum(susceptibility < 0.1) > 0  || sum(susceptibility > 0.9) > 0  ||
       start_day < min_start_day || start_day > max_start_day) {
      return(-Inf)
    }
    
    pop_vaccinated <- vaccinated_input$effectively_vaccinated_population
    init_infected_num <- known_pars$init_infected
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
    
    ## add start date
    start_of_year <- as.Date(paste0('01-01-',year(epidemic_to_fit$week_start[1])), format = '%d-%m-%Y')
    time_series <- time_series[, date := start_of_year + (start_day - 1) + t] ## add date
    time_series[, t := NULL]
    
    if(any(is.na(time_series))) return(-Inf)
    
    # Take rounded value of OBSERVED infections (need an integer)
    coverage_rates$imd_quintile <- factor(coverage_rates$imd_quintile)
    time_series <- time_series[coverage_rates, on = c('age_grp','imd_quintile','risk_level')]
    time_series_fit <- time_series[, .(infections = round(sum(OS_COVERAGE*infections))),
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
    
    # Convert epidemic_to_fit to data.table if not already
    epidemic_dt <- as.data.table(epidemic_to_fit)
    
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
    reporting_dt <- as.data.table(reporting_rates)
    reporting_long <- melt(reporting_dt,
                           measure.vars = c('primary_care', 'secondary_care'),
                           variable.name = 'setting',
                           value.name = 'rate')
    
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

    time_series_shifted %>% ## plot {observations} against {infections x reporting rates}
      ggplot() + theme_bw() +
      geom_line(aes(week_start,rate*infections, lty=risk_level, group=interaction(setting, risk_level))) +
      geom_line(aes(week_start,observations, lty=risk_level, group=interaction(setting, risk_level)),col=2) +
      facet_grid(age_grp ~ imd_quintile, scales = 'free')

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
    
    # pars[1] is transmissibility, pars[2:4] are susceptibility (children, adults, older adults), 
    # pars[3] is start date
    
    if(
      pars[1] < 0.05 || pars[1] > 0.4 || 
      pars[2] < 0.1 || pars[2] > 0.9 ||
      pars[3] < 0.1 || pars[3] > 0.9 ||
      pars[4] < 0.1 || pars[4] > 0.9 ||
      pars[5] < min_start_day || pars[5] > max_start_day 
    ) {return(-Inf)}
    
    lprob <- 0
    
    R0 <- R0_func(susceptibility = c(rep(pars[2], 3), rep(pars[3], 4), rep(pars[4], 2)),
                  inf_period = epid_periods[2],
                  beta_in = pars[1],
                  cm_in = cm_input,
                  per_capita = T,
                  population_vector = demography_input$population)
    
    if(R0 < 1 || R0 > 3) { return(-Inf) }
    
    # Uniform prior on R0 between 1 and 2
    lprob <- lprob + dunif(R0, min = 1, max = 3, log = TRUE)
    # Uniform prior on transmissibility
    lprob <- lprob + dunif(unname(pars[1]), min = 0.05, max = 0.4, log = TRUE)
    # Uniform prior on susceptibility, x3
    lprob <- lprob + dunif(unname(pars[2]), min = 0.1, max = 0.9, log = TRUE)
    lprob <- lprob + dunif(unname(pars[3]), min = 0.1, max = 0.9, log = TRUE)
    lprob <- lprob + dunif(unname(pars[4]), min = 0.1, max = 0.9, log = TRUE)
    
    return(lprob)
  }
  
  # first week with at least 200 gp visits
  week_first_obs <- (epidemic_to_fit %>% 
    group_by(week_start) %>% summarise(prim = sum(primary_care)) %>% 
    filter(prim>20))$week_start[1]
  # give a two week cushion so that the MLE isn't on the prior boundary
  max_start_day <<- 14 + as.numeric(week_first_obs - as.Date(paste0('01-01-',year(week_first_obs)), format='%d-%m-%Y'))
  # go back up to 6 weeks 
  min_start_day <<- max(c(0, max_start_day - 7*8)) 
  
  ## set up sampler
  lower_vals <- c(0.05, rep(0.1, 3), min_start_day)
  upper_vals <- c(0.4, rep(0.9, 3), max_start_day)
  
  sampler <- function(n = 1){
    out <- matrix(NA, nrow = n, ncol = length(initial_parameters))
    # for(i in 1:length(initial_parameters)){
    #   out[, i] <- runif(n, lower_vals[i], upper_vals[i])
    # }
    for(j in 1:n){
      R0 <- 0
      while(R0 < 1 | R0 > 2){
        trans <- runif(1, lower_vals[1], upper_vals[1])
        susc_1  <- runif(1, lower_vals[2], upper_vals[2])
        susc_2  <- runif(1, lower_vals[3], upper_vals[3])
        susc_3  <- runif(1, lower_vals[4], upper_vals[4])
        R0 <- R0_func(
          susceptibility    = c(rep(susc_1, 3), rep(susc_2, 4), rep(susc_3, 2)),
          inf_period        = epid_periods[2],
          beta_in           = trans,
          cm_in             = cm_input,
          per_capita        = TRUE,
          population_vector = demography_input$population
        )
      }
      out[j, ] <- c(trans, susc_1, susc_2, susc_3, runif(1, lower_vals[5], upper_vals[5]))
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
    iterations = nburn + n_samples*thinning, 
    burnin = 0,
    thin = 1,
    message = T, nrChains=n_chains, parallel = F
  ) ## setup to save all (pre-thinning etc.)
  
  # cat(llprior(initial_parameters))
  
  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = 'DEzs', settings = settings)
  
  #plot(out); summary(out)
  
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
  data %>%
    mutate(epidemic = paste0("Epidemic ", epidemic)) %>% 
    pivot_longer(!c(iteration,epidemic,chain)) %>%
    filter(name == var) %>%
    ggplot() +
    geom_line(aes(x = iteration, y = value, col = as.factor(chain), group = chain)) +
    geom_hline(data = epid_pars %>% mutate(epidemic = paste0("Epidemic ", epidemic)) %>% filter(name==var),
               aes(yintercept = value), lty=2) +
    facet_grid(.~epidemic) +
    theme_bw() + labs(y = ifelse(var=='R0', 'R0 (calculated after)', var)) +
    # scale_color_manual(values = var_cols) +
    theme(legend.position = 'none')
}

