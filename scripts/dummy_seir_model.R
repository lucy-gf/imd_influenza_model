## SEIR MODEL FOR DUMMY DATA PRODUCTION ##

suppressPackageStartupMessages(library(odin))

#### SEIR MODEL ####

seeiir_risk_odin <- odin::odin({
  
  # User supplied parameters
  
  no_groups <- user()          # 80 demographic groups
  pop[] <- user()              # length = 2 * no_groups
  I0[] <- user()               # initial infectious (distributed into I1)
  vacc[] <- user()             # vaccinated at t0 (enter R)
  
  trans <- user()
  lat_per <- user()
  inf_per <- user()
  
  cij[,] <- user()             # (2*no_groups x 2*no_groups)
  
  # Derived quantities
  
  lat_rate <- 2 / lat_per      # two E stages
  inf_rate <- 2 / inf_per      # two I stages
  
  I_tot[] <- I1[i] + I2[i]
  
  sij[,] <- cij[i,j] * I_tot[j] 
  lambda[] <- trans * sum(sij[i,])
  
  newInf[] <- lambda[i] * S[i]
  
  progE1[] <- lat_rate * E1[i]
  progE2[] <- lat_rate * E2[i]
  progI1[] <- inf_rate * I1[i]
  progI2[] <- inf_rate * I2[i]
  
  # ODEs
  
  deriv(S[])  <- - newInf[i]
  
  deriv(E1[]) <- newInf[i] - progE1[i]
  deriv(E2[]) <- progE1[i] - progE2[i]
  
  deriv(I1[]) <- progE2[i] - progI1[i]
  deriv(I2[]) <- progI1[i] - progI2[i]
  
  deriv(R[])  <- progI2[i]
  
  deriv(cumI[]) <- newInf[i]
  
  # Initial conditions
  
  initial(S[])  <- pop[i] - I0[i] - vacc[i]
  
  initial(E1[]) <- 0
  initial(E2[]) <- 0
  
  initial(I1[]) <- I0[i]
  initial(I2[]) <- 0
  
  initial(R[])  <- vacc[i]
  
  initial(cumI[]) <- 0
  
  # Dimensions
  
  dim(pop) <- no_groups
  dim(I0) <- no_groups
  dim(vacc) <- no_groups
  
  dim(lambda) <- no_groups
  dim(newInf) <- no_groups
  
  dim(progE1) <- no_groups
  dim(progE2) <- no_groups
  dim(progI1) <- no_groups
  dim(progI2) <- no_groups
  
  dim(I_tot) <- no_groups
  
  dim(S) <- no_groups
  dim(E1) <- no_groups
  dim(E2) <- no_groups
  dim(I1) <- no_groups
  dim(I2) <- no_groups
  dim(R) <- no_groups
  dim(cumI) <- no_groups
  
  dim(cij) <- c(no_groups, no_groups)
  dim(sij) <- c(no_groups, no_groups)
  
})

# seir model
seir_odin <- odin::odin({
  
  no_groups <- user()
  
  pop[] <- user()
  I0[] <- user()
  
  trans <- user()
  lat_per <- user()
  inf_per <- user()
  
  cij[,] <- user()
  
  sij[,] <- cij[i,j]*I[j]/sum(pop[]) # transmission matrix
  lambda[] <- trans*sum(sij[i,]) # FOI
  newInf[] <- lambda[i] * S[i] # Newly infected
  onsets[] <- E[i]/lat_per
  removal[] <- I[i]/inf_per
  
  # ODEs
  deriv(S[]) <- - newInf[i]
  deriv(E[]) <- newInf[i] - onsets[i]
  deriv(I[]) <- onsets[i] - removal[i]
  deriv(R[]) <- removal[i]
  deriv(cumI[]) <- newInf[i]
  
  # Initial values
  initial(S[]) <- pop[i] - I0[i]
  initial(E[]) <- 0
  initial(I[]) <- I0[i]
  initial(R[]) <- 0
  initial(cumI[]) <- 0
  
  # Dimensions
  dim(pop) <- no_groups
  dim(I0) <- no_groups
  dim(lambda) <- no_groups
  dim(newInf) <- no_groups
  dim(onsets) <- no_groups
  dim(removal) <- no_groups
  dim(S) <- no_groups
  dim(E) <- no_groups
  dim(I) <- no_groups
  dim(R) <- no_groups
  dim(cumI) <- no_groups
  
  dim(cij) <- c(no_groups, no_groups)
  dim(sij) <- c(no_groups, no_groups)
  
})

#### FORMATTING EPI DATA FUNCTIONS ####

make_group_lookup <- function() {
  
  base_df <- expand.grid(
    age_grp = age_labels,
    imd_quintile = paste0(1:nimd)
  ) %>%
    arrange(imd_quintile, age_grp)
  
  # 80 demographic groups
  base_df$demo_index <- 1:nrow(base_df)
  
  # Duplicate for risk groups
  low  <- base_df %>% mutate(risk_level = "low",  group = demo_index)
  high <- base_df %>% mutate(risk_level = "high", group = demo_index + ng)
  
  bind_rows(low, high) %>%
    arrange(group) %>% 
    mutate(group=as.character(group))
}

run_model <- function(
    pop,
    I0,
    vacc,
    cm,
    trans,
    lat_per,
    inf_per,
    t_end = 365
) {
  
  model <- seeiir_risk_odin$new(
    no_groups = ndim,
    pop = pop,
    I0 = I0,
    vacc = vacc,
    trans = trans,
    lat_per = lat_per,
    inf_per = inf_per,
    cij = cm
  )
  
  times <- seq(0, t_end, by = 0.1)
  
  out <- model$run(times)
  
  out_formatted <- tidy_output(data.table(out))
  
  return(out_formatted)
}

tidy_output <- function(out_dt) {
  
  lookup <- make_group_lookup()
  
  long_dt <- melt.data.table(out_dt, id.vars='t')
  
  long_dt[, compartment := gsub("(.+)\\[(\\d+)\\]", "\\1", long_dt$variable)]
  long_dt[, group := gsub("(.+)\\[(\\d+)\\]", "\\2", long_dt$variable)]
  long_dt[, variable := NULL]
  
  long_dt[compartment %in% c('E1','E2'), compartment := "E"]
  long_dt[compartment %in% c('I1','I2'), compartment := "I"]
  
  long_dt <- long_dt[, lapply(.SD, sum), by = c('t','compartment','group')]
  
  long_dt <- long_dt[lookup, on = "group"]
  
  return(long_dt[, c('t','age_grp','imd_quintile','risk_level','compartment','value')])
}

#### R0 ####

R0_func <- function(susceptibility,
                    inf_period,
                    beta_in,
                    cm_in,
                    R0assumed = NULL,
                    return_beta = F){
  
  ng    = dim(cm_in)[1]
  ngm   = cm_in
  
  for (k in 1:ng){ 
    for (j in 1:ng){
      ngm[j,k] = beta_in*susceptibility*cm_in[j,k]*inf_period 
    }
  }
  
  # max EV
  EVs = eigen(ngm)$values
  R00 = max(Re(EVs[which(Im(EVs)==0)]))
  
  if(return_beta){
    return(R0assumed/(R00/beta_in))  
  }else{
    return(R00)
  }
  
}






