## SEIR MODEL FOR DUMMY DATA PRODUCTION ##

# suppressPackageStartupMessages(library(odin))
Rcpp::sourceCpp('scripts/seir_model.cpp')

#### RUN SEIR MODEL ####

run_model <- function(
    pop,
    I0,
    vacc,
    cm,
    trans,
    susc,
    lat_per,
    inf_per,
    t_end = 250
) {
  
  # susceptibility is age-speciifc
  if(length(susc) != dim(cm)[1]){ 
    susc <- rep(susc, dim(cm)[1]/length(susc))
  }
  
  ## OLD CODE FROM ODIN SEIR MODEL
  # model <- seeiir_risk_odin$new(
  #   no_groups = ndim,
  #   pop = pop,
  #   I0 = I0,
  #   vacc = vacc,
  #   trans = trans,
  #   susc = susc,
  #   lat_per = lat_per,
  #   inf_per = inf_per,
  #   cij = cm
  # )
  # 
  # times <- seq(0, t_end, by = 0.1)
  # 
  # out <- data.table(model$run(times, method = "euler"))
  ## using euler solver for speed in the MCMC
  # # keep only cumI columns
  # cumI_cols <- c('t', grep('^cumI', colnames(out), value = TRUE))
  # out <- out[, ..cumI_cols]
  # 
  # out <- out[t %% 1 == 0] # keep only integer time points 
  
  # call C++ solver
  raw <- run_seir_cpp(
    pop    = pop,
    I0     = I0,
    vacc   = vacc,
    trans  = trans,
    susc   = susc,
    lat_per = lat_per,
    inf_per = inf_per,
    cij    = cm,
    t_end  = t_end,
    dt     = 0.1
  )
  
  # raw is matrix: col 1 = t, cols 2:(ng+1) = cumI per group
  ng  <- ncol(raw) - 1
  out <- data.table(raw)
  setnames(out, c('t', paste0('cumI[', 1:ng, ']')))
  
  out_formatted <- tidy_output(out)
  
  ## infections = cumulative(t) - cumulative(t-1)
  setorder(out_formatted, age_grp, imd_quintile, risk_level, t)
  out_formatted[, infections := value - shift(value, fill = 0), 
                   by = .(age_grp, imd_quintile, risk_level)]
  
  out_formatted[, c('compartment','value') := NULL]
  
  return(out_formatted)
}

#### FORMATTING EPI DATA FUNCTIONS ####

make_group_lookup <- function() {
  
  base_df <- expand.grid(
    age_grp = age_labels,
    imd_quintile = paste0(1:nimd)
  ) %>%
    arrange(imd_quintile, age_grp)
  
  # demographic groups
  base_df$demo_index <- 1:nrow(base_df)
  
  # Duplicate for risk groups
  low  <- base_df %>% mutate(risk_level = "low",  group = demo_index)
  high <- base_df %>% mutate(risk_level = "high", group = demo_index + nrow(base_df))
  
  bind_rows(low, high) %>%
    arrange(group) %>% 
    mutate(group=as.character(group))
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

#### EXPAND CONTACT MATRIX ####

expand_contact_matrix <- function(c45) {
  
  ndim <- nrow(c45)
  
  stopifnot(nrow(c45) == ncol(c45))
  
  # Create empty 90x90 matrix
  c90 <- matrix(0, nrow = 2*ndim, ncol = 2*ndim)
  
  # Fill all 4 quadrants with same matrix
  c90[1:ndim, 1:ndim]  <- c45   # low → low
  c90[1:ndim, (ndim+1):(2*ndim)]   <- c45   # low → high
  c90[(ndim+1):(2*ndim), 1:ndim]   <- c45   # high → low
  c90[(ndim+1):(2*ndim), (ndim+1):(2*ndim)] <- c45   # high → high
  
  return(c90)
}

#### R0 ####

R0_func <- function(susceptibility,
                    inf_period,
                    beta_in,
                    cm_in,
                    population_vector = NULL,
                    per_capita = F,
                    R0assumed = NULL,
                    return_beta = F){
  
  ng    = dim(cm_in)[1]
  ngm   = cm_in
  
  if(!is.null(colnames(cm_in))){
    if(substr(colnames(cm_in)[1],1,1) != substr(colnames(cm_in)[2],1,1)){
      stop('Contact matrix not arranged by IMD first!')
    }
  }
  
  if(per_capita){
    if(is.null(population)){stop('Need population sizes')}
    cm_in <- t(t(cm_in)*population_vector) 
  }
  
  if(length(susceptibility) > 1){ # age-specific susceptibility
    susceptibility <- rep(susceptibility, ng/length(susceptibility))
    for (k in 1:ng){ 
      for (j in 1:ng){
        ngm[j,k] = beta_in*susceptibility[k]*cm_in[j,k]*inf_period 
      }
    }  
  }else{
    for (k in 1:ng){ 
      for (j in 1:ng){
        ngm[j,k] = beta_in*susceptibility*cm_in[j,k]*inf_period 
      }
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

#### FIND LAST MONDAY ####

last_monday <- function(dates){
  
  dates <- as.Date(dates)
  as.Date(dates - ((as.integer(format(dates, "%u")) - 1) %% 7))
  
}


# #### OLD SEIR MODEL (ODIN) ####
# 
# seeiir_risk_odin <- odin::odin({
#   
#   # User supplied parameters
#   
#   no_groups <- user()          # demographic groups
#   pop[] <- user()              
#   I0[] <- user()               # initial infectious (distributed into I1)
#   vacc[] <- user()             # vaccinated at t0 (enter R)
#   
#   trans <- user()
#   susc[] <- user()               # length = 2 * no_groups
#   lat_per <- user()
#   inf_per <- user()
#   
#   cij[,] <- user()             # (2*no_groups x 2*no_groups)
#   
#   # Derived quantities
#   
#   lat_rate <- 2 / lat_per      # two E stages
#   inf_rate <- 2 / inf_per      # two I stages
#   
#   I_tot[] <- I1[i] + I2[i]
#   
#   sij[,] <- cij[i,j] * I_tot[j] 
#   lambda[] <- susc[i] * trans * sum(sij[i,])
#   
#   newInf[] <- lambda[i] * S[i]
#   
#   progE1[] <- lat_rate * E1[i]
#   progE2[] <- lat_rate * E2[i]
#   progI1[] <- inf_rate * I1[i]
#   progI2[] <- inf_rate * I2[i]
#   
#   # ODEs
#   
#   deriv(S[])  <- - newInf[i]
#   
#   deriv(E1[]) <- newInf[i] - progE1[i]
#   deriv(E2[]) <- progE1[i] - progE2[i]
#   
#   deriv(I1[]) <- progE2[i] - progI1[i]
#   deriv(I2[]) <- progI1[i] - progI2[i]
#   
#   deriv(R[])  <- progI2[i]
#   
#   deriv(V[])  <- 0
#   
#   deriv(cumI[]) <- newInf[i]
#   
#   # Initial conditions
#   
#   initial(S[])  <- pop[i] - I0[i] - vacc[i]
#   
#   initial(E1[]) <- 0 
#   initial(E2[]) <- 0
#   
#   initial(I1[]) <- I0[i]
#   initial(I2[]) <- 0
#   
#   initial(R[])  <- 0
#   initial(V[])  <- vacc[i]
#   
#   initial(cumI[]) <- 0
#   
#   # Dimensions
#   
#   dim(pop) <- no_groups
#   dim(I0) <- no_groups
#   dim(vacc) <- no_groups
#   dim(susc) <- no_groups
#   dim(lambda) <- no_groups
#   dim(newInf) <- no_groups
#   
#   dim(progE1) <- no_groups
#   dim(progE2) <- no_groups
#   dim(progI1) <- no_groups
#   dim(progI2) <- no_groups
#   
#   dim(I_tot) <- no_groups
#   
#   dim(S) <- no_groups
#   dim(E1) <- no_groups
#   dim(E2) <- no_groups
#   dim(I1) <- no_groups
#   dim(I2) <- no_groups
#   dim(R) <- no_groups
#   dim(V) <- no_groups
#   dim(cumI) <- no_groups
#   
#   dim(cij) <- c(no_groups, no_groups)
#   dim(sij) <- c(no_groups, no_groups)
#   
# })




