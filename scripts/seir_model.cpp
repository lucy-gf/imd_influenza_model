// seir_model.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix run_seir_cpp(
    NumericVector pop,
    NumericVector I0,
    NumericVector vacc,
    double trans,
    NumericVector susc,
    double lat_per,
    double inf_per,
    NumericMatrix cij,
    double t_end = 250.0,
    double dt = 0.1
) {
  int ng = pop.size();
  
  // Rates — identical to odin (two E and two I stages)
  double lat_rate = 2.0 / lat_per;
  double inf_rate = 2.0 / inf_per;
  
  // Initialise compartments
  NumericVector S(ng), E1(ng), E2(ng), I1(ng), I2(ng), R(ng), V(ng), cumI(ng);
  for (int i = 0; i < ng; i++) {
    S[i]    = pop[i] - I0[i] - vacc[i];
    E1[i]   = 0.0;
    E2[i]   = 0.0;
    I1[i]   = I0[i];
    I2[i]   = 0.0;
    R[i]    = 0.0;
    V[i]    = vacc[i];
    cumI[i] = 0.0;
  }
  
  // Output: one row per integer time point, columns = t + ng cumI values
  int n_steps     = (int)(t_end / dt) + 1;
  int n_out       = (int)(t_end) + 1;
  NumericMatrix out(n_out, ng + 1);
  int out_row = 0;
  
  // Store t=0
  out(out_row, 0) = 0.0;
  for (int i = 0; i < ng; i++) out(out_row, i + 1) = cumI[i];
  out_row++;
  
  // Temporary derivative vectors
  NumericVector I_tot(ng), lambda(ng), newInf(ng);
  NumericVector progE1(ng), progE2(ng), progI1(ng), progI2(ng);
  
  double t = 0.0;
  for (int step = 0; step < n_steps - 1; step++) {
    
    // I_tot[i] = I1[i] + I2[i]
    for (int i = 0; i < ng; i++) I_tot[i] = I1[i] + I2[i];
    
    // lambda[i] = susc[i] * trans * sum_j(cij[i,j] * I_tot[j])
    for (int i = 0; i < ng; i++) {
      double s = 0.0;
      for (int j = 0; j < ng; j++) s += cij(i, j) * I_tot[j];
      lambda[i] = susc[i] * trans * s;
    }
    
    // Flows
    for (int i = 0; i < ng; i++) {
      newInf[i] = lambda[i] * S[i];
      progE1[i] = lat_rate * E1[i];
      progE2[i] = lat_rate * E2[i];
      progI1[i] = inf_rate * I1[i];
      progI2[i] = inf_rate * I2[i];
    }
    
    // Euler update — identical to odin deriv() equations
    for (int i = 0; i < ng; i++) {
      S[i]    += dt * (-newInf[i]);
      E1[i]   += dt * (newInf[i]  - progE1[i]);
      E2[i]   += dt * (progE1[i]  - progE2[i]);
      I1[i]   += dt * (progE2[i]  - progI1[i]);
      I2[i]   += dt * (progI1[i]  - progI2[i]);
      R[i]    += dt * (progI2[i]);
      cumI[i] += dt * (newInf[i]);
      // V is constant (deriv = 0), no update needed
    }
    
    t += dt;
    
    // Save at integer time points
    double t_next = t + dt;
    if (std::fmod(t + 1e-9, 1.0) < dt && out_row < n_out) {
      out(out_row, 0) = std::round(t);
      for (int i = 0; i < ng; i++) out(out_row, i + 1) = cumI[i];
      out_row++;
    }
  }
  
  return out;
}