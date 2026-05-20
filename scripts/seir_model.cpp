// seir_model.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix run_seir_cpp(
    NumericVector pop,
    NumericVector I0,
    NumericVector vacc_cov,
    NumericVector ve_inf,
    double trans,
    NumericVector susc,
    double lat_per,
    double inf_per,
    NumericMatrix cij,
    double t_end = 250.0,
    double dt = 0.1
) {
  int ng = pop.size();
  
  // Rates (two E and two I stages)
  double lat_rate = 2.0 / lat_per;
  double inf_rate = 2.0 / inf_per;
  
  // Initialise compartments
  NumericVector S(ng), E1(ng), E2(ng), I1(ng), I2(ng), R(ng), V(ng), cumI(ng);
  NumericVector Sv(ng), E1v(ng), E2v(ng), I1v(ng), I2v(ng), Rv(ng), cumIv(ng);
  for (int i = 0; i < ng; i++) {
    S[i]    = pop[i] - vacc_cov[i] - I0[i];
    E1[i]   = 0.0;
    E2[i]   = 0.0;
    I1[i]   = I0[i];
    I2[i]   = 0.0;
    R[i]    = 0.0;
    Sv[i]    = (1 - ve_inf[i])*vacc_cov[i];
    E1v[i]   = 0.0;
    E2v[i]   = 0.0;
    I1v[i]   = 0.0;
    I2v[i]   = 0.0;
    Rv[i]    = 0.0;
    V[i]    = ve_inf[i]*vacc_cov[i]; // those effectively vaccinated
    cumI[i] = 0.0;
    cumIv[i] = 0.0;
  }
  
  // Output: one row per integer time point, columns = t + ng cumI values
  int n_steps = (int)(t_end / dt) + 1;
  int n_out   = (int)(t_end) + 1;
  int steps_per_unit = (int)std::round(1.0 / dt); // robust integer-time check
  NumericMatrix out(n_out, 2 * ng + 1);            
  int out_row = 0;
  
  // Store t=0
  out(out_row, 0) = 0.0;
  for (int i = 0; i < ng; i++) {
    out(out_row, i + 1)      = cumI[i];
    out(out_row, i + 1 + ng) = cumIv[i];           
  }
  out_row++;
  
  // Temporary derivative vectors
  NumericVector I_tot(ng), lambda(ng), newInf(ng), newInfv(ng);
  NumericVector progE1(ng), progE2(ng), progI1(ng), progI2(ng);
  NumericVector progE1v(ng), progE2v(ng), progI1v(ng), progI2v(ng);
  
  double t = 0.0;
  for (int step = 0; step < n_steps - 1; step++) {
    
    // I_tot[i] = I1[i] + I2[i] + I1v[i] + I2v[i]
    for (int i = 0; i < ng; i++) I_tot[i] = I1[i] + I2[i] + I1v[i] + I2v[i];
    
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
      newInfv[i] = lambda[i] * Sv[i];
      progE1v[i] = lat_rate * E1v[i];
      progE2v[i] = lat_rate * E2v[i];
      progI1v[i] = inf_rate * I1v[i];
      progI2v[i] = inf_rate * I2v[i];
    }
    
    // Euler update (identical to odin deriv() equations)
    for (int i = 0; i < ng; i++) {
      S[i]    += dt * (-newInf[i]);
      E1[i]   += dt * (newInf[i]  - progE1[i]);
      E2[i]   += dt * (progE1[i]  - progE2[i]);
      I1[i]   += dt * (progE2[i]  - progI1[i]);
      I2[i]   += dt * (progI1[i]  - progI2[i]);
      R[i]    += dt * (progI2[i]);
      Sv[i]    += dt * (-newInfv[i]);
      E1v[i]   += dt * (newInfv[i]  - progE1v[i]);
      E2v[i]   += dt * (progE1v[i]  - progE2v[i]);
      I1v[i]   += dt * (progE2v[i]  - progI1v[i]);
      I2v[i]   += dt * (progI1v[i]  - progI2v[i]);
      Rv[i]    += dt * (progI2v[i]);
      cumI[i] += dt * (newInf[i]);
      cumIv[i] += dt * (newInfv[i]);
      // V is constant (deriv = 0), no update needed
    }
    
    t += dt;
    
    // Save at integer time points 
    if ((step + 1) % steps_per_unit == 0 && out_row < n_out) {
      out(out_row, 0) = std::round(t);
      for (int i = 0; i < ng; i++) {
        out(out_row, i + 1)      = cumI[i];
        out(out_row, i + 1 + ng) = cumIv[i];       // NEW
      }
      out_row++;
    }
  }
  
  return out;
}