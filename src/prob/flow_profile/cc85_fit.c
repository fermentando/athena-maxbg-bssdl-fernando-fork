#include "cc85_fit.h"

#include <math.h>

#define VELOCITY_CONV 1.1192359020509211e-03 // from km/s to code units

// Initial offset (in units of R*)
static Real r0 = 1.7;

// For velocity profile
static Real vinf     = 1278.01660112; // km/s
static Real lgvinf = log(vinf * VELOCITY_CONV);
static Real v_params[] = {1.52693591, 0.92286245}; // alpha, x0

// Density & pressure
static Real d_params[] = {exp(6.9118067), -2.19446264,   0.9560074}; // y0, m, alpha
static Real P_params[] = {exp(5.96088173), -3.65743773,  0.9560074}; // y0, m, alpha


void flow_profile_init() {
  return;
}

Real flow_profile_velocity(Real r) {
  return exp(lgvinf * (1 - v_params[1] * exp(-v_params[0] * log(r + r0))));
}


Real flow_profile_density(Real r) {
  return d_params[0] * exp(d_params[1] * pow(log(r + r0), d_params[2]))
}


Real flow_profile_pressure(Real r) {
  return P_params[0] * exp(P_params[1] * pow(log(r + r0), P_params[2]))
}



#undef VELOCITY_CONV
