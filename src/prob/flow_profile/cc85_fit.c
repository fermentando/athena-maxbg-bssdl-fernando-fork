#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../defs.h"
#include "../../athena.h"
#include "../../globals.h"
#include "../../prototypes.h"

#include "cc85_fit.h"


#define VELOCITY_CONV 1.1192359020509211e-03 // from km/s to code units

// Initial offset (in units of R*)
static Real r0 = 2.0;

// For velocity profile
static Real vinf     = 1278.01660112; // km/s
static Real lgvinf;
static Real v_params[] = {1.53539389, 0.54344914};

// Density & pressure
static Real d_params[] = { 6.89385868, -2.17459827,  0.96205687};
static Real P_params[] = { 5.93096836, -3.62433045,  0.96205687};


void flow_profile_init() {
  //lgvinf = log(vinf * VELOCITY_CONV);
  lgvinf = log(1.846641123303628);
  return;
}

Real flow_profile_velocity_norm(Real r) {
  return exp(lgvinf * (1 - v_params[1] * exp(-v_params[0] * log(r))));
}

Real flow_profile_r(Real x1, Real x2, Real x3) {
  return sqrt((x1 + r0) * (x1 + r0) + x2 * x2 + x3 * x3);
}

Real flow_profile_velocity_x(Real x1, Real x2, Real x3) {
  Real r = flow_profile_r(x1, x2, x3);
  return (x1 + r0) / r * flow_profile_velocity_norm(r);
}

Real flow_profile_velocity_y(Real x1, Real x2, Real x3) { 
  Real r = flow_profile_r(x1, x2, x3);
  return x2 / r * flow_profile_velocity_norm(r);
}

Real flow_profile_velocity_z(Real x1, Real x2, Real x3) {
  Real r = flow_profile_r(x1, x2, x3);
  return x3 / r * flow_profile_velocity_norm(r);
}


Real flow_profile_density(Real x1, Real x2, Real x3) {
  Real r = flow_profile_r(x1, x2, x3);
  return exp(d_params[0]) * exp(d_params[1] * pow(log(r), d_params[2]));
}


Real flow_profile_pressure(Real x1, Real x2, Real x3) {
  Real r = flow_profile_r(x1, x2, x3);
  return exp(P_params[0]) * exp(P_params[1] * pow(log(r), P_params[2]));
}



#undef VELOCITY_CONV
