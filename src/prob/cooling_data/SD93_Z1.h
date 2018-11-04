/*
  Original fit from Mike to the Sutherland & Dopita (1993), Z~1 cooling curve
  Temperatures in keV
  Lambda in 1e-23 erg/s
*/

#define nfit_cool_d 1
#define nfit_cool_T 7

// Temperature bins
static Real sdT[nfit_cool_T] = {
  1.0e-5,
  0.0017235,
  0.02,
  0.13,
  0.7,
  5.0,
  100.0};

#if nfit_cool_d == 1

// Density bins
static Real sdd[nfit_cool_d] = {
  1e-24
};

// Lambda_k
static const Real sdL[nfit_cool_d][nfit_cool_T] = {{
  5.890872e-13,
  15.438249,
  66.831473,
  2.773501,
  1.195229,
  1.842056,
  6.10541
}};

// alpha_k
static const Real sdexpt[nfit_cool_d][nfit_cool_T] = {{
  6.0,
  0.6,
  -1.7,
  -0.5,
  0.22,
  0.4,
  0.4
}};


// for testing density dependent cooling function
#elif nfit_cool_d == 2
static Real sdd[nfit_cool_d] = {
  1e-24,
  5e4
};
static const Real sdL[nfit_cool_d][nfit_cool_T] = {
  {5.890872e-13, 15.438249, 66.831473, 2.773501, 1.195229, 1.842056, 6.10541},
  {5.890872e-13, 15.438249, 66.831473, 2.773501, 1.195229, 1.842056, 6.10541}
};
static const Real sdexpt[nfit_cool_d][nfit_cool_T] = {
  {6.0, 0.6, -1.7, -0.5, 0.22, 0.4, 0.4},
  {6.0, 0.6, -1.7, -0.5, 0.22, 0.4, 0.4}
};
#else
#error "Invalid `nfit_cool_d`."
#endif
