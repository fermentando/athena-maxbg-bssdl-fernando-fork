// Original fit from Mike to the Sutherland & Dopita (1993), Z~1 cooling curve
// Temperatures in keV
// Lambda in 1e-23 erg/s

#define nfit_cool 7

static Real sdT[nfit_cool] = {
  1.0e-5,
  0.0017235,
  0.02,
  0.13,
  0.7,
  5.0,
  100.0};

static const Real sdL[nfit_cool] = {
  5.890872e-13,
  15.438249,
  66.831473,
  2.773501,
  1.195229,
  1.842056,
  6.10541
};

static const Real sdexpt[nfit_cool] = {
  6.0,
  0.6,
  -1.7,
  -0.5,
  0.22,
  0.4,
  0.4
};
