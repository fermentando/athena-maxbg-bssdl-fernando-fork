// Fit from Mike to the Sutherland & Dopita (1993), Z~1 cooling curve
// with all temperatures scaled by * 10
// Temperatures in keV
// Lambda in 1e-23 erg/s

#define nfit_cool 7

static Real sdT[nfit_cool] = {
  1.0000e-04, 1.7235e-02, 2.0000e-01, 1.3000e+00, 7.0000e+00,
  5.0000e+01, 1.0000e+03
};

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
