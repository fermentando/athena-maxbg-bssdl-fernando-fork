/*
  Powerlaw cooling curve with alpha = 1.5
  Temperatures in keV
  Lambda in 1e-23 erg/s
*/

#define nfit_cool_d 1
#define nfit_cool_T 1

// Temperature bins
static Real sdT[nfit_cool_T] = {
  1.0e-5,
};

#if nfit_cool_d == 1

// Density bins
static Real sdd[nfit_cool_d] = {
  1e-24
};

// Lambda_k
static const Real sdL[nfit_cool_d][nfit_cool_T] = {{
    0.003713567542367185
}};

// alpha_k
static const Real sdexpt[nfit_cool_d][nfit_cool_T] = {{
    1.5
}};


#elif nfit_cool_d == 2
#error "Invalid `nfit_cool_d`."
#endif


