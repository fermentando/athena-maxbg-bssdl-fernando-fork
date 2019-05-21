/*
  Powerlaw cooling curve with alpha = 1.0
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
/*
In order to have the same Lambda(Tcl) use

sdL = Lambda_old(Tcl) * (Tcl / T0)^(-1)

where
T0 = sdT[0] from above
Lambda_old(Tcl) = 0.00030361776545242194 [with default values]
Tcl = 0.0006666666666666676 [default values]
 */
static const Real sdL[nfit_cool_d][nfit_cool_T] = {{
    4.554266481786323e-06
}};

// alpha_k
static const Real sdexpt[nfit_cool_d][nfit_cool_T] = {{
    1.0
}};


#elif nfit_cool_d == 2
#error "Invalid `nfit_cool_d`."
#endif


