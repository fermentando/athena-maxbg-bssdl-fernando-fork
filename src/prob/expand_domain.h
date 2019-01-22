/*
Defining the scaling exponents for domain rescaling.
*/


// Exponents of expansion. Enter in this order
struct ed_exp_s {
  Real rho;
  Real pressure;
  Real vx, vy, vz;
};

/* ======================================== */
/*         Uniformly expanding box          */
/* ======================================== */
// exp_P = Gamma * ed_exp_rho
#define ED_UNIFORM {-3, -5, -1, -1, -1}

/* ======================================== */
/*          Radially expanding box          */
/* ======================================== */
// exp_P = Gamma * ed_exp_rho
#define ED_RADIAL {-2, -10/3., 0, -1, -1}

/* ======================================== */
/*    Radially expanding with T=const.      */
/* ======================================== */
// P = T * rho
#define ED_RADIAL_ISOTHERMAL {-2, -2, 0, -1, -1}

/* ----------------------------------------*/
/*        Define active modes here         */
/* ----------------------------------------*/
static const struct ed_exp_s ed_exp[2] = {ED_RADIAL,ED_RADIAL_ISOTHERMAL};

#undef ED_RADIAL
#undef ED_RADIAL_ISOTHERMAL
#undef ED_UNIFORM
