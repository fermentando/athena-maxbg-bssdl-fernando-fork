/*
Defining the scaling exponents for domain rescaling.
*/


/* ======================================== */
/*         Uniformly expanding box          */
/* ======================================== */
/*
static const Real ed_exp_rho = -3;
static const Real ed_exp_pressure = -5; // Gamma * ed_exp_rho

static const Real ed_exp_vx = -1;
static const Real ed_exp_vy = -1;
static const Real ed_exp_vz = -1;
*/

/* ======================================== */
/*          Radially expanding box          */
/* ======================================== */
static const Real ed_exp_rho = -2;
static const Real ed_exp_pressure = -10/3.; // Gamma * ed_exp_rho

static const Real ed_exp_vx = 0.0;
static const Real ed_exp_vy = -1;
static const Real ed_exp_vz = -1;

