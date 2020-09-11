#include "copyright.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "prob/math_functions.h"

/* ----------- Define options here --------- */
#define REPORT_NANS          // Verbose
#define ENERGY_COOLING       // Cooling
//#define CUSTOM_BC 2          // If defined use custom boundary conditions.
                             //     1 = shifted periodic bcs
                             //     2 = outflow with forced background
//#define INSTANTCOOL
#define ENERGY_HEATING 0    // 0 = no heating, 1 = heat what cooled, 2 = constant heating
/* ----------------------------------------- */

#ifdef MPI_PARALLEL
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

Real randomreal2(Real min, Real max);
static Real RandomNormal2(Real mu, Real sigma);

/* unit vectors */
Real3Vect get_e1(Real theta, Real phi);
Real3Vect get_e2(Real theta, Real phi);
Real3Vect get_e3(Real theta, Real phi);

void init_hst_enrolls();

/* custom hst quantities */
static Real hst_m13(const GridS *pG, const int i, const int j, const int k);
static Real hst_m110(const GridS *pG, const int i, const int j, const int k);
static Real hst_mT2(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mx13(const GridS *pG, const int i, const int j, const int k);
static Real hst_mTfl2(const GridS *pG, const int i, const int j, const int k);

static Real hst_Erad(const GridS *pG, const int i, const int j, const int k);

static Real hst_V13(const GridS *pG, const int i, const int j, const int k);
static Real hst_VT2(const GridS *pG, const int i, const int j, const int k);
static Real hst_VTfl2(const GridS *pG, const int i, const int j, const int k);

/* dye-weighted hst quantities */
#if (NSCALARS > 0)
static Real hst_c(const GridS *pG, const int i, const int j, const int k);

static Real hst_cE(const GridS *pG, const int i, const int j, const int k);
static Real hst_c_sq(const GridS *pG, const int i, const int j, const int k);

static Real hst_cx1(const GridS *pG, const int i, const int j, const int k);

static Real hst_cvx(const GridS *pG, const int i, const int j, const int k);
static Real hst_cvy(const GridS *pG, const int i, const int j, const int k);
static Real hst_cvz(const GridS *pG, const int i, const int j, const int k);

static Real hst_cvx_sq(const GridS *pG, const int i, const int j, const int k);
static Real hst_cvy_sq(const GridS *pG, const int i, const int j, const int k);
static Real hst_cvz_sq(const GridS *pG, const int i, const int j, const int k);

#ifdef ENERGY_COOLING
static Real hst_cstcool(const GridS *pG, const int i, const int j, const int k);
#endif
static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k);
#endif  /* NSCALARS */

#ifdef CUSTOM_BC
#if CUSTOM_BC == 1 /* shifted periodic */
static int _shift_index(int i, int is, int ie, int ish);
static void bc_shifted_periodic_ix1(GridS *pGrid);
static void bc_shifted_periodic_ox1(GridS *pGrid);
static int ishift, jshift, kshift;
#elif CUSTOM_BC == 2
static void bc_outflowmod_ix1(GridS *pGrid);
static void bc_outflowmod_ox1(GridS *pGrid);
#else
#error "Unknown CUSTOM_BC."
#endif             /* end shifted periodic */
#endif // CUSTOM_BC


#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix);
static OutputS nan_dump;
static int nan_dump_count;
#endif  /* REPORT_NANS */

#ifdef ENERGY_COOLING
/* global definitions for the SD cooling curve using the
   Townsend (2009) exact integration scheme */

#include "prob/cooling_data/SD93_Z1.h"
//#include "prob/cooling_data/denstest.h"
//#include "prob/cooling_data/WSS09_z0_Z1.h"
//#include "prob/cooling_data/powerlaw/alpha1.5.h"

static Real Yk[nfit_cool_d][nfit_cool_T];
/* -- end piecewise power-law fit */

/* must call init_cooling() in both problem() and read_restart() */
static void init_cooling();
static void test_cooling();

static Real sdLambda(const Real d, const Real T);
static Real tcool(const Real d, const Real T);

static Real Y(const Real T, const int id);
static Real Yinv(const Real Y1, const int id);

static Real newtemp_townsend(const Real d, const Real T, const Real dt_hydro);

static void integrate_cooling(GridS *pG);
#if ENERGY_HEATING == 1
static void radiate_energy(MeshS *pM);
#endif

#endif  /* ENERGY_COOLING */

#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3);

static Real maxnu;
#endif  /* VISCOSITY */
#if ENERGY_HEATING == 2 // fixed heating
static Real heating_rate;
#endif /* ENERGY_HEATING == 2 */


static Real drat, dr,dp,tnotcool, r_cloud, acc, pressfac_bkg, pressfac_cl;

static Real tfloor, tceil, rhofloor, betafloor, tfloor_cooling; /* Used in nancheck*/
static Real dens_conv;

static Real dtmin;
static Real pro(Real r, Real rcloud)
{
  return (r/rcloud - log(cosh(r/rcloud))) / log(2);
}


static Real get_pressure(ConsS *u) {
  Real E0 = 0.5 * (SQR(u->M1) + SQR(u->M2) + SQR(u->M3)) / u->d;
  return (u->E - E0) * Gamma_1;
}


/* 
 Helper function to initialize all the hst enrolls.
 Is called for init & restart.
 */
void init_hst_enrolls() {
    dump_history_enroll(hst_m13, "m13");
  //dump_history_enroll(hst_m110, "m110");
  dump_history_enroll(hst_mT2, "mT2");
  dump_history_enroll(hst_mTfl2, "mTfl2");
  dump_history_enroll(hst_Mx13, "Mx13");

  dump_history_enroll(hst_VT2, "VT2");
  dump_history_enroll(hst_V13, "V13");
  dump_history_enroll(hst_VTfl2, "VTfl2");

  dump_history_enroll(hst_Erad, "Erad");

#if (NSCALARS > 0)
  dump_history_enroll(hst_c,    "<c>");
  dump_history_enroll(hst_c_sq, "<c^2>");

  dump_history_enroll(hst_cE,  "<c * E>");
  //dump_history_enroll(hst_cx1, "<c * x1>");

  dump_history_enroll(hst_cvx, "<c * Vx>");
  dump_history_enroll(hst_cvy, "<c * Vy>");
  //dump_history_enroll(hst_cvz, "<c * Vz>");

  dump_history_enroll(hst_cvx_sq, "<(c * Vx)^2>");
  dump_history_enroll(hst_cvy_sq, "<(c * Vy)^2>");
  //dump_history_enroll(hst_cvz_sq, "<(c * Vz)^2>");

  dump_history_enroll(hst_Sdye, "dye entropy");
  /*
#ifdef ENERGY_COOLING
  dump_history_enroll(hst_cstcool, "cs*tcool");
#endif
  */
#endif /* NSCALARS */

}


/*==============================================================================
 * INITIAL CONDITION:
 *
 *----------------------------------------------------------------------------*/
void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0, ii = 0;
  int is,ie,js,je,ks,ke;
  int il,iu,jl,ju,kl,ku;
  Real x,y,z;

  int iseed;
  int nterms,nx1,nx2,nx3,ierr;
  Real scal[5];
  Real my_scal[5];

  Real x1min, x1max, x2min, x2max, x3min, x3max, tmp;

  int cloud_geometry;
#if (NSCALARS > 0)
  Real dye;
#endif

  const Real tconv = 6e4 / ((Gamma_1 - 0.6) / 1e2);

  Real chi   = par_getd("problem", "chi");
  Real X     = par_getd("problem", "X");

  tfloor_cooling = par_getd_def("problem", "tfloor", 2e4) / tconv;
  Real tcold = X * tfloor_cooling;
  tfloor = 0.9 * tfloor_cooling;
  Real thot = chi * tcold;
  tnotcool = 0.6 * thot;

  tceil = par_getd_def("problem", "tceil", 8e10) / tconv;
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-7);
  dtmin = par_getd_def("problem", "dtmin", 1.e-7);

  Real Lcold = -par_getd("domain1", "x1min");

  Real perturb_sigma = par_getd_def("problem", "perturb_sigma", -1);
  Real perturb_max = par_getd_def("problem", "perturb_max", 0.03);

  Real rho_h = 1.;
  Real rho_c = chi * rho_h;

  iseed = -10;
#ifdef MPI_PARALLEL
  iseed -= myID_Comm_world;
#endif
  srand(iseed);


#ifdef ENERGY_COOLING
  init_cooling();
  /* test_cooling(); */
#endif

  // Set up all the hst output
  init_hst_enrolls();

  ath_pout(0, "[init_problem] chi = %g, X = %g, T_hot = %e, t_cool_cold = %e, t_cool_hot = %e, "
           "t_sc_floor = %e\n",
           chi, X, thot, 
           tcool(rho_c, tcold), tcool(rho_h, thot),
           Lcold / sqrt(tcold / X)
           );


#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif


  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  int iprint = 0;

  Real bound_amp = 0.01 * Lcold;
  Real bound_lambda = 0.1 * Lcold;

  /* Begin cell loop */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x,&y,&z);


        Real dr = 0.025 * Lcold;
        //Real dens = (rho_c-rho_h)*exp(alpha_exp*f)*(0.5 - atan(x/dr)/PI) + rho_h;
        Real x_boundary = bound_amp * sin(PI*y / bound_lambda); //-5.0+5.0*sin(0.1*PI*y)
        Real dens = (rho_c-rho_h)* 1.0 * (0.5 - atan((x-x_boundary)/dr)/PI) + rho_h;

        // Constant pressure everywhere
        Real press = rho_h*thot;

        // Cold already at low pressure
        //Real press_h = rho_h*thot;
        //Real press_c = press_h / X;
        //Real press = (press_c-press_h)* 1.0 * (0.5 - atan((x-x_boundary)/dr)/PI) + press_h;

	// Cold gas overpressured (but at floor)
        //Real press_h = rho_h*thot / X;
        //Real press_c = press_h * X;
        //Real press = (press_c-press_h)* 1.0 * (0.5 - atan((x-x_boundary)/dr)/PI) + press_h;


        /* write values to the grid */
        pGrid->U[k][j][i].d = dens;
        pGrid->U[k][j][i].E = press / Gamma_1;
      }
    }
  } /* end grid loops */



  /* seed a perturbation */
  Real fact;
  if(perturb_sigma > 0) {
    ath_pout(0, "[setup] Seeding perturbations with perturb_sigma = %e and sigma_max = %e\n",
	     perturb_sigma, perturb_max);
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          fact = -1.0;
          while (fabs(fact) > perturb_max)
            fact = (RandomNormal(0.0, perturb_sigma));
          if(fabs(fact) < perturb_max)
            pGrid->U[k][j][i].d *= (1.0+fact);
        }
      }
    }
  }

  ath_pout(0, "[setup] done.");

}




/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{

  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;



#ifdef ENERGY_COOLING
  init_cooling();
#endif

#ifdef INSTANTCOOL
  //  CoolingFunc = instant_cool;
#endif

  // Re-enroll hst dumps after restart
  init_hst_enrolls();

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_before_loop(MeshS *pM)
{
  int nl, nd, ntot;

  /* report nans first, so we can fix them before they propagate into
     the following functions. */
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef REPORT_NANS
        ntot = report_nans(pM, &(pM->Domain[nl][nd]),1);
        //if(ntot > 0)
        //report_nans(pM, &(pM->Domain[nl][nd]),1);
#endif
#ifdef INSTANTCOOL
        after_cool(pM, &(pM->Domain[nl][nd]),1);
#endif
      }
    }
  }

  return;
}

void Userwork_in_loop(MeshS *pM)
{
  int nl, nd, ntot;


  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef ENERGY_COOLING
        integrate_cooling(pM->Domain[nl][nd].Grid);
#endif
      }
    }
  }

#if ENERGY_HEATING == 1
  radiate_energy(pM);
#endif


  return;
}

void Userwork_after_loop(MeshS *pM)
{

  return;
}



/*==============================================================================
 * PHYSICS FUNCTIONS:
 * boost_frame()         - boost simulation frame by a velocity increment
 * ed_exp_mode()         - returns current expansion mode
 * expand_domain()       - expands domain by scalefactor
 * report_nans()         - apply a ceiling and floor to the temperature
 * cloud_velocity()      - find the mass-weighted velocity of the cloud.
 * nu_fun()              - kinematic viscosity (i.e., cm^2/s)
 * cooling_func()        - cooling function for the energy equation
 *----------------------------------------------------------------------------*/

#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix)
{
#ifndef ISOTHERMAL
  int i, j, k;
  int is,ie,js,je,ks,ke;
  Real x1, x2, x3;
  int V=1; //verbose off = 0
  int NO = 5;
  Real KE, rho, press, temp;
  int nanpress=0, nanrho=0, nanv=0, nnan;   /* nan count */
  int npress=0,   nrho=0,   nv=0,   nfloor; /* floor count */
  Real scal[8];
  Real beta = 0;
#ifdef MPI_PARALLEL
  Real my_scal[8];
  int ierr;
#endif

  /*Real tfloor    = 1.0e-2 / drat;
  Real tceil     = 100.0;
  Real rhofloor  = 1.0e-2;
  Real betafloor = 3.0e-3;
  */
  GridS *pGrid = pDomain->Grid;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        rho = pGrid->U[k][j][i].d;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        KE = (SQR(pGrid->U[k][j][i].M1) +
              SQR(pGrid->U[k][j][i].M2) +
              SQR(pGrid->U[k][j][i].M3)) /
          (2.0 * rho);

        press = pGrid->U[k][j][i].E - KE;
        press *= Gamma_1;
        temp = press / rho;

        if (press != press) {
          nanpress++;
          if(V && nanpress < NO) printf("bad press %e R %e  %e %e %e  %d %d %d %e %e %e %e\n",press, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta);
          if(fix)
            temp = tfloor;
        } else if (temp < tfloor) {
          npress++;
          if(V && npress < NO) printf("bad tempF %e R %e  %e %e %e  %d %d %d %e %e %e %e\n",temp, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta);
          if(fix)
            temp = tfloor;
        } else if (temp > tceil && beta > 10.0 * betafloor) {
          npress++;
          if(V && npress < NO) printf("bad tempC %e R %e  %e %e %e  %d %d %d %e %e %e %e\n",temp, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta);
          if(fix)
            temp = tceil;
        }

        if (rho != rho) {
          nanrho++;
          if(V&&nanrho < NO) printf("bad rho %e R %e  %e %e %e  %d %d %d\n",rho, 1.0/pGrid->dx1, x1, x2, x3, i, j, k);
          if(fix)
            rho = rhofloor;
        } else if (rho < rhofloor) {
          nrho++;
          if(V&& nrho < NO) printf("bad rho %e R %e  %e %e %e  %d %d %d\n",rho, 1.0/pGrid->dx1, x1, x2, x3, i, j, k);
          if(fix)
            rho = rhofloor;
        }

        if (pGrid->U[k][j][i].M1 != pGrid->U[k][j][i].M1) {
          nanv++;
          if(fix)
            pGrid->U[k][j][i].M1 = 0.0;
        }
        if (pGrid->U[k][j][i].M2 != pGrid->U[k][j][i].M2) {
          nanv++;
          if(fix)
            pGrid->U[k][j][i].M2 = 0.0;
        }
        if (pGrid->U[k][j][i].M3 != pGrid->U[k][j][i].M3) {
          nanv++;
          if(fix)
            pGrid->U[k][j][i].M3 = 0.0;
        }

        /* write values back to the grid */
        if(fix) {
          pGrid->U[k][j][i].d  = rho;
          KE = (SQR(pGrid->U[k][j][i].M1) +
                SQR(pGrid->U[k][j][i].M2) +
                SQR(pGrid->U[k][j][i].M3)) /
            (2.0 * rho);

          pGrid->U[k][j][i].E = temp*rho/Gamma_1 + KE;
        }
      }
    }
  }

  /* synchronize over grids */
#ifdef MPI_PARALLEL
  my_scal[0] = nanpress;
  my_scal[1] = nanrho;
  my_scal[2] = nanv;
  my_scal[4] = npress;
  my_scal[5] = nrho;
  my_scal[6] = nv;

  ierr = MPI_Allreduce(&my_scal, &scal, 8, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);

  nanpress = scal[0];
  nanrho   = scal[1];
  nanv     = scal[2];
  npress   = scal[4];
  nrho     = scal[5];
  nv       = scal[6];
#endif  /* MPI_PARALLEL */


  /* sum up the # of bad cells and report */
  nnan = nanpress+nanrho+nanv;

  /* sum up the # of floored cells and report */
  nfloor = npress+nrho+nv;


  if(nfloor > 0) {
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v.\n",
             nfloor, npress, nrho, nv);
  }


  //  if ((nnan > 0 || nfloor -nmag > 30) && fix == 0) {
    if (nnan > 0 ){// && fix == 0) {
    ath_pout(0, "[report_nans]: found %d nan cells: %d P, %d d, %d v.\n",
             nnan, nanpress, nanrho, nanv);

    nan_dump.n      = 100;
    nan_dump.dt     = HUGE_NUMBER;
    nan_dump.t      = pM->time;
    nan_dump.num    = 1000 + nan_dump_count;
    nan_dump.out    = "prim";
    nan_dump.nlevel = -1;       /* dump all levels */


    dump_vtk(pM, &nan_dump);
    if(nnan) nan_dump_count++;
    if (nan_dump_count > 10)
      ath_error("[report_nans]: too many nan'd timesteps.\n");

    if (nfloor > 1000)
      ath_error("[report_nans]: Too many floored cells.\n");
  }


#endif  /* ISOTHERMAL */

  return nfloor+nnan;
}
#endif  /* REPORT_NANS */


#ifdef ENERGY_COOLING
/* ================================================================ */
/* cooling routines */

static void init_cooling()
{
  int i, k, n=nfit_cool_T-1;
  Real term;
  const Real mu = 0.62, mu_e = 1.17;

  /* convert T in the cooling function from keV to code units */
  for (k=0; k<=n; k++) {
    sdT[k] /= (8.197 * mu);
    if(sdT[k] <= 0) ath_error("sdT[%d]=%e. Has to be > 0.", k, sdT[k]);
  }

  if(tfloor_cooling < sdT[0])
    ath_error("Cooling floor is smaller than first entry of cooling function (%e vs %e).",
              tfloor_cooling, sdT[0]);

  /* populate Yk following equation A6 in Townsend (2009) */
  for(i = 0; i < nfit_cool_d; i++) {
    Yk[i][n] = 0.0;
    for (k=n-1; k>=0; k--){
      if(sdL[i][k] <= 0) ath_error("sdL[%d][%d]=%e. Has to be > 0.", i, k, sdL[i][k]);
      term = (sdL[i][n]/sdL[i][k]) * (sdT[k]/sdT[n]);

      if (sdexpt[i][k] == 1.0)
        term *= log(sdT[k]/sdT[k+1]);
      else
        term *= ((1.0 - pow(sdT[k]/sdT[k+1], sdexpt[i][k]-1.0)) / (1.0-sdexpt[i][k]));

      Yk[i][k] = Yk[i][k+1] - term;

      if(isnan(Yk[i][k]))
        ath_error("Error initializing cooling. nan in Yk[%d][%d]", i, k);
    }
  }
  return;
}

/* piecewise power-law fit to the cooling curve with temperature in
   keV and L in 1e-23 erg cm^3 / s */
static Real sdLambda(const Real d0, const Real T)
{
  int iT, id; // bin indices for T,d
  Real L1, L2;
  const Real conv_fac = 1.311e-5; // from units of 1e-23 erg cm^3 /s to code units.
  const Real d = d0 * dens_conv;
  int interpolate = (nfit_cool_d > 1);

  /* first find the temperature bin */
  for(iT=nfit_cool_T-1; iT>=0; iT--){
    if (T >= sdT[iT])
      break;
  }
  if(iT < 0)
    ath_error("[sdLambda] T %e d %e %e %d\n", T, d, sdT[0], iT);

  /* Find the density bin */
  if(d <= sdd[0]) {
    interpolate = 0;
    id = 0;
  } else if(d >= sdd[nfit_cool_d - 1]) {
    interpolate = 0;
    id = nfit_cool_d - 1;
  } else {
    for(id=nfit_cool_d - 1; id>=0; id--) {
      if(d >= sdd[id])
        break;
    }
  }

  /* piecewise power-law; see equation A4 of Townsend (2009) */
  L1 = conv_fac * sdL[id][iT] * pow(T/sdT[iT], sdexpt[id][iT]);

  if(!interpolate) 
    return L1;

  L2 = conv_fac * sdL[id+1][iT] * pow(T/sdT[iT], sdexpt[id+1][iT]);

  // Linear interpolation
  return L1 + (L2 - L1) / (sdd[id+1] - sdd[id]) * d;

}

static Real tcool(const Real d, const Real T)
{
  const Real mu = 0.62, mu_e = 1.17;

  /* equation 13 of Townsend (2009) */
  return (SQR(mu_e) * T) / (Gamma_1 * d * sdLambda(d,T));
}


/* see sdLambda() or equation A1 of Townsend (2009) for the
   definition */
static Real Y(const Real T, const int id)
{
  int iT;
  int nT = nfit_cool_T - 1;
  int nd = nfit_cool_d - 1;
  Real term;

  /* first find the temperature bin */
  for(iT=nT; iT>=0; iT--){
    if (T >= sdT[iT])
      break;
  }

  /* calculate Y using equation A5 in Townsend (2009) */
  term = (sdL[id][nT]/sdL[id][iT]) * (sdT[iT]/sdT[nT]);

  if (sdexpt[id][iT] == 1.0)
    term *= log(sdT[iT]/T);
  else
    term *= ((1.0 - pow(sdT[iT]/T, sdexpt[id][iT]-1.0)) / (1.0-sdexpt[id][iT]));

  return (Yk[id][iT] + term);
}


static Real Yinv(const Real Y1, const int id) {
  //int iT,id;
  int nT = nfit_cool_T - 1;
  int nd = nfit_cool_d - 1;
  int iT;
  Real term;

  /* find the bin i in which the final temperature will be */
  for(iT=nT; iT>0; iT--){       // use iT>0 instead of iT>=0 to force min(iT)=0
    if (Y(sdT[iT], id) >= Y1)   // this means that newT<sdT[0] are wrong
      break;                    // but since we have Tfloor>sdT[0] we're good.
  }

  /* calculate Yinv using equation A7 in Townsend (2009) */
  term = (sdL[id][iT]/sdL[id][nT]) * (sdT[nT]/sdT[iT]);
  term *= (Y1 - Yk[id][iT]);

  if (sdexpt[id][iT] == 1.0)
    term = exp(-1.0*term);
  else{
    term = pow(1.0 - (1.0-sdexpt[id][iT])*term,
               1.0/(1.0-sdexpt[id][iT]));
  }


  if(!isfinite(term)) {
    ath_error("nan detected in Yinv (term=%e, Y1=%e).\n", term, Y1);
  }

  return (sdT[iT] * term);
}


static Real newtemp_townsend(const Real d0, const Real T, const Real dt_hydro)
{
  int id;
  Real term1, Tref, dref;
  Real T1, T2;
  const Real d = d0 * dens_conv;
  int interpolate = nfit_cool_d > 1;

  if(T <= tfloor_cooling)
    return tfloor_cooling;

  Tref = sdT[nfit_cool_T-1];
  dref = sdd[nfit_cool_d-1] / dens_conv;

  /* Find the density bin */
  if(d <= sdd[0]) {
    interpolate = 0;
    id = 0;
  } else if(d >= sdd[nfit_cool_d - 1]) {
    interpolate = 0;
    id = nfit_cool_d - 1;
  } else {
    for(id=nfit_cool_d - 1; id>=0; id--) {
      if(d >= sdd[id])
        break;
    }
  }
  assert(id>=0);
  assert(id<nfit_cool_d);

  /* Calculate new T */
  term1 = (T/Tref) * (sdLambda(d0, Tref)/sdLambda(d0, T)) * (dt_hydro/tcool(d0, T));
  T1 = Yinv(Y(T,id) + term1, id);
  if(isnan(T1)) {
    printf("[newtemp] d=%.3e, id=%d, T=%.3e --> %.3e, dt_hydro=%e, tcool=%e, Tref=%e, sdLambda=%e, term1=%e, Y(T,id)=%e\n", d, id, T, T1,dt_hydro, tcool(d0,T),Tref, sdLambda(d0, T), term1, Y(T,id));
    assert(!isnan(T1));
  }
  if(!interpolate)
    return T1;

  T2 = Yinv(Y(T,id+1) + term1, id+1);
  assert(!isnan(T2));

  /* Linear interpolation along density */
  return T1 + (T2 - T1) / (sdd[id+1] - sdd[id]) * d;
}


static void integrate_cooling(GridS *pG)
{
  int i, j, k, iprint = 0;
  int is, ie, js, je, ks, ke;

  PrimS W;
  ConsS U;
  // Changed this for tfloor!!
  Real temp, tempold, heat;

  /* ath_pout(0, "integrating cooling using Townsend (2009) algorithm.\n"); */

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        W = Cons_to_Prim(&(pG->U[k][j][i]));
        pG->U[k][j][i].Erad = 0;

        /* find temp in keV */
        temp = W.P/W.d;
        tempold = temp;

        /* do not cool above a certain threshold */
        if( (tnotcool > 0) && (temp > tnotcool) )
          continue;

        temp = newtemp_townsend(W.d, temp, pG->dt);

        /* apply a temperature floor (nans tolerated) */
        if (isnan(temp) || temp < tfloor_cooling)
          temp = tfloor_cooling;

        W.P = W.d * temp;
        U = Prim_to_Cons(&W);

        /* record cooled energy */
	// pG->U[k][j][i].Erad += (pG->U[k][j][i].E - U.E);
        pG->U[k][j][i].Erad = (pG->U[k][j][i].E - U.E) / pG->dt;

        pG->U[k][j][i].E = U.E;

      }
    }
  }

  return;

}


#if ENERGY_HEATING == 1
/*
  Radiate cooled energy over whole domain
 */
static void radiate_energy(MeshS *pM) {
  //  ath_error("Does not work right now.");
  Real Erad_total = 0;
  GridS *pG;
  int i, j, k, is, ie, js, je, ks, ke;
  int nl, nd;
  Real dV;
  Real V = (pM->RootMaxX[2] - pM->RootMinX[2]) * \
    (pM->RootMaxX[1] - pM->RootMinX[1]) * \
    (pM->RootMaxX[0] - pM->RootMinX[0]);

#ifdef MPI_PARALLEL
  int ierr;
  Real my_Etot;
#endif

  // Calculate total energy
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG = pM->Domain[nl][nd].Grid;

        is = pG->is;  ie = pG->ie;
        js = pG->js;  je = pG->je;
        ks = pG->ks;  ke = pG->ke;

        for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
              Erad_total += pG->U[k][j][i].Erad;
              /*
                Let's not set this to zero here but instead after the output so we
                can visualize Erad.
               */
              // pG->U[k][j][i].Erad = 0.0;
            }
          }
        }
      }
    }
  }
  //  ath_pout(0,"Local E %e\n", Erad_total);

#ifdef MPI_PARALLEL
  my_Etot = Erad_total;

  ierr = MPI_Allreduce(&my_Etot, &Erad_total, 1, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[radiate_energy]: MPI_Allreduce returned error %d\n", ierr);
#endif

  // Distribute energy over grid
  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG = pM->Domain[nl][nd].Grid;

        is = pG->is;  ie = pG->ie;
        js = pG->js;  je = pG->je;
        ks = pG->ks;  ke = pG->ke;

        dV = pG->dx1 * pG->dx2 * pG->dx3;
        for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
              pG->U[k][j][i].E += Erad_total * dV / V;
            }
          }
        }
      }
    }
  }
  ath_pout(0, "[radiate_energy] Distributed a total energy of %e in fractional "
           "volumes of %e.\n", Erad_total, dV / V);

}
#endif /* ENERGY_HEATING */


static void test_cooling()
{
  int i, j, npts=100;
  Real logt, temp, tc, logdt, dt;
  Real err;
  Real dens = 1.0;

  FILE *outfile;

  /* this sometimes crashes, so let's try it */
  newtemp_townsend(sdd[nfit_cool_d - 1] - 1e-4, 7.173e-04, 8.050312e-03);

  /* Now outputting some information from Townsend (2009) */
  outfile = fopen("lambda.dat", "w");
  for(i=0; i<npts; i++){
    logt = log(1.0e-4) + (log(5.0)-log(1.0e-4))*((double) i/(npts-1));
    temp = exp(logt);

    fprintf(outfile, "%e\t%e\t%e\n", temp, sdLambda(dens,temp), sdLambda(dens * drat,temp));
  }
  fclose(outfile);


  outfile = fopen("Yk.dat", "w");
  for(i=0; i<nfit_cool_T; i++) {
    fprintf(outfile, "%d\t", i);
    for(j = 0; j < nfit_cool_d; j++)
      fprintf(outfile, "%e\t", Yk[j][i]);
    fprintf(outfile, "\n");
  }
  fclose(outfile);

  temp = 10.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-10kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);
    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 3.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-3kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 1.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-1kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 0.3;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.3kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 0.1;
  tc = tcool(1.0, temp);


  outfile = fopen("townsend-fig1-0.1kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 0.1;
  tc = tcool(100.0, temp);

  outfile = fopen("townsend-fig1-0.1kev-100.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(100.0, temp, dt));
  }
  fclose(outfile);


  ath_error("check cooling stuff done.\n");

  return;
}
/* end cooling routines */
/* ================================================================ */
#endif  /* ENERGY_COOLING */

#ifdef INSTANTCOOL
static Real instant_cool(const Real rho, const Real P, const Real dt)
{
  /*returns cooling to decrease temperature to initial cloud temperature */
  Real temp = P/rho;
  Real tcloud = Gamma_1 / drat;
  Real Edot = 0.0;
  /* cool to tcloud, but no further */
  if((temp > tcloud) && (temp < 10.*tcloud)){
    Edot = (temp - tcloud)*rho/dt/Gamma_1;
  }
  //Edot = (temp < tcloud) ? 0.0 : (temp-tcloud)/dt;


  return Edot;
}



static int after_cool(MeshS *pM, DomainS *pDomain, int fix)
{
  int i, j, k;
  int is,ie,js,je,ks,ke;
  Real x1, x2, x3;
  int V=0; //verbose off
  int NO = 2;
  Real KE, rho, press, temp;

#ifdef MHD
  Real ME;
  int nanmag=0;
  int nmag=0;
#endif  /* MHD */

  GridS *pGrid = pDomain->Grid;
  Real tcloud = Gamma_1 / drat;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        rho = pGrid->U[k][j][i].d;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        KE = (SQR(pGrid->U[k][j][i].M1) +
              SQR(pGrid->U[k][j][i].M2) +
              SQR(pGrid->U[k][j][i].M3)) /
          (2.0 * rho);

        press = pGrid->U[k][j][i].E - KE;
#ifdef MHD
        ME = (SQR(pGrid->U[k][j][i].B1c) +
              SQR(pGrid->U[k][j][i].B2c) +
              SQR(pGrid->U[k][j][i].B3c)) * 0.5;
        press -= ME;
#endif  /* MHD */

        press *= Gamma_1;
        temp = press / rho;


        if((temp  > 1.1*tcloud) && (temp < 10. * tcloud) && (pGrid->U[k][j][i].s[0] >= 0.1)){
          temp = 1.1*tcloud;

          pGrid->U[k][j][i].E = temp*rho/Gamma_1 + KE;
#ifdef MHD
          pGrid->U[k][j][i].E += ME;
#endif  /* MHD */
        }
      }
    }
  }
  return 0;
}

#endif /* INSTANTCOOL */


/*==============================================================================
 * HISTORY OUTPUTS:
 *

 *----------------------------------------------------------------------------*/

static Real _hst_mcut(const GridS *pG, const int i, const int j, const int k, const Real frac)
{
#ifdef EXPAND_DOMAIN
  Real s = pow(scalefac, ed_exp[ed_exp_mode()].rho);
#else
  Real s = 1.0;
#endif
  if(pG->U[k][j][i].d < frac * drat * s * pressfac_cl)
    return 0;
  return pG->U[k][j][i].d;
}


static Real hst_m13(const GridS *pG, const int i, const int j, const int k)
{
  return _hst_mcut(pG, i, j, k, 1/3.);
}

static Real hst_mT2(const GridS *pG, const int i, const int j, const int k)
{
  Real temp = get_pressure(&(pG->U[k][j][i])) / pG->U[k][j][i].d;
  const Real Tcl = (Gamma_1 + dp) / drat;
  if(temp > 2 * Tcl)
    return 0;
  return pG->U[k][j][i].d;
}


static Real hst_mTfl2(const GridS *pG, const int i, const int j, const int k)
{
  Real temp = get_pressure(&(pG->U[k][j][i])) / pG->U[k][j][i].d;
  const Real Tfloor = MAX(tfloor, tfloor_cooling);
  if(temp > 2 * Tfloor)
    return 0;
  return pG->U[k][j][i].d;
}


static Real hst_m110(const GridS *pG, const int i, const int j, const int k)
{
  return _hst_mcut(pG, i, j, k, 0.1);
}


static Real hst_Mx13(const GridS *pG, const int i, const int j, const int k)
{
  Real s = 1.0;
  if(pG->U[k][j][i].d < pressfac_cl * drat / 3. * s)
    return 0;
  return pG->U[k][j][i].M1;
}

static Real hst_Erad(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].Erad;
}

static Real hst_VTfl2(const GridS *pG, const int i, const int j, const int k)
{
  Real temp = get_pressure(&(pG->U[k][j][i])) / pG->U[k][j][i].d;
  const Real Tfloor = MAX(tfloor, tfloor_cooling);
  if(temp > 2 * Tfloor)
    return 0;
  return 1.0;
}


static Real hst_VT2(const GridS *pG, const int i, const int j, const int k)
{
  Real temp = get_pressure(&(pG->U[k][j][i])) / pG->U[k][j][i].d;
  const Real Tcl = (Gamma_1 + dp) / drat;
  if(temp > 2 * Tcl)
    return 0;
  return 1.0;
}

static Real hst_V13(const GridS *pG, const int i, const int j, const int k)
{
#ifdef EXPAND_DOMAIN
  Real s = pow(scalefac, ed_exp[ed_exp_mode()].rho);
#else
  Real s = 1.0;
#endif
  if(pG->U[k][j][i].d < drat * s * pressfac_cl /3.)
    return 0;
  return 1.0;
}



/* dye-weighted hst quantities */
#if (NSCALARS > 0)
static Real hst_c(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0]);
}

static Real hst_c_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0]));
}

static Real hst_cE(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].E);
}

static Real hst_cx1(const GridS *pG, const int i, const int j, const int k)
{
  Real x1, x2, x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return (pG->U[k][j][i].s[0] * x1);
}

static Real hst_cvx(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M1  /  pG->U[k][j][i].d);
}

static Real hst_cvy(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M2  /  pG->U[k][j][i].d);
}

static Real hst_cvz(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M3 /  pG->U[k][j][i].d);
}


static Real hst_cvx_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0] * pG->U[k][j][i].M1  /  pG->U[k][j][i].d));
}

static Real hst_cvy_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0] * pG->U[k][j][i].M2  /  pG->U[k][j][i].d));
}

static Real hst_cvz_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0] * pG->U[k][j][i].M3 /  pG->U[k][j][i].d));
}

static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k)
{
  Real dye = pG->U[k][j][i].s[0] / pG->U[k][j][i].d;
  Real rho = pG->U[k][j][i].d;

  if (dye < TINY_NUMBER)
    return 0.0;

  return (-1.0 * rho * dye * log(dye));
}

#ifdef ENERGY_COOLING
static Real hst_cstcool(const GridS *pG, const int i, const int j, const int k)
{
  PrimS W = Cons_to_Prim(&(pG->U[k][j][i]));
  Real temp = W.P / W.d;
  Real cs = sqrt(Gamma * temp);
  return cs * tcool(W.d, temp);
  //return 1.;
}
#endif


#endif  /* NSCALARS */



Real randomreal2(Real min, Real max)
{
  Real eta = ((Real)rand()/(Real)RAND_MAX);
  return min + eta * (max-min);
}

static Real RandomNormal2(Real mu, Real sigma)
/* Implements the box-muller routine.  Gives a mean of mu, and a
   standard deviation sigma.  */
{
  Real x1, x2, w, y1;
  static Real y2;
  static int use_last = 0;

  if (use_last){ /* use value from previous call */
    y1 = y2;
    use_last = 0;
  }
  else {
    do {
      x1 = randomreal2(-1.0, 1.0);
      x2 = randomreal2(-1.0, 1.0);
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }

  return (mu + y1 * sigma);
}
