#include "copyright.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#include "prob/math_functions.h"

#ifdef MPI_PARALLEL
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

#ifdef PARTICLES
#include "particles/particle.h"
#endif

#define FOLLOW_CLOUD
#define REPORT_NANS
#define ENERGY_COOLING
//#define FLOW_PROFILE   // Uncomment this line for changing v(r), rho(r),...
//#define ENERGY_HEATING // Uncomment this line for heating
/* #define INSTANTCOOL */
static void bc_ix1(GridS *pGrid);
static void bc_ox1(GridS *pGrid);
//static void bc_ix2(GridS *pGrid);
//static void bc_ox2(GridS *pGrid);
//static void bc_ix3(GridS *pGrid);
//static void bc_ox3(GridS *pGrid);

static void check_div_b(GridS *pGrid);


/* add a force-free term to B */
void add_term(Real3Vect ***A, GridS *pG,
              Real theta, Real phi,
              Real alpha, Real beta, Real amp);

Real randomreal2(Real min, Real max);
static Real RandomNormal2(Real mu, Real sigma);

/* unit vectors */
Real3Vect get_e1(Real theta, Real phi);
Real3Vect get_e2(Real theta, Real phi);
Real3Vect get_e3(Real theta, Real phi);

/* custom hst quantities */
static Real hst_m13(GridS *pG, int i, int j, int k);
static Real hst_m110(GridS *pG, int i, int j, int k);
static Real hst_Mx13(GridS *pG, int i, int j, int k);

static Real hst_Erad(GridS *pG, int i, int j, int k);


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

#ifdef MHD
static Real hst_cBx(const GridS *pG, const int i, const int j, const int k);
static Real hst_cBy(const GridS *pG, const int i, const int j, const int k);
static Real hst_cBz(const GridS *pG, const int i, const int j, const int k);
#endif  /* MHD */

static Real hst_Sdye(const GridS *pG, const int i, const int j, const int k);
#endif  /* NSCALARS */

#ifdef PARTICLES
static int parnumproc;

void init_particles(DomainS *pDomain);
void add_new_particles(MeshS *pM);
#endif


#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix);
static OutputS nan_dump;
static int nan_dump_count;
#endif  /* REPORT_NANS */

#ifdef FLOW_PROFILE
#include "prob/flow_profile/cc85_fit.h"
#endif /* FLOW_PROFILE */

#ifdef ENERGY_COOLING
/* global definitions for the SD cooling curve using the
   Townsend (2009) exact integration scheme */

#include "prob/cooling_data/SD93_Z1.h"
//#include "prob/cooling_data/SD93_Z1_S10.h"
//#include "prob/cooling_data/WSS09_n1_Z1.h"
//#include "prob/cooling_data/WSS09_CIE_Z1.h"

static Real Yk[nfit_cool];
/* -- end piecewise power-law fit */


/* must call init_cooling() in both problem() and read_restart() */
static void init_cooling();
static void test_cooling();

static Real sdLambda(const Real T);
static Real tcool(const Real d, const Real T);

static Real Y(const Real T);
static Real Yinv(const Real Y1);

static Real newtemp_townsend(const Real d, const Real T, const Real dt_hydro);

static void integrate_cooling(GridS *pG);
#ifdef ENERGY_HEATING
static void radiate_energy(MeshS *pM);
#endif
#endif  /* ENERGY_COOLING */

#ifdef INSTANTCOOL
static Real instant_cool(const Real rho, const Real P, const Real dt);
static int after_cool(MeshS *pM, DomainS *pDomain, int fix);
#endif  /* INSTANTCOOL */

#ifdef FOLLOW_CLOUD
static Real cloud_mass_weighted_velocity(MeshS *pM);
static void boost_frame(DomainS *pDomain, Real dv);
static Real x_shift;

static Real hst_xshift(GridS *pG, int i, int j, int k);
static Real hst_vflow(GridS *pG, int i, int j, int k);
#endif

static Real drat, vflow, vflow0, betain, betaout_y, betaout_z,dr,dp,tnotcool, r_cloud;

static Real tfloor, tceil, rhofloor, betafloor, tfloor_cooling; /* Used in nancheck*/

#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3);

static Real nu;
#endif  /* VISCOSITY */

static Real dtmin;
static Real pro(Real r, Real rcloud)
{
  return (r/rcloud - log(cosh(r/rcloud))) / log(2);
}



/*==============================================================================
 * INITIAL CONDITION:
 *
 *----------------------------------------------------------------------------*/
void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int il,iu,jl,ju,kl,ku;
  Real x1,x2,x3, r;
  Real rho, vx, vy, vz, v_turb, v_rot;
  Real fact, l_cloud;

  int iseed;
  Real3Vect ***A;
  int nterms,nx1,nx2,nx3,ierr;
  Real theta, phi, alpha, beta, amp;
  Real scal[5];
  Real my_scal[5];

  Real x1min, x1max, x2min, x2max, x3min, x3max, tmp;

  Real jmag, bmag, JdB, JcBforce, norm, Brms, Bs, Bmax, Bmax_cloud, ncells;
  Real Bin, Bout_z, Bout_y;
  Real bscale;
  int tangled;
#if (NSCALARS > 0)
  Real dye;
#endif

  drat   = par_getd("problem", "drat");
#ifdef FLOW_PROFILE
  vflow = 0.0;
#else
  vflow  = par_getd("problem", "vflow"); // TODO: change here for FLOW_PROFILE
#endif
  vflow0 = vflow;
#ifdef FOLLOW_CLOUD
  x_shift = 0.0;
#endif

#ifdef MHD
  tangled = par_getd_def("problem","tangled",1);
  betain = par_getd("problem", "betain");
  betaout_y = par_getd_def("problem", "betaout_y", 1e20);
  betaout_z = par_getd("problem", "betaout_z");
  betafloor = par_getd_def("problem", "betafloor", 3.e-3);
#endif

  v_turb = par_getd_def("problem", "v_turb", 0.0);
  v_rot  = par_getd_def("problem", "v_rot",  0.0);

  r_cloud = par_getd_def("problem", "r_cloud", 0.25);
  l_cloud = par_getd_def("problem", "l_cloud", 0.0);
  dr = par_getd_def("problem", "dr", 0.0);

  /*Real tfloor    = 1.0e-2 / drat;
  Real tceil     = 100.0;
  Real rhofloor  = 1.0e-2;
  Real betafloor = 3.0e-3;
  */
  // note that there's another tfloor in the cooling function 
  tfloor = par_getd_def("problem", "tfloor", 1.e-2/drat);
  tceil = par_getd_def("problem", "tceil", 100.);
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-2);
  d_MIN = par_getd_def("problem", "d_MIN", 1.e-4);
  dtmin = par_getd_def("problem", "dtmin", 1.e-7);

  dp = par_getd_def("problem", "dp", 0.0);
  tnotcool = par_getd_def("problem", "tnotcool", -1.0);
  tfloor_cooling = par_getd_def("problem", "tfloor_cooling", (Gamma_1 + dp) / drat);

#ifdef VISCOSITY
  nu   = par_getd("problem","nu");
  NuFun_i = NULL;
  NuFun_a = nu_fun;
#endif

#ifdef THERMAL_CONDUCTION
  kappa_iso = par_getd_def("problem","kappa_iso",0.0);
  kappa_aniso = par_getd_def("problem","kappa_aniso",0.0);
#endif

  iseed = -10;
#ifdef MPI_PARALLEL
  iseed -= myID_Comm_world;
#endif
  srand(iseed);

#ifdef FOLLOW_CLOUD
#if (NSCALARS == 0)
  ath_error("[problem]: requires NSCALARS > 0.\n");
#endif
#endif

#ifdef FLOW_PROFILE
  flow_profile_init();
#endif

#ifdef ENERGY_COOLING
  init_cooling();
  /* test_cooling(); */
#endif

#ifdef INSTANTCOOL
  //  CoolingFunc = instant_cool;
#endif

  dump_history_enroll(hst_m13, "m13");
  dump_history_enroll(hst_m110, "m110");
  dump_history_enroll(hst_Mx13, "Mx13");

  dump_history_enroll(hst_Erad, "Erad");

#ifdef FOLLOW_CLOUD
  dump_history_enroll_alt(hst_xshift, "x_shift");
  dump_history_enroll_alt(hst_vflow,  "v_flow");
#endif

#if (NSCALARS > 0)
  dump_history_enroll(hst_c,    "<c>");
  dump_history_enroll(hst_c_sq, "<c^2>");

  dump_history_enroll(hst_cE,  "<c * E>");
  dump_history_enroll(hst_cx1, "<c * x1>");

  dump_history_enroll(hst_cvx, "<c * Vx>");
  dump_history_enroll(hst_cvy, "<c * Vy>");
  dump_history_enroll(hst_cvz, "<c * Vz>");

  dump_history_enroll(hst_cvx_sq, "<(c * Vx)^2>");
  dump_history_enroll(hst_cvy_sq, "<(c * Vy)^2>");
  dump_history_enroll(hst_cvz_sq, "<(c * Vz)^2>");
#ifdef MHD
  dump_history_enroll(hst_cBx, "<c * Bx>");
  dump_history_enroll(hst_cBy, "<c * By>");
  dump_history_enroll(hst_cBz, "<c * Bz>");
#endif /* MHD */
  dump_history_enroll(hst_Sdye, "dye entropy");
#ifdef ENERGY_COOLING
  dump_history_enroll(hst_cstcool, "cs*tcool");
#endif
#endif /* NSCALARS */

#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif

  /* Initialize grid loading */
  int Nx, Ny, Nz, ndim_file;
  int ii, jj, kk;
  Nx = pDomain->Nx[0]; Ny = pDomain->Nx[1]; Nz = pDomain->Nx[2];
  char* fn_rho = par_gets_def("problem", "load_grid_rho", "");
  char* fn_vx;
  Real *cloud_dat_rho, *cloud_dat_vx;
  FILE* fp_rho = NULL; FILE* fp_vx = NULL;
  if((char)fn_rho[0]) {
#ifdef FLOW_PROFILE
    ath_error("FLOW_PROFILE and grid loading not implemented yet.");
#endif
    fn_vx = par_gets("problem", "load_grid_vx");
    ndim_file = par_geti("problem", "load_grid_ndim");
    ath_pout(0, "Loading grid `%s` / `%s` with %d cells (grid cells: [%d,%d,%d])\n",
             fn_rho, fn_vx, ndim_file, Nx, Ny, Nz);
    if((ndim_file > Nx) || (ndim_file > Ny) || (ndim_file > Nz))
      ath_error("[read_grid]: Grid in file larger than ATHENA grid (%d versus %d, %d, %d).\n",
                ndim_file, Nx, Ny, Nz);
    fp_rho = fopen(fn_rho, "r");
    if(fp_rho == NULL) ath_error("[read_grid] Problem loading `%s`.", fn_rho);
    fp_vx = fopen(fn_vx, "r");
    if(fp_vx == NULL) ath_error("[read_grid] Problem loading `%s`.", fn_vx);
    cloud_dat_rho = (Real*)malloc(sizeof(Real) * ndim_file * ndim_file * ndim_file);
    if(cloud_dat_rho == NULL) ath_error("[read_grid] Error allocating memory for rho.");
    fread(cloud_dat_rho, sizeof(Real), ndim_file * ndim_file * ndim_file, fp_rho);
    fread(&rho, sizeof(Real), 1, fp_rho);
    if(!feof(fp_rho)) ath_error("[read_grid] Not eof for rho.");
    cloud_dat_vx = (Real*)malloc(sizeof(Real) * ndim_file * ndim_file * ndim_file);
    if(cloud_dat_vx == NULL) ath_error("[read_grid] Error allocating memory for vx.");
    fread(cloud_dat_vx, sizeof(Real), ndim_file * ndim_file * ndim_file, fp_vx);
    fread(&vx, sizeof(Real), 1, fp_vx);
    if(!feof(fp_vx)) ath_error("[read_grid] Not eof for vx.");
  }   /* end grid loading */


  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  int iprint = 0;
  vx = vy = vz = 0.0;

  /* Begin cell loop */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        //r = sqrt(x1*x1+x2*x2+x3*x3);
        r = sqrt(pow(MIN(MAX(x1 - l_cloud, 0), x1), 2) + x2*x2 + x3*x3);

#ifdef FLOW_PROFILE
        rho = flow_profile_density(x1, x2, x3);
        vx  = flow_profile_velocity_x(x1, x2, x3);
        vy  = flow_profile_velocity_y(x1, x2, x3);
        vz  = flow_profile_velocity_z(x1, x2, x3);


        /* if((j == js) && (k== ks)) */
        /*   printf("%.2f --> (%e, %e)\n", x1, rho, vx); */

#else // Static inflow
        rho  = 1.0;
        // Fix pressure so that temperature is normed to what it was with drat=1e3
        // dp = (drat / 1000. - 1.) * Gamma_1; // --> does not keep t_cool constant!
        vx   = vflow;
#endif

#if (NSCALARS > 0)
        dye = 0.0;
#endif

        /* Initialize cloud */
        if(fp_rho != NULL) { // read grid from file
          ii = (int)(x1 / pGrid->dx1) + ndim_file / 2; //i - ((Nx - ndim_file) / 2);
          jj = (int)(x2 / pGrid->dx2) + ndim_file / 2; //j - ((Ny - ndim_file) / 2);
          kk = (int)(x3 / pGrid->dx3) + ndim_file / 2; //k - ((Nz - ndim_file) / 2);
          if((ii > 0) && (ii < ndim_file) &&
             (jj > 0) && (jj < ndim_file) &&
             (kk > 0) && (kk < ndim_file)) {
            rho = cloud_dat_rho[kk * ndim_file * ndim_file + jj * ndim_file + ii];
            vx  =  cloud_dat_vx[kk * ndim_file * ndim_file + jj * ndim_file + ii];
            if(rho < 0) ath_error("Found rho < 0 in grid file!");
            rho = rho * (drat - 1.0);
            dye = rho;
            rho += 1.0; // also add wind density
            vx  = vx * vflow;
          }
        } else { // do not read grid from file
          if (r < r_cloud) {
#ifdef FLOW_PROFILE
            if(iprint == 0) {
              ath_pout(-1, "[init_problem] t_cc = %g\tt_cool,cl = %g\n",
                       sqrt(drat) * r_cloud / vx,
                       tcool(rho, flow_profile_pressure(x1, x2, x3) / (rho * drat)));
              iprint = 1;
            }
#endif
            vx   = -(x2/r_cloud) * v_rot;
            vy   =  (x1/r_cloud) * v_rot;

            rho  *= drat;
#if (NSCALARS > 0)
            dye = drat ; //1.0;
#endif
          }
          if (dr > 0.0){
#ifdef FLOW_PROFILE
            ath_error("dr > 0 with flow profile not (yet) supported.\n");
#endif
            rho = (1.0 + drat*0.5*(1.0+tanh((r_cloud-r)/(dr*r_cloud))));
            vx = vflow*0.5*(1.0+tanh((-(1.0+dr)*r_cloud+r)/(dr*r_cloud)))/rho;
            vx   += -(x2/r_cloud) * v_rot;
            vy   =  (x1/r_cloud) * v_rot;
#if (NSCALARS > 0)
            // This line means that the dye does not follow the density in the boundary (dr) region
            if(r < r_cloud){ 
              dye = rho;
            }
#endif
          }
        }

        /* write values to the grid */
        pGrid->U[k][j][i].d = rho;
        pGrid->U[k][j][i].M1 = rho * vx;
        pGrid->U[k][j][i].M2 = rho * vy;
        pGrid->U[k][j][i].M3 = rho * vz;

        /* if (r < r_cloud && v_turb > 0.0) {
          pGrid->U[k][j][i].M1 += rho * RandomNormal(0.0, v_turb);
          pGrid->U[k][j][i].M2 += rho * RandomNormal(0.0, v_turb);
          if (pGrid->Nx[2] > 1)
            pGrid->U[k][j][i].M3 += rho * RandomNormal(0.0, v_turb);
        }
        */

        // Defining the pressure implicitly through the "+ 1.0" --> P_init = Gamma - 1
#ifndef ISOTHERMAL
#ifdef FLOW_PROFILE
        pGrid->U[k][j][i].E = flow_profile_pressure(x1, x2, x3) / Gamma_1;
#else
        pGrid->U[k][j][i].E = 1.0 + dp / Gamma_1 ;
#endif /* FLOW_PROFILE */
        pGrid->U[k][j][i].E += 0.5 * rho * (SQR(vx)+SQR(vy)+SQR(vz));
#endif  /* ISOTHERMAL */

#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = dye;
#endif

#ifdef ENERGY_COOLING
        pGrid->U[k][j][i].Erad = 0;
#endif

        /* printf("1337 %g %g %g %g %g\n", x1, pGrid->U[k][j][i].d, vx, vy, vz); */


      }
    }
  } /* end grid loops */

  if(fp_rho != NULL) {
    fclose(fp_rho);
    fclose(fp_vx);
    free(cloud_dat_rho);
    free(cloud_dat_vx);
  }

  if(tangled){
    A = (Real3Vect***) calloc_3d_array(nx3, nx2, nx1, sizeof(Real3Vect));


    for (k=0; k<nx3; k++) {
      for (j=0; j<nx2; j++) {
        for (i=0; i<nx1; i++) {
          A[k][j][i].x3 = A[k][j][i].x2 = A[k][j][i].x1 = 0.0;
        }
      }
    }


    nterms = par_getd_def("problem", "nterms",10);
    alpha = par_getd_def("problem","alpha",50.0);
    for (i=0; i<nterms; i++) {
      if (myID_Comm_world == 0) {
        phi = randomreal2(0.0, 2.0*PI);
        theta   = acos(2.0*randomreal2(0.0, 1.0)-1.0);
        beta  = randomreal2(0.0, 2.0*PI);

        amp   = RandomNormal2(1.0, 0.25);

      } else {
        theta = phi = beta = amp = alpha = 0.0;
      }

#ifdef MPI_PARALLEL
      my_scal[0] = theta;
      my_scal[1] = phi;
      my_scal[2] = beta;
      my_scal[3] = amp;
      my_scal[4] = alpha;


      ierr = MPI_Allreduce(&my_scal, &scal, 5, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
      if (ierr)
        ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

      theta = scal[0];
      phi   = scal[1];
      beta  = scal[2];
      amp   = scal[3];
      alpha = scal[4];
#endif

      add_term(A, pGrid, theta, phi, alpha, beta, amp);
    } /* end loop over nterms */ 

    Bin = sqrt(2.0 * (Gamma_1 + dp) / betain);

    for (k=0; k<nx3; k++) {
      for (j=0; j<nx2; j++) {
        for (i=0; i<nx1; i++) {
          A[k][j][i].x1 *= Bin / sqrt(nterms);
          A[k][j][i].x2 *= Bin / sqrt(nterms);
          A[k][j][i].x3 *= Bin / sqrt(nterms);

          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
          r = sqrt(x1*x1 + x2*x2 + x3*x3);
          if (r > r_cloud)
            A[k][j][i].x1 = A[k][j][i].x2 = A[k][j][i].x3 = 0.0;
        }
      }
    }


  } /*end of if(tangled) */


  /* take curl here */
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        pGrid->B1i[k][j][i] = 0.0;
        if(tangled){
          pGrid->B1i[k][j][i] = (A[k][j+1][i].x3 - A[k][j][i].x3)/pGrid->dx2 -
            (A[k+1][j][i].x2 - A[k][j][i].x2)/pGrid->dx3;
        }
      }
    }
  }

  Bout_z = sqrt(2.0 * (Gamma_1 + dp) / betaout_z);
  Bout_y = sqrt(2.0 * (Gamma_1 + dp) / betaout_y);

  ju = (pGrid->Nx[1] > 1) ? je+1 : je;  //so we don't go beyond array in 2d 
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][j][i] = 0.0;
        if(tangled){
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
          bscale = 0.5+0.5*tanh((-r_cloud*3.-x1)*3./r_cloud);
          pGrid->B2i[k][j][i] = (A[k+1][j][i].x1 - A[k][j][i].x1)/pGrid->dx3 -
            (A[k][j][i+1].x3 - A[k][j][i].x3)/pGrid->dx1 + bscale * Bout_y;
          //printf("1234.1234 %g %g\n", x1, pGrid->B2i[k][j][i]);
        }
        else{
          pGrid->B2i[k][j][i] = Bout_y;
        }
      }
    }
  }

  ku = (pGrid->Nx[2] > 1) ? ke+1 : ke; //so we don't go beyond array in 2d
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[k][j][i] = 0.0;
        if(tangled){
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
           bscale = 0.5+0.5*tanh((-r_cloud*3.-x1)*3./r_cloud);
          pGrid->B3i[k][j][i] = (A[k][j][i+1].x2 - A[k][j][i].x2)/pGrid->dx1 -
            (A[k][j+1][i].x1 - A[k][j][i].x1)/pGrid->dx2 + bscale * Bout_z;
        }
        else{
          pGrid->B3i[k][j][i] = Bout_z;
        }
      }
    }
  }


#endif
   if(tangled){
     free_3d_array((void***) A);
   }
  /* cell-centered magnetic field */
  /*   derive this from interface field to be internally consistent
       with athena */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef MHD
        pGrid->U[k][j][i].B1c =
          0.5 * (pGrid->B1i[k][j][i] + pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = pGrid->B2i[k][j][i];
        pGrid->U[k][j][i].B3c = pGrid->B3i[k][j][i];

        if (pGrid->Nx[1] > 1)
          pGrid->U[k][j][i].B2c =
            0.5 * (pGrid->B2i[k][j][i] + pGrid->B2i[k][j+1][i]);
        if (pGrid->Nx[2] > 1)
          pGrid->U[k][j][i].B3c =
            0.5 * (pGrid->B3i[k][j][i] + pGrid->B3i[k+1][j][i]);

#ifndef ISOTHERMAL
        /* add magnetic energy to the total energy */
        pGrid->U[k][j][i].E +=
          0.5 * (SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B2c)+
                 SQR(pGrid->U[k][j][i].B3c));
#endif  /* ISOTHERMAL */
#endif  /* MHD */
      }
    }
  }

#ifdef MHD


  /* -- Some checks -- */
  Brms = Bmax = Bmax_cloud = 0.0;
  ncells = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        Bs = (SQR(pGrid->U[k][j][i].B1c) +
              SQR(pGrid->U[k][j][i].B2c) +
              SQR(pGrid->U[k][j][i].B3c));

        if(Bs > 1.e-6){
          Brms += Bs;
          ncells += 1.0;
        }
        Bmax = MAX(Bmax, Bs);

        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
        r = sqrt(x1*x1 + x2*x2 + x3*x3);
        if (r < r_cloud)
          Bmax_cloud = MAX(Bmax_cloud, Bs);
      }
    }
  }

#ifdef MPI_PARALLEL
  my_scal[0] = Brms;
  my_scal[1] = ncells;

  ierr = MPI_Allreduce(&my_scal, &scal, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

  Brms   = scal[0];
  ncells = scal[1];

  my_scal[0] = Bmax;
  my_scal[1] = Bmax_cloud;

  ierr = MPI_Allreduce(&my_scal, &scal, 2, MPI_RL, MPI_MAX, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[problem]: MPI_Allreduce returned error %d\n", ierr);

  Bmax   = scal[0];
  Bmax_cloud   = scal[1];
#endif

  Brms = sqrt(Brms/ncells);
  Bmax = sqrt(Bmax);
  Bmax_cloud = sqrt(Bmax_cloud);

  ath_pout(0, "Brms = %f, Bmax = %f, Bmax_cloud = %f, "
           "beta_rms = %f, beta_max = %f, beta_cloud_max = %f\n",
           Brms, Bmax, Bmax_cloud,
           2 * (Gamma_1 + dp) / SQR(Brms),
           2 * (Gamma_1 + dp) / SQR(Bmax), 2 * (Gamma_1 + dp) / SQR(Bmax_cloud));
  /* ath_error("Brms = %f\tBmax = %f\n", Brms, Bmax); */


#endif //MHD


#ifdef MHD
  check_div_b(pGrid);
#endif //MHD

  if (pDomain->Disp[0] == 0)
    bvals_mhd_fun(pDomain, left_x1,  bc_ix1);
  /* if (pDomain->Disp[1] == 0) */
  /*   bvals_mhd_fun(pDomain, left_x2,  bc_ix2); */
  /* if (pDomain->Disp[2] == 0) */
  /*   bvals_mhd_fun(pDomain, left_x3,  bc_ix3); */
  if (pDomain->MaxX[0] == pDomain->RootMaxX[0])
    bvals_mhd_fun(pDomain, right_x1, bc_ox1);
  /* if (pDomain->MaxX[1] == pDomain->RootMaxX[1]) */
  /*   bvals_mhd_fun(pDomain, right_x2, bc_ox2); */
  /* if (pDomain->MaxX[2] == pDomain->RootMaxX[2]) */
  /*   bvals_mhd_fun(pDomain, right_x3, bc_ox3); */


  /* seed a perturbation */
   for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        fact = -1.0;
        while (fabs(fact) > 0.03)
          fact = (RandomNormal(0.0, 0.01));
        if(fabs(fact) < .03)
          pGrid->U[k][j][i].d *= (1.0+fact);
      }
    }
   }

#ifdef PARTICLES
   init_particles(pDomain);
#endif

  return;
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
#ifdef FOLLOW_CLOUD
  fwrite(&x_shift, sizeof(Real), 1, fp);
  fwrite(&vflow,   sizeof(Real), 1, fp);
#endif

  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  for (nl=0; nl<(pM->NLevels); nl++) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
      if (pM->Domain[nl][nd].Disp[0] == 0)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x1,  bc_ix1);
       /* if (pM->Domain[nl][nd].Disp[1] == 0) */
       /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  bc_ix2); */
       /* if (pM->Domain[nl][nd].Disp[2] == 0) */
       /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3,  bc_ix3); */
       if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0])
         bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, bc_ox1);
       /* if (pM->Domain[nl][nd].MaxX[1] == pM->Domain[nl][nd].RootMaxX[1]) */
       /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, bc_ox2); */
       /* if (pM->Domain[nl][nd].MaxX[2] == pM->Domain[nl][nd].RootMaxX[2]) */
       /*   bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, bc_ox3); */

    }
  }

  drat  = par_getd("problem", "drat");
  vflow = par_getd("problem", "vflow");
  vflow0 = vflow;

  tfloor = par_getd_def("problem", "tfloor", 1.e-2/drat);
  tceil = par_getd_def("problem", "tceil", 100.);
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-2);
  d_MIN = par_getd_def("problem", "d_MIN", 1.e-4);
  dtmin = par_getd_def("problem", "dtmin", 1.e-7);

#ifdef MHD
  betain = par_getd("problem", "betain");
  betaout_y = par_getd_def("problem", "betaout_y", 1e20);
  betaout_z = par_getd("problem", "betaout_z");
  betafloor = par_getd_def("problem", "betafloor", 3.e-3);
#endif

  dp = par_getd_def("problem", "dp", 0.0);
  tnotcool = par_getd_def("problem", "tnotcool", -1.0);
  tfloor_cooling = par_getd_def("problem", "tfloor_cooling", (Gamma_1 + dp) / drat);

#ifdef VISCOSITY
  nu   = par_getd("problem","nu");
  NuFun_i = NULL;
  NuFun_a = nu_fun;
#endif

#ifdef PARTICLES
  parnumproc = (int)(par_geti("particle","parnumproc"));
#endif


#ifdef ENERGY_COOLING
  init_cooling();
#endif

#ifdef INSTANTCOOL
  //  CoolingFunc = instant_cool;
#endif

  // Re-enroll hst dumps after restart
  dump_history_enroll(hst_m13, "m13");
  dump_history_enroll(hst_m110, "m110");
  dump_history_enroll(hst_Mx13, "Mx13");

  dump_history_enroll(hst_Erad, "Erad");


#ifdef FOLLOW_CLOUD
  dump_history_enroll_alt(hst_xshift, "x_shift");
  dump_history_enroll_alt(hst_vflow,  "v_flow");
#endif

#if (NSCALARS > 0)
  dump_history_enroll(hst_c,    "<c>");
  dump_history_enroll(hst_c_sq, "<c^2>");

  dump_history_enroll(hst_cE,  "<c * E>");
  dump_history_enroll(hst_cx1, "<c * x1>");

  dump_history_enroll(hst_cvx, "<c * Vx>");
  dump_history_enroll(hst_cvy, "<c * Vy>");
  dump_history_enroll(hst_cvz, "<c * Vz>");

  dump_history_enroll(hst_cvx_sq, "<(c * Vx)^2>");
  dump_history_enroll(hst_cvy_sq, "<(c * Vy)^2>");
  dump_history_enroll(hst_cvz_sq, "<(c * Vz)^2>");
#ifdef MHD
  dump_history_enroll(hst_cBx, "<c * Bx>");
  dump_history_enroll(hst_cBy, "<c * By>");
  dump_history_enroll(hst_cBz, "<c * Bz>");
#endif /* MHD */
  dump_history_enroll(hst_Sdye, "dye entropy");
#endif  /* NSCALARS */

#ifdef ENERGY_COOLING
  dump_history_enroll(hst_cstcool, "cs*tcool");
#endif


  /* DANGER: make sure the order here matches the order in write_restart() */
#ifdef FOLLOW_CLOUD
  fread(&x_shift, sizeof(Real), 1, fp);
  fread(&vflow,   sizeof(Real), 1, fp);
#endif

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef MHD
static void check_div_b(GridS *pGrid)
{
  int i,j,k;
  int is,js,ie,je,ks,ke;
  int n;
  Real divb[5],divcell,x,y,z;
#ifdef MPI_PARALLEL
  int ierr;
  Real my_divb[5];
#endif
  int three_d = (pGrid->Nx[2] > 1) ? 1 : 0;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;


  divb[0] = 0.0;
  n = 0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        divcell = (pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
          + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2;
        if (three_d)
          divcell+= (pGrid->B3i[k+1][j][i]-pGrid->B3i[k][j][i])/pGrid->dx3;
        divb[0] += fabs(divcell);
#if (NSCALARS > 1)
        pGrid->U[k][j][i].s[1] = divcell;
#endif

        /*if(fabs(divcell) > 1.e-13){
          cc_pos(pGrid,i,j,k,&x,&y,&z);

          printf("bdv %e %d %d %d  %e  %e    %e %e %e\n", pGrid->time, i,j,k,  divcell, 1./pGrid->dx1, x,y,z);
        }
        */
        n++;
      }
    }
  }
  //divb[0] /= n;

  /* check divB along ix1: */
  divb[1] = 0.0;
  n=0;
  k = ks;
  i = is;
  for (j=js; j<=je; j++) {
    divb[1] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[1] /= n;

  /* check divB along ox1: */
  divb[2] = 0.0;
  n=0;
  k = ks;
  i = ie;
  for (j=js; j<=je; j++) {
    divb[2] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[2] /= n;


  /* check divB along ix2: */
  divb[3] = 0.0;
  n=0;
  k = ks;
  j = js;
  for (i=is; i<=ie; i++) {
    divb[3] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[3] /= n;


  /* check divB along ox2: */
  divb[4] = 0.0;
  n=0;
  k = ks;
  j = je;
  for (i=is; i<=ie; i++) {
    divb[4] +=
      fabs((pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1
           + (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2);
    n++;
  }
  divb[4] /= n;


#ifdef MPI_PARALLEL
  for (i=0; i<=4; i++)
    my_divb[i] = divb[i];
  ierr = MPI_Allreduce(&my_divb, &divb, 5, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[check_div_b]: MPI_Allreduce returned error %d\n", ierr);
#endif



  if (three_d) {
    ath_pout(0, "divb = %e\n", divb[0]);
  } else {
    ath_pout(0, "divb = %e\t%e\t%e\t%e\t%e\n",
             divb[0], divb[1], divb[2], divb[3], divb[4]);
  }

  return;
}
#endif  /* MHD */

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
#ifdef FOLLOW_CLOUD
  Real dvx, newdvx, expt;
#endif


#ifdef FOLLOW_CLOUD
  dvx = cloud_mass_weighted_velocity(pM);
  if(fabs(dvx) > 100.01 || dvx < 0.0){ // TODO: some maximum shift is defined here...
    ath_pout(0,"[bad dvx:] %0.15e setting to 0.\n",dvx);
    dvx = 0.0;
  }

  /*  
  if(dvx > 0.01) {
    dvx = 0.01;
    ath_pout(0,"[bad dvx:] %0.15e setting to 0.01\n",dvx);
  }
  */

  if(dvx > 0.0){
    expt = floor(log10(dvx));
    newdvx = dvx / pow(10, expt);
    newdvx = floor(newdvx * 1.0e4) / 1.0e4;
    newdvx = newdvx * pow(10.0, expt);
    dvx = newdvx;
  }

#ifndef FLOW_PROFILE
  if(vflow - dvx < 0.00){ // does not allow vflow < 0
    dvx = vflow;
  }
#endif /* not FLOW_PROFILE */
  ath_pout(0,"[dvx:]  %0.15e [vflow:] %.15e\n",dvx,vflow);
  vflow -= dvx;
#endif /* FOLLOW_CLOUD */

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef FOLLOW_CLOUD
        boost_frame(&(pM->Domain[nl][nd]), dvx);
#endif
      }
    }
  }

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef MHD
        check_div_b(pM->Domain[nl][nd].Grid);
#endif
      }
    }
  }

  for (nl=0; nl<=(pM->NLevels)-1; nl++) {
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef ENERGY_COOLING
        integrate_cooling(pM->Domain[nl][nd].Grid);
#endif
      }
    }
  }
#ifdef ENERGY_HEATING
  radiate_energy(pM);
#endif

  if (pM->dt < dtmin){
    data_output(pM,1);
    ath_error("dt too small\n");

  }

#ifdef PARTICLES
  add_new_particles(pM);
#endif

  return;
}

void Userwork_after_loop(MeshS *pM)
{
#ifdef FOLLOW_CLOUD
  ath_pout(0, "[follow_cloud]: shifted the domain by a total amount %e\n",
           x_shift);
#endif

  return;
}



/*==============================================================================
 * PHYSICS FUNCTIONS:
 * boost_frame()         - boost simulation frame by a velocity increment
 * report_nans()         - apply a ceiling and floor to the temperature
 * cloud_velocity()      - find the mass-weighted velocity of the cloud.
 * nu_fun()              - kinematic viscosity (i.e., cm^2/s)
 * cooling_func()        - cooling function for the energy equation
 *----------------------------------------------------------------------------*/
#ifdef FOLLOW_CLOUD
static void boost_frame(DomainS *pDomain, Real dvx)
{
  int i, j, k;
  int is,ie,js,je,ks,ke;

  Real d;

  GridS *pGrid = pDomain->Grid;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        d = pGrid->U[k][j][i].d;

#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += 0.5 * d * SQR(dvx);
        pGrid->U[k][j][i].E -= dvx * pGrid->U[k][j][i].M1;
#endif  /* ISOTHERMAL */
        pGrid->U[k][j][i].M1 -= dvx * d;

      }
    }
  }

  //  x_shift -= dvx * pDomain->Grid->dt;
  x_shift -= vflow * pDomain->Grid->dt;

  return;
}
#endif  /* FOLLOW_CLOUD */


#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix)
{
#ifndef ISOTHERMAL
  int i, j, k;
  int is,ie,js,je,ks,ke;
  Real x1, x2, x3;
  int V=0; //verbose off = 0
  int NO = 10;
  Real KE, rho, press, temp;
  int nanpress=0, nanrho=0, nanv=0, nnan;   /* nan count */
  int npress=0,   nrho=0,   nv=0,   nfloor; /* floor count */
#ifdef MHD
  Real ME;
  int nanmag=0;
  int nmag=0;
#endif  /* MHD */
  Real beta;
  Real scal[8];
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
#ifdef MHD
        ME = (SQR(pGrid->U[k][j][i].B1c) +
              SQR(pGrid->U[k][j][i].B2c) +
              SQR(pGrid->U[k][j][i].B3c)) * 0.5;
        press -= ME;
#endif  /* MHD */

        press *= Gamma_1;
        temp = press / rho;
#ifdef MHD
        beta = press / ME;
#else
        beta = fabs(200.*betafloor);
#endif  /* MHD */

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

#ifdef MHD
        if (ME != ME) {
          nanmag++;
          /* TODO: apply a fix to B? */
        } else if (beta < betafloor && betafloor > 0.0) {
          nmag++;
          if(V && nmag < NO) printf("bad mag %e R %e  %e %e %e  %d %d %d %e %e %e %e %e\n",beta, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta, ME);
          if(fix){
	    // Old fix --> via new T
            rho  = MAX(rho, sqrt(betafloor * ME));
            temp = betafloor * ME / rho;
            temp = MAX(temp,tfloor);

	    /*
	    // New fix --> via new B 
	    temp = press / rho; // Otherwise beta sometimes < 0
	    beta = press / ME; 
	    if(beta < 0)
	      ath_error("[floor_beta]: beta < 0 (%g), rho = %g, temp = %g, press = %g, ME = %g\n", beta, rho, temp, press, ME);
	    pGrid->U[k][j][i].B1c *= sqrt(beta / betafloor);
	    pGrid->U[k][j][i].B2c *= sqrt(beta / betafloor);
	    pGrid->U[k][j][i].B3c *= sqrt(beta / betafloor);
	    ME *= beta / betafloor;
	    */

            if(V && nmag < NO) printf("bad magf %e R %e  %e %e %e  %d %d %d %e %e %e %e %e\n",beta, 1.0/pGrid->dx1, x1, x2, x3, i, j, k,rho, press, temp, beta,ME);
          }
        }
#endif  /* MHD */

        /* write values back to the grid */
        /* TODO: what about B??? */
        if(fix) {
          pGrid->U[k][j][i].d  = rho;
          KE = (SQR(pGrid->U[k][j][i].M1) +
                SQR(pGrid->U[k][j][i].M2) +
                SQR(pGrid->U[k][j][i].M3)) /
            (2.0 * rho);

          pGrid->U[k][j][i].E = temp*rho/Gamma_1 + KE;
#ifdef MHD
          pGrid->U[k][j][i].E += ME;
#endif  /* MHD */
        }
      }
    }
  }

  /* synchronize over grids */
#ifdef MPI_PARALLEL
  my_scal[0] = nanpress;
  my_scal[1] = nanrho;
  my_scal[2] = nanv;
#ifdef MHD
  my_scal[3] = nanmag;
#endif  /* MHD */
  my_scal[4] = npress;
  my_scal[5] = nrho;
  my_scal[6] = nv;
#ifdef MHD
  my_scal[7] = nmag;
#endif  /* MHD */

  ierr = MPI_Allreduce(&my_scal, &scal, 8, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[report_nans]: MPI_Allreduce returned error %d\n", ierr);

  nanpress = scal[0];
  nanrho   = scal[1];
  nanv     = scal[2];
#ifdef MHD
  nanmag   = scal[3];
#endif  /* MHD */
  npress   = scal[4];
  nrho     = scal[5];
  nv       = scal[6];
#ifdef MHD
  nmag     = scal[7];
#endif  /* MHD */
#endif  /* MPI_PARALLEL */


  /* sum up the # of bad cells and report */
  nnan = nanpress+nanrho+nanv;

#ifdef MHD
  nnan += nanmag;
#endif  /* MHD */

  /* sum up the # of floored cells and report */
  nfloor = npress+nrho+nv;

#ifdef MHD
  nfloor += nmag;
#endif  /* MHD */

  // if (fix == 0) {
#ifdef MHD
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v, %d beta.\n",
      nfloor, npress, nrho, nv, nmag);
#else
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v.\n",
      nfloor, npress, nrho, nv);
#endif  /* MHD */
    // }


  //  if ((nnan > 0 || nfloor -nmag > 30) && fix == 0) {
    if (nnan > 0 ){// && fix == 0) {
#ifdef MHD
    ath_pout(0, "[report_nans]: found %d nan cells: %d P, %d d, %d v, %d B.\n",
             nnan, nanpress, nanrho, nanv, nanmag);
#else
    ath_pout(0, "[report_nans]: found %d nan cells: %d P, %d d, %d v.\n",
             nnan, nanpress, nanrho, nanv);
#endif  /* MHD */


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

    if (nfloor > 50000)
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
  int k, n=nfit_cool-1;
  Real term;
  const Real mu = 0.62, mu_e = 1.17;

  /* convert T in the cooling function from keV to code units */
  for (k=0; k<=n; k++)
    sdT[k] /= (8.197 * mu);

  if(tfloor_cooling < sdT[0])
    ath_error("Cooling floor is smaller than first entry of cooling function.");

  /* populate Yk following equation A6 in Townsend (2009) */
  Yk[n] = 0.0;
  for (k=n-1; k>=0; k--){
    term = (sdL[n]/sdL[k]) * (sdT[k]/sdT[n]);

    if (sdexpt[k] == 1.0)
      term *= log(sdT[k]/sdT[k+1]);
    else
      term *= ((1.0 - pow(sdT[k]/sdT[k+1], sdexpt[k]-1.0)) / (1.0-sdexpt[k]));

    Yk[k] = Yk[k+1] - term;

    if(isnan(Yk[k]))
      ath_error("Error initializing cooling. nan in Yk[%d]", k);
  }

  return;
}

/* piecewise power-law fit to the cooling curve with temperature in
   keV and L in 1e-23 erg cm^3 / s */
static Real sdLambda(const Real T)
{
  int k, n=nfit_cool-1;

  /* first find the temperature bin */
  for(k=n; k>=0; k--){
    if (T >= sdT[k])
      break;
  }

  /* piecewise power-law; see equation A4 of Townsend (2009) */
  /* (factor of 1.311e-5 takes lambda from units of 1e-23 erg cm^3 /s
     to code units.) */
  return (1.311e-5 * sdL[k] * pow(T/sdT[k], sdexpt[k]));
}

static Real tcool(const Real d, const Real T)
{
  const Real mu = 0.62, mu_e = 1.17;

  /* equation 13 of Townsend (2009) */
  return (SQR(mu_e) * T) / (Gamma_1 * d * sdLambda(T));
}

/* see sdLambda() or equation A1 of Townsend (2009) for the
   definition */
static Real Y(const Real T)
{
  int k, n=nfit_cool-1;
  Real term;

  /* first find the temperature bin */
  for(k=n; k>=0; k--){
    if (T >= sdT[k])
      break;
  }

  /* calculate Y using equation A5 in Townsend (2009) */
  term = (sdL[n]/sdL[k]) * (sdT[k]/sdT[n]);

  if (sdexpt[k] == 1.0)
    term *= log(sdT[k]/T);
  else
    term *= ((1.0 - pow(sdT[k]/T, sdexpt[k]-1.0)) / (1.0-sdexpt[k]));

  return (Yk[k] + term);
}

static Real Yinv(const Real Y1)
{
  int k, n=nfit_cool-1;
  Real term;

  /* find the bin i in which the final temperature will be */
  for(k=n; k>=0; k--){
    if (Y(sdT[k]) >= Y1)
      break;
  }


  /* calculate Yinv using equation A7 in Townsend (2009) */
  term = (sdL[k]/sdL[n]) * (sdT[n]/sdT[k]);
  term *= (Y1 - Yk[k]);

  if (sdexpt[k] == 1.0)
    term = exp(-1.0*term);
  else{
    term = pow(1.0 - (1.0-sdexpt[k])*term,
               1.0/(1.0-sdexpt[k]));
  }

  return (sdT[k] * term);
}

static Real newtemp_townsend(const Real d, const Real T, const Real dt_hydro)
{
  Real term1, Tref;
  int n=nfit_cool-1;

  Tref = sdT[n];

  term1 = (T/Tref) * (sdLambda(Tref)/sdLambda(T)) * (dt_hydro/tcool(d, T));

  return Yinv(Y(T) + term1);
}

static void integrate_cooling(GridS *pG)
{
  int i, j, k;
  int is, ie, js, je, ks, ke;

  PrimS W;
  ConsS U;
  // Changed this for tfloor!!
  Real temp;

  /* ath_pout(0, "integrating cooling using Townsend (2009) algorithm.\n"); */

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        W = Cons_to_Prim(&(pG->U[k][j][i]));

        /* find temp in keV */
        temp = W.P/W.d;

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
        pG->U[k][j][i].Erad += (pG->U[k][j][i].E - U.E);

        pG->U[k][j][i].E = U.E;

      }
    }
  }

  return;

}


#ifdef ENERGY_HEATING
/*
  Radiate cooled energy over whole domain
 */
static void radiate_energy(MeshS *pM) {
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
              pG->U[k][j][i].Erad = 0.0;
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
  int i, npts=100;
  Real logt, temp, tc, logdt, dt;
  Real err;

  FILE *outfile;

  outfile = fopen("lambda.dat", "w");
  for(i=0; i<npts; i++){
    logt = log(1.0e-4) + (log(5.0)-log(1.0e-4))*((double) i/(npts-1));
    temp = exp(logt);

    fprintf(outfile, "%e\t%e\n", temp, sdLambda(temp));
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

  temp = 3.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-3kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 1.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-1kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 0.3;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.3kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }

  temp = 0.1;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.1kev.dat", "w");
  for(i=0; i<npts; i++){
    logdt = log(0.1) + (log(2.0)-log(0.1))*((double) i / (npts-1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt/tc, newtemp_townsend(1.0, temp, dt));
  }


  ath_error("check cooling stuff.\n");

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



/* Return the center-of-mass velocity of the cloud */
/*   the scalar s obeys the same equation as density, so we weight
 *   the average as s*v. */
/*   for now, compute the average using the level 1 domain (i.e, the
 *   first refined level) */
#ifdef FOLLOW_CLOUD
static Real cloud_mass_weighted_velocity(MeshS *pM)
{
  GridS *pG;
  int i, j, k, is, ie, js, je, ks, ke;
  int nl, nd;

  Real s, d, scal[2], tmp, x1, x2, x3;
#ifdef MPI_PARALLEL
  Real my_scal[2];
  int ierr;
#endif

  /* do the integral over level-1 domains, if they exist */
  nl = (pM->NLevels > 1) ? 1 : 0;

  scal[0] = scal[1] = 0.0;
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {

      pG = pM->Domain[nl][nd].Grid;
      is = pG->is;  ie = pG->ie;
      js = pG->js;  je = pG->je;
      ks = pG->ks;  ke = pG->ke;

      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            d = pG->U[k][j][i].d;
#if(NSCALARS > 0)
            s = pG->U[k][j][i].s[0];
#endif
            cc_pos(pG,i,j,k,&x1,&x2,&x3);
            tmp = s * pG->U[k][j][i].M1 / d;
            if (tmp == tmp && x1 < 0.0) {

              scal[0] += tmp*d;
              scal[1] += s*d;
            }

          }
        }
      }

    }
  }

#ifdef MPI_PARALLEL
  my_scal[0] = scal[0];
  my_scal[1] = scal[1];

  ierr = MPI_Allreduce(&my_scal, &scal, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (ierr)
    ath_error("[cloud_velocity]: MPI_Allreduce returned error %d\n", ierr);
#endif

  return scal[0] / scal[1];
}
#endif  /* FOLLOW_CLOUD */


#ifdef VISCOSITY
static Real nu_fun(const Real d, const Real T,
                   const Real x1, const Real x2, const Real x3)
{
  Real newnu;
  return nu;
  /* newnu = nu*pow(T*1.5,2.5); */
  /* if(newnu != newnu || newnu < 0.0) */
  /*   newnu = TINY_NUMBER; */
  /* if(newnu > 3.*nu) */
  /*   newnu = 3.*nu; */

  /* return newnu; */
}
#endif  /* VISCOSITY */




/*==============================================================================
 * HISTORY OUTPUTS:
 *

 *----------------------------------------------------------------------------*/

static Real _hst_mcut(GridS *pG, int i, int j, int k, const Real frac)
{
  if(pG->U[k][j][i].d < frac * drat)
    return 0;
  return pG->U[k][j][i].d;
}


static Real hst_m13(GridS *pG, int i, int j, int k)
{
  return _hst_mcut(pG, i, j, k, 1/3.);
}

static Real hst_m110(GridS *pG, int i, int j, int k)
{
  return _hst_mcut(pG, i, j, k, 0.1);
}


static Real hst_Mx13(GridS *pG, int i, int j, int k)
{
  if(pG->U[k][j][i].d < drat / 3.)
    return 0;
  return pG->U[k][j][i].M1;
}

static Real hst_Erad(GridS *pG, int i, int j, int k)
{
  return pG->U[k][j][i].Erad;
}


#ifdef FOLLOW_CLOUD
static Real hst_xshift(GridS *pG, int i, int j, int k)
{
  return x_shift;
}
#endif

#ifdef FOLLOW_CLOUD
static Real hst_vflow(GridS *pG, int i, int j, int k)
{
  return vflow;
}
#endif


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

#ifdef MHD
static Real hst_cBx(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].B1c);
}

static Real hst_cBy(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].B2c);
}

static Real hst_cBz(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].B3c);
}
#endif  /* MHD */

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


/*==============================================================================
 * BOUNDARY CONDITIONS:
 *
 *----------------------------------------------------------------------------*/

static void bc_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real vx, vy, vz, rho, cx1;
  Real x1, x2, x3, r;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  vx = vy = vz = 0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        /*
         This line was uncommented before...which seems to change the rest of
         the function. Seems better that way, though...?
         *but* caused MHD error!
        */
        // pGrid->U[k][j][is-i] = pGrid->U[k][j][is];

#if (NSCALARS > 0)
        pGrid->U[k][j][is - i].s[0] = 0.0;
#endif

#ifdef FLOW_PROFILE
        cc_pos(pGrid,is - i,j,k,&x1,&x2,&x3);
        cx1 = x1 + x_shift; // x_shift > 0
        vx = flow_profile_velocity_x(cx1, x2, x3) + vflow; // v_flow < 0
        vy = flow_profile_velocity_y(cx1, x2, x3);
        vz = flow_profile_velocity_z(cx1, x2, x3);
        rho = flow_profile_density(cx1, x2, x3);
        /* if((k == ks) && (j == js)) */
        /*   printf("1338 %g %g %g %g\n", cx1, rho, vx - vflow, vflow); */
#else
        vx = vflow;
        rho = 1.0;
#endif
        pGrid->U[k][j][is-i].d  = rho;
        pGrid->U[k][j][is-i].M1 = rho * vx;
        pGrid->U[k][j][is-i].M2 = rho * vy;
        pGrid->U[k][j][is-i].M3 = rho * vz;

#ifndef ISOTHERMAL
#ifdef FLOW_PROFILE
        pGrid->U[k][j][is-i].E = flow_profile_pressure(cx1, x2, x3) / Gamma_1;
#else
        pGrid->U[k][j][is-i].E = 1.0 + dp / Gamma_1 ;
#endif /* FLOW_PROFILE */
        pGrid->U[k][j][is-i].E += 0.5 * rho * (SQR(vx) + SQR(vy) + SQR(vz));
#endif  /* ISOTHERMAL */

#ifdef MHD
        pGrid->U[k][j][is-i].B1c = 0.0;
        pGrid->U[k][j][is-i].B2c = sqrt(2.0 * (Gamma_1 + dp) / betaout_y);
        pGrid->U[k][j][is-i].B3c = sqrt(2.0 * (Gamma_1 + dp) / betaout_z);
        if(i == 1)
          pGrid->U[k][j][is-i].B1c = 0.5*pGrid->B1i[k][j][is];

        pGrid->U[k][j][is-i].E  += 0.5*(SQR(pGrid->U[k][j][is-i].B1c)
                                        +SQR(pGrid->U[k][j][is-i].B2c)
                                        +SQR(pGrid->U[k][j][is-i].B3c));
#endif
      }
    }
  } // End loop over grid cells


#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<nghost; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x1 -= 0.5 * pGrid->dx1;
        r = sqrt(x1*x1+x2*x2);

        pGrid->B1i[k][j][i] = 0.0;
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=0; i<nghost; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        x2 -= 0.5 * pGrid->dx2;
        r = sqrt(x1*x1+x2*x2);

        pGrid->B2i[k][j][i] = sqrt(2.0 * (Gamma_1 + dp) / betaout_y);
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=0; i<nghost; i++) {
        pGrid->B3i[k][j][i] = sqrt(2.0 * (Gamma_1 + dp) / betaout_z);
      }
    }
  }
#endif /* MHD */

  return;
}

static void bc_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  int V=0;
  int NO=10;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
        if(pGrid->U[k][j][ie+i].M1 < 0.0){
	  if(V && (NO > 0)) {
	    printf("bc_ox1 %d %d %d %e\n",
		   i, j, k, pGrid->U[k][j][ie+i].M1);
	    NO--;
	  }
          pGrid->U[k][j][ie+i].E -= 0.5*SQR(pGrid->U[k][j][ie+i].M1)/pGrid->U[k][j][ie+i].d;
          pGrid->U[k][j][ie+i].M1 = 0.0;
        }
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}


void add_term(Real3Vect ***A, GridS *pG,
              Real theta, Real phi,
              Real alpha, Real beta, Real amp)
{
  int i, j ,k;
  Real phase;
  Real3Vect e1, e2, e3, kv, r;

  int is, ie, ks, ke, js, je;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  e1 = get_e1(theta, phi);
  e2 = get_e2(theta, phi);
  e3 = get_e3(theta, phi);

  /* e3 x e1 = e2 */
  /* think of e3 as x-axis, e1 as y-axis, e2 as z-axis */
  /* then use this: */
  /* A =   cos(a x)/a z^ +   sin(a x)/a y^ */
  /* B =   cos(a x)   z^ +   sin(a x)   y^ */
  /* j = a*cos(a x)   z^ + a*sin(a x)   y^ */

  /* think of k as x^, e1 as y^, e2 as z^ */
  /* k cross y = z, so this is consistent  */

  /* e3.x1 = 1.0;   e3.x2 = 0.0;   e3.x3 = 0.0; */
  /* e2.x1 = 0.0;   e2.x2 = 0.0;   e2.x3 = 1.0; */
  /* e1.x1 = 0.0;   e1.x2 = 1.0;   e1.x3 = 0.0; */
  /* alpha = 2.0 * PI; */

  kv.x1 = alpha * e3.x1;
  kv.x2 = alpha * e3.x2;
  kv.x3 = alpha * e3.x3;
  if(ie+nghost >= (ie-is)+1+2*nghost){
    ath_error("overrunning A array!\n");
  }


  for (k=0; k<=ke+nghost; k++) {
    for (j=0; j<=je+nghost; j++) {
      for (i=0; i<=ie+nghost; i++) {
        cc_pos(pG, i, j, k, &r.x1, &r.x2, &r.x3);
        r.x1 -= 0.5 * pG->dx1;
        r.x2 -= 0.5 * pG->dx2;
        r.x3 -= 0.5 * pG->dx3;
        phase = r.x1 * kv.x1 + r.x2 * kv.x2 + r.x3 * kv.x3 + beta;

        A[k][j][i].x1 += amp * (e2.x1 * cos(phase) + e1.x1 * sin(phase))/alpha;
        A[k][j][i].x2 += amp * (e2.x2 * cos(phase) + e1.x2 * sin(phase))/alpha;
        A[k][j][i].x3 += amp * (e2.x3 * cos(phase) + e1.x3 * sin(phase))/alpha;
      }
    }
  }

  return;
}




Real3Vect get_e1(Real theta, Real phi)
{
  Real3Vect ret;
  ret.x1 = -1.0 * sin(theta);
  ret.x2 = cos(theta) * cos(phi);
  ret.x3 = cos(theta) * sin(phi);

  return ret;
}

Real3Vect get_e2(__attribute__((unused))Real theta, Real phi)
{
  Real3Vect ret;
  ret.x1 = 0.0;
  ret.x2 = -1.0 * sin(phi);
  ret.x3 = cos(phi);

  return ret;
}

Real3Vect get_e3(Real theta, Real phi)
{
  Real3Vect ret;
  ret.x1 = cos(theta);
  ret.x2 = sin(theta) * cos(phi);
  ret.x3 = sin(theta) * sin(phi);

  return ret;
}

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





/*
Below follow the functions if one introduces particles
 */
#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
               Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2,
                        const Real x3, const Real v1, const Real v2, const Real v3)
{
  return;
}


void place_particle(GrainS *p, Real pos[]) {
  p->x1 = pos[0];
  p->x2 = pos[1];
  p->x3 = pos[2];

  p->pos = 1; /*!< position: 0: ghost; 1: grid; >=10: cross out/in; */

  if(SQR(p->x1) + SQR(p->x2) + SQR(p->x3) < r_cloud * r_cloud)
    p->v1 = 0.0;
  else
    p->v1 = vflow;
  p->v2 = 0.0;
  p->v3 = 0.0;

#ifdef MPI_PARALLEL
  p->init_id = myID_Comm_world;
#endif

}

/* Initialize particles

Place `npart` particles which are `parnumproc` for each processor except
proc 0 where it is `parnumproc` + `npart_cloud`
 */
void init_particles(DomainS *pDomain) {
  GridS *pGrid = pDomain->Grid;
  int i, j, fail;
  int n = 0;
  GrainS *p;
  int npart_cloud = (int)(par_geti("particle","parnumcloud"));
  // Define how close to the cloud the `parnumcloud` particles should be
  Real part_drcloud = par_getd_def("particle", "part_drcloud", 1.1);
  Real pos[3];
  int npart;
  parnumproc = (int)(par_geti("particle","parnumproc"));

  npart = parnumproc + npart_cloud;

  if (npart +2 > pGrid->arrsize)
    particle_realloc(pGrid, npart + 2);
  tstop0[0] = 0.0;

  for(i=0;i<npart;i++) {
    fail = 0;
    for(j = 0; j < 3; j++) {
      if(i >= parnumproc)
        // put inside or close to cloud
        pos[j] =  randomreal2(-r_cloud * part_drcloud, r_cloud * part_drcloud); 
      else
        pos[j] =  randomreal2(pGrid->MinX[j], pGrid->MaxX[j]); // uniformly on grid

      if((pos[j] > pGrid->MaxX[j]) || (pos[j] < pGrid->MinX[j]))
        fail = 1;
    }
    if(fail)
      continue;

    p = &(pGrid->particle[n]);
    if(i >= parnumproc) // cloud
      p->my_id = 1000000 + n - parnumproc + 1 + npart_cloud * myID_Comm_world;
    else
      p->my_id = n + myID_Comm_world * parnumproc;

    place_particle(p, pos);
    n += 1;
  }

  pGrid->nparticle = n;
  ath_pout(-1, "[init_particles] Initialized %d particles on processor %d (%g to %g).\n",
           pGrid->nparticle, myID_Comm_world, pGrid->MinX[0], pGrid->MaxX[0]);
}



/* Place new particles upstream */
void add_new_particles(MeshS *pM) {
  GridS *pGrid = pM->Domain[0][0].Grid;
  int i,j;
  GrainS *p;
  Real minX = pM->RootMinX[0];
  static long int cur_pid = 0;
  int id_offset = parnumproc * (myID_Comm_world + 4096) + 10000; // TODO: FIXME better
  Real pos[3];

  if(minX != pGrid->MinX[0])
    return; // don't feed new particles if not on upstream boundary

  // Put some new particles in grid
  for(i=pGrid->nparticle; i<parnumproc;i++) {
    pos[0] = minX;
    for(j = 1; j < 3; j++)
      pos[j] = randomreal2(pGrid->MinX[j], pGrid->MaxX[j]);
    p = &(pGrid->particle[i]);
    place_particle(p, pos);

    pGrid->nparticle++;
    p->my_id = id_offset + cur_pid; 
    cur_pid++;

    ath_pout(-1, "[add_new_particle] Added particle %ld (%.2g, %.2g) on processor %d\n",
             p->my_id, p->x2, p->x3, myID_Comm_world);
  }

}

#endif /* PARTICLES */
