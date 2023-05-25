#include "copyright.h"
#include <float.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "prob/math_functions.h"
#include <stdbool.h>

/* ----------- Define options here --------- */
#define FOLLOW_CLOUD     // Moving reference frame
#define REPORT_NANS      // Verbose
#define ENERGY_COOLING   // Cooling
#define ENERGY_HEATING 0 // 0 = no heating, 1 = heat what cooled, 2 = constant heating
//#define PRINTCOOLING
#define INPUTFILE         // Reads centerpositions of clouds from a file in the same directory

// #define INSTANTCOOL        //Test cooling
/* ----------------------------------------- */

/* ----------- Parallel options --------- */
#ifdef MPI_PARALLEL
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */
/* ----------------------------------------- */

/* ----------- Define variables and functions --------- */

/* Boundary terms*/
static void bc_ix1(GridS *pGrid);
static void bc_ox1(GridS *pGrid);

static void check_div_b(GridS *pGrid);

/* custom hst quantities */
static Real hst_m13(const GridS *pG, const int i, const int j, const int k);
static Real hst_m110(const GridS *pG, const int i, const int j, const int k);
static Real hst_mT2(const GridS *pG, const int i, const int j, const int k);
static Real hst_Mx13(const GridS *pG, const int i, const int j, const int k);
static Real hst_Erad(const GridS *pG, const int i, const int j, const int k);

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
#endif /* NSCALARS */

#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix);
static OutputS nan_dump;
static int nan_dump_count;
#endif /* REPORT_NANS */

/* ----------- Cooling functions --------- */

#ifdef ENERGY_COOLING
/* global definitions for the SD cooling curve using the
   Townsend (2009) exact integration scheme */
/* Choose cooling curve here */
#include "prob/cooling_data/SD93_Z1.h"
// #include "prob/cooling_data/denstest.h"
// #include "prob/cooling_data/WSS09_z0_Z1.h"
// #include "prob/cooling_data/powerlaw/alpha1.5.h"

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
#endif /* ENERGY_COOLING */

#if ENERGY_HEATING == 2 // fixed heating
static Real heating_rate;
#endif /* ENERGY_HEATING == 2 */

#ifdef INSTANTCOOL
static Real instant_cool(const Real rho, const Real P, const Real dt);
static int after_cool(MeshS *pM, DomainS *pDomain, int fix);
static Real hst_xshift(const GridS *pG, const int i, const int j, const int k);
#endif /* INSTANTCOOL */

/* ----------------------------------------- */

#ifdef FOLLOW_CLOUD
static Real cloud_mass_weighted_velocity(MeshS *pM);
static void boost_frame(DomainS *pDomain, Real dv);
static Real x_shift;

static Real hst_xshift(const GridS *pG, const int i, const int j, const int k);
static Real hst_v_wind(const GridS *pG, const int i, const int j, const int k);
#endif

static Real scalefac = 1.0;

/* ----------- Problem specific variables --------- */

static Real drat, dr, r_cloud, acc, v_cloud;

// Newly introduced

static Real v_wind, v_wind0, M_wind, T_cloud, rho_hot, Press, T_ceil_cool, scaling_fac;

//-------------------------------------------------------------

static Real T_floor, T_ceil, rhofloor, betafloor, T_floor_cooling; /* Used in nancheck*/
static Real dens_conv;

/* ----------------------------------------- */

static Real dtmin;
static Real pro(Real r, Real rcloud)
{
  return (r / rcloud - log(cosh(r / rcloud))) / log(2);
}

static Real get_pressure(ConsS *u)
{
  Real E0 = 0.5 * (SQR(u->M1) + SQR(u->M2) + SQR(u->M3)) / u->d;

  return (u->E - E0) * Gamma_1;
}

/*==============================================================================
 * INITIAL CONDITION:
 *
 *----------------------------------------------------------------------------*/

void problem(DomainS *pDomain)
{

  GridS *pGrid = pDomain->Grid;
  int i = 0, j = 0, k = 0;
  int is, ie, js, je, ks, ke;
  int il, iu, jl, ju, kl, ku;
  Real x1, x2, x3, r;
  Real rho, vx, vy, vz;
  Real fact;

  int iseed;
  Real3Vect ***A;
  int nterms, nx1, nx2, nx3, ierr;
  Real theta, phi, alpha, beta, amp;
  Real scal[5];
  Real my_scal[5];

  Real x1min, x1max, x2min, x2max, x3min, x3max, tmp;

  Real jmag, bmag, JdB, JcBforce, norm, Brms, Bs, Bmax, Bmax_cloud, ncells;
  Real Bin, Bout_z, Bout_y;
  Real bscale, ascale;
  int tangledloud = 0;
#if (NSCALARS > 0)
  Real dye;
#endif

  /*==============================================================================
   * Retrieve from input file                                                   */

  drat = par_getd("problem", "drat");              // rho_cold / rho_hot
  T_cloud = par_getd("problem", "T_cloud");        // Temperature of the cloud
  rho_hot = par_getd_def("problem", "rho_hot", 1); // set rho_hot default value

  r_cloud = par_getd_def("problem", "r_cloud", 0.25); // size of cloud
  dr = par_getd_def("problem", "dr", 0.0);            // boundary layer for density profile

  M_wind = par_getd("problem", "M_wind");               // Get the Mach Number of the hot wind
  v_wind = M_wind * sqrt(drat) * sqrt(Gamma * T_cloud); // Defining speed of the hot wind
  v_wind0 = v_wind;

  dtmin = par_getd_def("problem", "dtmin", 1.e-7);

  scaling_fac = par_getd("problem", "scaling_fac");
  
  Press = T_cloud * drat;

  //Real centerpos[3][3] = {{0,0,0}, {-6,-6,0}, {-6,6,0}};

  char *filename = par_gets_def("problem", "centerpos","");
  ath_pout(0, "%s", filename);

  // v_wind = par_getd("problem", "v_wind"); // TODO: change here for FLOW_PROFILE
  /* Uncomment this for constantly outflowing (fully entrained / comoving) test */
  // v_wind = 0.1 * v_wind;
  /*==============================================================================
   * Could be retrieved, but are predefined                                     */

#ifdef FOLLOW_CLOUD
  x_shift = 0.0;
#endif

  acc = par_getd_def("problem", "acceleration", 0.0); // accelerated wind
  v_cloud = par_getd_def("problem", "v_cloud", 0.0);  // extra cloud velocity

  T_floor = par_getd_def("problem", "T_floor", 0.01);     // cooling floor.
                                                          // note that there's another T_floor in the cooling function
  T_ceil = par_getd_def("problem", "T_ceil", 10. * drat); // no T above this
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-2);

  T_ceil_cool = par_getd_def("problem", "T_ceil_cool", 0.6 * T_cloud * drat); // no cooling above this

  // temp floor in cooling routine. effectively max(T_floor,T_floor_cooling) is
  // the temperture floor
  T_floor_cooling = par_getd_def("problem", "T_floor_cooling", T_cloud);

  dens_conv = par_getd_def("problem", "dens_conv", 1.0); // for density
                                                         // dependent cooling
  /*============================================================================ */

#if ENERGY_HEATING == 2                          // fixed heating
  heating_rate = par_getd("problem", "heating"); // fixed heating rate
#endif                                           /* ENERGY_HEATING == 2 */

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

#ifdef ENERGY_COOLING
  init_cooling();
  /* test_cooling(); */
#endif

#ifdef PRINTCOOLING
  Real T_start, T_end, P_const, t_crush, temp_loop;

  T_start = 1;
  T_end = 20000;

  P_const = T_cloud * drat;

  t_crush = sqrt(drat) * r_cloud / v_wind;

  ath_pout(0, "Data of cooling function with T_cloud=%f, T_hot=%f, t_cc=%f\n", T_cloud, T_cloud * drat, t_crush);
  ath_pout(0, "Going from 1 to 20.000 in steps of 0.5 Temperature Code Units\n");

  for (k = T_start; k <= T_end; k++)
  {
    temp_loop = k;
    temp_loop *= 0.5;
    ath_pout(0, "%f\t%f\n", temp_loop, sdLambda(drat, temp_loop));
  }

  ath_pout(0, "Data of cooling time with T_cloud=%f, T_hot=%f, t_cc=%f, P=%f, rho_cloud=%f", T_cloud, T_cloud * drat, t_crush, T_cloud * drat, drat);
  T_start = 1;
  T_end = 20000;

  for (k = T_start; k <= T_end; k++)
  { 
    temp_loop = k;
    temp_loop *= 0.5;
    ath_pout(0, "%f\t%f\n", temp_loop, temp_loop / (drat * sdLambda(drat, temp_loop)));
  }

  ath_pout(0, "Data of cooling time Townsend with T_cloud=%f, T_hot=%f, t_cc=%f, P=%f, rho_cloud=%f", T_cloud, T_cloud * drat, t_crush, T_cloud * drat, drat);
  T_start = 1;
  T_end = 20000;

  for (k = T_start; k <= T_end; k++)
  {
    temp_loop = k;
    temp_loop *= 0.5;
    ath_pout(0, "%f\t%f\n", temp_loop, tcool(drat, temp_loop));
  }
#endif

#ifdef INPUTFILE

#define COLUMNS 3

  int getSize(char *filename)
  {
    FILE *file;

    file = fopen(filename, "r");
    if (file == NULL)
    {
      ath_pout(0, "Error opening the file.\n Setting centerpos to default 0,0,0\n");
      Real centerpos[1][3] = {{0,0,0}};
    }

    int linecount = 1;
    char c;

    do
    {
      c = fgetc(file);
      if (c == '\n')
        linecount++;
    } while (c != EOF);

    fclose(file);

    return linecount;
  }
  
  FILE *file;

  int rows;

  rows = getSize(filename);
  float centerpos[rows][COLUMNS];

  file = fopen(filename, "r");

  while (!feof(file))
  {
    if (ferror(file))
    {
      ath_pout(0, "Error opening the file.\n Setting centerpos to default 0,0,0\n");
      Real centerpos[1][3] = {{0,0,0}};
    }
    for (int j = 0; j < rows; j++)
    {
      for (int i = 0; i < 3; i++)
      {
        if (fscanf(file, "%f", &centerpos[j][i]) == EOF)
          break;
      }
    }
  }
  fclose(file);

  ath_pout(0, "Contents of the array:\n");

  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < COLUMNS; j++)
      ath_pout(0, "%f ", centerpos[i][j]);
    printf("\n");
  }

#else

  Real centerpos[1][3] = {{0,0,0}};

#endif



#ifdef INSTANTCOOL
  //  CoolingFunc = instant_cool;
#endif

  /* Add additional output in hst file */
  dump_history_enroll(hst_m13, "m13");   // mass(rho > rho_cl / 3)
  dump_history_enroll(hst_m110, "m110"); // mass(rho > rho_cl / 10)
  dump_history_enroll(hst_mT2, "mT2");   // mass(T < T_cl /2)
  dump_history_enroll(hst_Mx13, "Mx13"); // Momentum(rho > rho_cl / 3)

  dump_history_enroll(hst_Erad, "Erad"); // Total radiated energy

#ifdef FOLLOW_CLOUD
  dump_history_enroll_alt(hst_xshift, "x_shift"); // Total shift
  dump_history_enroll_alt(hst_v_wind, "v_wind");  // Current inflow velocity
#endif

#if (NSCALARS > 0)
  // Various concentrations
  dump_history_enroll(hst_c, "<c>");
  dump_history_enroll(hst_c_sq, "<c^2>");

  dump_history_enroll(hst_cE, "<c * E>");
  dump_history_enroll(hst_cx1, "<c * x1>");

  dump_history_enroll(hst_cvx, "<c * Vx>");
  dump_history_enroll(hst_cvy, "<c * Vy>");
  dump_history_enroll(hst_cvz, "<c * Vz>");

  dump_history_enroll(hst_cvx_sq, "<(c * Vx)^2>");
  dump_history_enroll(hst_cvy_sq, "<(c * Vy)^2>");
  dump_history_enroll(hst_cvz_sq, "<(c * Vz)^2>");
  dump_history_enroll(hst_Sdye, "dye entropy");
  /*
#ifdef ENERGY_COOLING
  dump_history_enroll(hst_cstcool, "cs*tcool");
#endif
  */
#endif /* NSCALARS */

#ifdef REPORT_NANS
  nan_dump_count = 0;
#endif

  /* Initialize grid loading */
  int Nx, Ny, Nz, ndim_file[3];
  int ii, jj, kk;
  Nx = pDomain->Nx[0];
  Ny = pDomain->Nx[1];
  Nz = pDomain->Nx[2];
  char *fn_rho = par_gets_def("problem", "load_grid_rho", "");
  Real *cloud_dat_rho;
  FILE *fp_rho = NULL;
  int nreadwrite = 0;
  if ((char)fn_rho[0])
  {
    ndim_file[0] = par_geti("problem", "load_grid_ndim_x1");
    ndim_file[1] = par_geti("problem", "load_grid_ndim_x2");
    ndim_file[2] = par_geti("problem", "load_grid_ndim_x3");
    ath_pout(0, "Loading grid `%s` with (%d,%d,%d) cells (grid cells: [%d,%d,%d])\n",
             fn_rho, ndim_file[0], ndim_file[1], ndim_file[2], Nx, Ny, Nz);
    if ((ndim_file[0] > Nx) || (ndim_file[1] > Ny) || (ndim_file[2] > Nz))
      ath_error("[read_grid]: Grid in file larger than ATHENA grid.");
    fp_rho = fopen(fn_rho, "r");
    if (fp_rho == NULL)
      ath_error("[read_grid] Problem loading `%s`.", fn_rho);
    cloud_dat_rho = (Real *)malloc(sizeof(Real) * ndim_file[2] * ndim_file[1] * ndim_file[0]);
    if (cloud_dat_rho == NULL)
      ath_error("[read_grid] Error allocating memory for rho.");
    fread(cloud_dat_rho, sizeof(Real), ndim_file[0] * ndim_file[1] * ndim_file[2], fp_rho);
    fread(&rho, sizeof(Real), 1, fp_rho);
    if (!feof(fp_rho))
      ath_error("[read_grid] Not eof for rho.");
  } /* end grid loading */

  is = pGrid->is;
  ie = pGrid->ie;
  js = pGrid->js;
  je = pGrid->je;
  ks = pGrid->ks;
  ke = pGrid->ke;
  nx1 = (ie - is) + 1 + 2 * nghost;
  nx2 = (je - js) + 1 + 2 * nghost;
  nx3 = (ke - ks) + 1 + 2 * nghost;
  int iprint = 0;
  vx = vy = vz = 0.0;

  int m,arrsize;
  float d; 
  arrsize = sizeof(centerpos)/sizeof(centerpos[0]);

  bool inCloud = false;

  /* Begin cell loop */
  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = is; i <= ie; i++)
      {
        // Convert (i,j,k) to real space (x1,x2,x3)
        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* old Radius Implementation

        // Compute radius in cloud
        r = sqrt(x1 * x1 + x2 * x2 + x3 * x3); // spherical

*/
    inCloud = false;

    for (m=0; m < 3; m++)
    {
      d = sqrt((x1 - centerpos[m][0])*(x1 - centerpos[m][0])+(x2 - centerpos[m][1])*(x2 - centerpos[m][1])+(x3 - centerpos[m][2])*(x3 - centerpos[m][2]));

      if(d < r_cloud){
        inCloud = true;
        break;
      } 
    }
  
    // Static inflow
    rho = rho_hot;
    vx = v_wind;

#if (NSCALARS > 0)
        dye = 0.0;
#endif
        if (inCloud == true)
        {
          vx = v_cloud;
          vy = 0.0;
          rho *= drat;

#if (NSCALARS > 0)
          dye = drat; // 1.0;
#endif
        }
        if (dr > 0.0)
        {
          rho = (1.0 + drat * 0.5 * (1.0 + tanh((r_cloud - d) / (dr * r_cloud))));
          vx = v_wind * 0.5 * (1.0 + tanh((-(1.0 + dr) * r_cloud + d) / (dr * r_cloud))) / rho;
          vx += v_cloud;
          vx += 0.0;
          vy = 0.0;
#if (NSCALARS > 0)
          // This line means that the dye does not follow the density in the boundary (dr) region
          if (inCloud == true)
          {
            dye = rho;
          }
#endif
        }

        /* write values to the grid */
        pGrid->U[k][j][i].d = rho; // write rho, in case of hot gas =1 otherwise not
        pGrid->U[k][j][i].M1 = rho * vx;
        pGrid->U[k][j][i].M2 = rho * vy;
        pGrid->U[k][j][i].M3 = rho * vz;

#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E = Press / Gamma_1;                            // Thermal energy
        pGrid->U[k][j][i].E += 0.5 * rho * (SQR(vx) + SQR(vy) + SQR(vz)); // kinetic energy
#endif                                                                    /* not ISOTHERMAL */

#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = dye;
#endif

#ifdef ENERGY_COOLING
        pGrid->U[k][j][i].Erad = 0;
#endif
      }
    }
  } /* end grid loops */

  // Some info printed

  Real t_cc, t_cool_cl, t_cool_hot, t_cool_mix, cool_mix_cc;

  t_cc = sqrt(drat) * r_cloud / v_wind;                           // Cloud crushing time
  t_cool_cl = tcool(drat, T_cloud);                               // Cloud cooling time
  t_cool_hot = tcool(rho_hot, drat * T_cloud);                    // Cooling time of hot gas
  t_cool_mix = tcool(rho_hot * sqrt(drat), sqrt(drat) * T_cloud); // Cooling time of mixed gas
  cool_mix_cc = t_cool_mix / t_cc;                                // Ration of cool,mix over cc

  ath_pout(0, "[init_problem] t_cc = %g, t_cool,cl = %g, T_cloud = %g, t_cool,mix = %g, t_cool,hot = %g, t_cool,mix/t_cc = %g\n",
           t_cc,
#ifdef ENERGY_COOLING
           t_cool_cl,
#else
           -1, -1
#endif
           T_cloud,
           t_cool_mix,
           t_cool_hot,
           cool_mix_cc);

  // close file if ICs are read from file
  if (fp_rho != NULL)
  {
    fclose(fp_rho);
    free(cloud_dat_rho);
    if (nreadwrite != (ndim_file[0] * ndim_file[1] * ndim_file[2]))
      ath_error("Not all cells correctly read in (%d vs %d).", nreadwrite,
                ndim_file[0] * ndim_file[1] * ndim_file[2]);
  }

  if (pDomain->Disp[0] == 0)
    bvals_mhd_fun(pDomain, left_x1, bc_ix1);
  if (pDomain->MaxX[0] == pDomain->RootMaxX[0])
    bvals_mhd_fun(pDomain, right_x1, bc_ox1);

  /* seed a perturbation */
  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = is; i <= ie; i++)
      {
        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
        fact = -1.0;
        while (fabs(fact) > 0.03)
          fact = (RandomNormal(0.0, 0.01));
        if (fabs(fact) < .03)
          pGrid->U[k][j][i].d *= (1.0 + fact);
      }
    }
  }

#if ENERGY_HEATING == 2
  ath_pout(0, "Heating mode is enabled with rate %.5e (Lambda(rho_cl,T_cl) = %e).\n",
           heating_rate, sdLambda(drat, T_cloud));
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
  fwrite(&v_wind, sizeof(Real), 1, fp);
#endif
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl, nd;

  for (nl = 0; nl < (pM->NLevels); nl++)
  {
    for (nd = 0; nd < (pM->DomainsPerLevel[nl]); nd++)
    {
      if (pM->Domain[nl][nd].Disp[0] == 0)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x1, bc_ix1);
      if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0])
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, bc_ox1);
    }
  }

  /*==============================================================================
   * Retrieve from input file                                                   */

  drat = par_getd("problem", "drat");              // rho_cold / rho_hot
  T_cloud = par_getd("problem", "T_cloud");        // Temperature of the cloud
  rho_hot = par_getd_def("problem", "rho_hot", 1); // set rho_hot default value

  r_cloud = par_getd_def("problem", "r_cloud", 0.25); // size of cloud
  dr = par_getd_def("problem", "dr", 0.0);            // boundary layer for density profile

  M_wind = par_getd("problem", "M_wind");               // Get the Mach Number of the hot wind
  v_wind = M_wind * sqrt(drat) * sqrt(Gamma * T_cloud); // Defining speed of the hot wind
  v_wind0 = v_wind;

  dtmin = par_getd_def("problem", "dtmin", 1.e-7);
  scaling_fac = par_getd("problem", "scaling_fac");

  Press = T_cloud * drat;

  // v_wind = par_getd("problem", "v_wind"); // TODO: change here for FLOW_PROFILE
  /* Uncomment this for constantly outflowing (fully entrained / comoving) test */
  // v_wind = 0.1 * v_wind;
  /*==============================================================================
   * Could be retrieved, but are predefined                                     */

#ifdef FOLLOW_CLOUD
  x_shift = 0.0;
#endif

  acc = par_getd_def("problem", "acceleration", 0.0); // accelerated wind
  v_cloud = par_getd_def("problem", "v_cloud", 0.0);  // extra cloud velocity

  T_floor = par_getd_def("problem", "T_floor", 0.01);     // cooling floor.
                                                          // note that there's another T_floor in the cooling function
  T_ceil = par_getd_def("problem", "T_ceil", 10. * drat); // no T above this
  rhofloor = par_getd_def("problem", "rhofloor", 1.e-2);

  T_ceil_cool = par_getd_def("problem", "T_ceil_cool", 10); // no cooling above this

  // temp floor in cooling routine. effectively max(T_floor,T_floor_cooling) is
  // the temperture floor
  T_floor_cooling = par_getd_def("problem", "T_floor_cooling", T_cloud);

  dens_conv = par_getd_def("problem", "dens_conv", 1.0); // for density
                                                         // dependent cooling
  /*============================================================================ */

#ifdef ENERGY_COOLING
  init_cooling();
#endif

#ifdef INSTANTCOOL
  //  CoolingFunc = instant_cool;
#endif

  // Re-enroll hst dumps after restart
  dump_history_enroll(hst_m13, "m13");
  dump_history_enroll(hst_m110, "m110");
  dump_history_enroll(hst_mT2, "mT2");
  dump_history_enroll(hst_Mx13, "Mx13");

  dump_history_enroll(hst_Erad, "Erad");

#ifdef FOLLOW_CLOUD
  dump_history_enroll_alt(hst_xshift, "x_shift");
  dump_history_enroll_alt(hst_v_wind, "v_wind");
#endif

#if (NSCALARS > 0)
  dump_history_enroll(hst_c, "<c>");
  dump_history_enroll(hst_c_sq, "<c^2>");

  dump_history_enroll(hst_cE, "<c * E>");
  dump_history_enroll(hst_cx1, "<c * x1>");

  dump_history_enroll(hst_cvx, "<c * Vx>");
  dump_history_enroll(hst_cvy, "<c * Vy>");
  dump_history_enroll(hst_cvz, "<c * Vz>");

  dump_history_enroll(hst_cvx_sq, "<(c * Vx)^2>");
  dump_history_enroll(hst_cvy_sq, "<(c * Vy)^2>");
  dump_history_enroll(hst_cvz_sq, "<(c * Vz)^2>");
  dump_history_enroll(hst_Sdye, "dye entropy");
#endif /* NSCALARS */

  /*
#ifdef ENERGY_COOLING
  dump_history_enroll(hst_cstcool, "cs*tcool");
#endif
  */

  /* DANGER: make sure the order here matches the order in write_restart() */
#ifdef FOLLOW_CLOUD
  fread(&x_shift, sizeof(Real), 1, fp);
  fread(&v_wind, sizeof(Real), 1, fp);
#endif
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}

void Userwork_before_loop(MeshS *pM)
{
  int nl, nd, ntot;

  /* report nans first, so we can fix them before they propagate into
     the following functions. */
  for (nl = 0; nl <= (pM->NLevels) - 1; nl++)
  {
    for (nd = 0; nd <= (pM->DomainsPerLevel[nl]) - 1; nd++)
    {
      if (pM->Domain[nl][nd].Grid != NULL)
      {
#ifdef REPORT_NANS
        ntot = report_nans(pM, &(pM->Domain[nl][nd]), 1);
        // if(ntot > 0)
        // report_nans(pM, &(pM->Domain[nl][nd]),1);
#endif
#ifdef INSTANTCOOL
        after_cool(pM, &(pM->Domain[nl][nd]), 1);
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
  // Artificial shift
  // dvx = 4e-4 * MAX(0, (1 - (r0 + x_shift * x_shift) / (4 * r0))) + 1e-3 / (1 + x_shift) + MIN(1e-9 * x_shift * x_shift, 1e-3);

  /* if(pM->time < 2)
     dvx = 0;
  */

  if ((dvx < 0.0) || isnan(dvx))
  {
    ath_pout(0, "[bad dvx:] %0.15e setting to 0.\n", dvx);
    dvx = 0.0;
  }

  /* Enforcing ceiling to  shift */
  if (dvx > 1e-2)
  {
    dvx = 1e-2;
    ath_pout(0, "[bad dvx:] %0.15e setting to 0.01\n", dvx);
  }

  if (dvx > 0.0)
  {
    expt = floor(log10(dvx));
    newdvx = dvx / pow(10, expt);
    newdvx = floor(newdvx * 1.0e4) / 1.0e4;
    newdvx = newdvx * pow(10.0, expt);
    dvx = newdvx;
  }

  if (v_wind - dvx < 0.00)
  { // does not allow v_wind < 0
    dvx = v_wind;
  }
  ath_pout(0, "[dvx:]  %0.10e [v_wind:] %.10e [xshift:] %.10e\n", dvx, v_wind, x_shift);
  v_wind -= dvx;
#endif /* FOLLOW_CLOUD */

  for (nl = 0; nl <= (pM->NLevels) - 1; nl++)
  {
    for (nd = 0; nd <= (pM->DomainsPerLevel[nl]) - 1; nd++)
    {
      if (pM->Domain[nl][nd].Grid != NULL)
      {
#ifdef FOLLOW_CLOUD
        boost_frame(&(pM->Domain[nl][nd]), dvx);
#endif
      }
    }
  }

  for (nl = 0; nl <= (pM->NLevels) - 1; nl++)
  {
    for (nd = 0; nd <= (pM->DomainsPerLevel[nl]) - 1; nd++)
    {
      if (pM->Domain[nl][nd].Grid != NULL)
      {
      }
    }
  }

  for (nl = 0; nl <= (pM->NLevels) - 1; nl++)
  {
    for (nd = 0; nd <= (pM->DomainsPerLevel[nl]) - 1; nd++)
    {
      if (pM->Domain[nl][nd].Grid != NULL)
      {
#ifdef ENERGY_COOLING
        integrate_cooling(pM->Domain[nl][nd].Grid);
#endif
      }
    }
  }
#if ENERGY_HEATING == 1
  radiate_energy(pM);
#endif

  // Accelerate
  v_wind += acc * pM->dt;

  if (pM->dt < dtmin)
  {
    data_output(pM, 1);
    ath_error("dt too small\n");
  }

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
 * ed_exp_mode()         - returns current expansion mode
 * expand_domain()       - expands domain by scalefactor
 * report_nans()         - apply a ceiling and floor to the temperature
 * cloud_velocity()      - find the mass-weighted velocity of the cloud.
 * nu_fun()              - kinematic viscosity (i.e., cm^2/s)
 * cooling_func()        - cooling function for the energy equation
 *----------------------------------------------------------------------------*/
#ifdef FOLLOW_CLOUD
static void boost_frame(DomainS *pDomain, Real dvx)
{
  int i, j, k;
  int is, ie, js, je, ks, ke;
  Real d;

  GridS *pGrid = pDomain->Grid;
  is = pGrid->is;
  ie = pGrid->ie;
  js = pGrid->js;
  je = pGrid->je;
  ks = pGrid->ks;
  ke = pGrid->ke;

  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = is; i <= ie; i++)
      {
        d = pGrid->U[k][j][i].d;

#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += 0.5 * d * SQR(dvx);
        pGrid->U[k][j][i].E -= dvx * pGrid->U[k][j][i].M1;
#endif /* ISOTHERMAL */
        pGrid->U[k][j][i].M1 -= dvx * d;
      }
    }
  }

  //  x_shift -= dvx * pDomain->Grid->dt;
  x_shift -= v_wind * pDomain->Grid->dt;
  return;
}
#endif /* FOLLOW_CLOUD */

#ifdef REPORT_NANS
static int report_nans(MeshS *pM, DomainS *pDomain, int fix)
{
#ifndef ISOTHERMAL
  int i, j, k;
  int is, ie, js, je, ks, ke;
  Real x1, x2, x3;
  int V = 1; // verbose off = 0
  int NO = 5;
  Real KE, rho, press, temp;
  int nanpress = 0, nanrho = 0, nanv = 0, nnan; /* nan count */
  int npress = 0, nrho = 0, nv = 0, nfloor;     /* floor count */
  Real beta;
  Real scal[8];
#ifdef MPI_PARALLEL
  Real my_scal[8];
  int ierr;
#endif

  /*Real T_floor    = 1.0e-2 / drat;
  Real T_ceil     = 100.0;
  Real rhofloor  = 1.0e-2;
  Real betafloor = 3.0e-3;
  */
  GridS *pGrid = pDomain->Grid;

  is = pGrid->is;
  ie = pGrid->ie;
  js = pGrid->js;
  je = pGrid->je;
  ks = pGrid->ks;
  ke = pGrid->ke;

  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = is; i <= ie; i++)
      {
        rho = pGrid->U[k][j][i].d;
        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
        KE = (SQR(pGrid->U[k][j][i].M1) +
              SQR(pGrid->U[k][j][i].M2) +
              SQR(pGrid->U[k][j][i].M3)) /
             (2.0 * rho);

        press = pGrid->U[k][j][i].E - KE;
        press *= Gamma_1;
        temp = press / rho;
        beta = fabs(200. * betafloor);

        if (press != press)
        {
          nanpress++;
          if (V && nanpress < NO)
            printf("bad press %e R %e  %e %e %e  %d %d %d %e %e %e %e\n", press, 1.0 / pGrid->dx1, x1, x2, x3, i, j, k, rho, press, temp, beta);
          if (fix)
            temp = T_floor;
        }
        else if (temp < T_floor)
        {
          npress++;
          if (V && npress < NO)
            printf("bad tempF %e R %e  %e %e %e  %d %d %d %e %e %e %e\n", temp, 1.0 / pGrid->dx1, x1, x2, x3, i, j, k, rho, press, temp, beta);
          if (fix)
            temp = T_floor;
        }
        else if (temp > T_ceil && beta > 10.0 * betafloor)
        {
          npress++;
          if (V && npress < NO)
            printf("bad tempC %e R %e  %e %e %e  %d %d %d %e %e %e %e\n", temp, 1.0 / pGrid->dx1, x1, x2, x3, i, j, k, rho, press, temp, beta);
          if (fix)
            temp = T_ceil;
        }

        if (rho != rho)
        {
          nanrho++;
          if (V && nanrho < NO)
            printf("bad rho %e R %e  %e %e %e  %d %d %d\n", rho, 1.0 / pGrid->dx1, x1, x2, x3, i, j, k);
          if (fix)
            rho = rhofloor;
        }
        else if (rho < rhofloor)
        {
          nrho++;
          if (V && nrho < NO)
            printf("bad rho %e R %e  %e %e %e  %d %d %d\n", rho, 1.0 / pGrid->dx1, x1, x2, x3, i, j, k);
          if (fix)
            rho = rhofloor;
        }

        if (pGrid->U[k][j][i].M1 != pGrid->U[k][j][i].M1)
        {
          nanv++;
          if (fix)
            pGrid->U[k][j][i].M1 = 0.0;
        }
        if (pGrid->U[k][j][i].M2 != pGrid->U[k][j][i].M2)
        {
          nanv++;
          if (fix)
            pGrid->U[k][j][i].M2 = 0.0;
        }
        if (pGrid->U[k][j][i].M3 != pGrid->U[k][j][i].M3)
        {
          nanv++;
          if (fix)
            pGrid->U[k][j][i].M3 = 0.0;
        }
        /* write values back to the grid */
        /* TODO: what about B??? */
        if (fix)
        {
          pGrid->U[k][j][i].d = rho;
          KE = (SQR(pGrid->U[k][j][i].M1) +
                SQR(pGrid->U[k][j][i].M2) +
                SQR(pGrid->U[k][j][i].M3)) /
               (2.0 * rho);

          pGrid->U[k][j][i].E = temp * rho / Gamma_1 + KE;
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
  nanrho = scal[1];
  nanv = scal[2];
  npress = scal[4];
  nrho = scal[5];
  nv = scal[6];
#endif /* MPI_PARALLEL */

  /* sum up the # of bad cells and report */
  nnan = nanpress + nanrho + nanv;
  /* sum up the # of floored cells and report */
  nfloor = npress + nrho + nv;
  if (nfloor > 0)
  {
    ath_pout(0, "[report_nans]: floored %d cells: %d P, %d d, %d v.\n",
             nfloor, npress, nrho, nv);
  }

  //  if ((nnan > 0 || nfloor -nmag > 30) && fix == 0) {
  if (nnan > 0)
  { // && fix == 0) {
    ath_pout(0, "[report_nans]: found %d nan cells: %d P, %d d, %d v.\n",
             nnan, nanpress, nanrho, nanv);

    nan_dump.n = 100;
    nan_dump.dt = HUGE_NUMBER;
    nan_dump.t = pM->time;
    nan_dump.num = 1000 + nan_dump_count;
    nan_dump.out = "prim";
    nan_dump.nlevel = -1; /* dump all levels */

    dump_vtk(pM, &nan_dump);
    if (nnan)
      nan_dump_count++;
    if (nan_dump_count > 10)
      ath_error("[report_nans]: too many nan'd timesteps.\n");

    if (nfloor > 1000)
      ath_error("[report_nans]: Too many floored cells.\n");
  }

#endif /* ISOTHERMAL */

  return nfloor + nnan;
}
#endif /* REPORT_NANS */

#ifdef ENERGY_COOLING
/* ================================================================ */
/* cooling routines */

static void init_cooling()
{
  int i, k, n = nfit_cool_T - 1;
  Real term;
  const Real mu = 0.62, mu_e = 1.17;

  const Real conv_fac = 8.61e-4; // conversion from 1e4 K to KeV !!!!!!!

  /* convert T in the cooling function from keV to code units */
  for (k = 0; k <= n; k++)
  {
    sdT[k] /= conv_fac;
    if (sdT[k] <= 0)
      ath_error("sdT[%d]=%e. Has to be > 0.", k, sdT[k]);
  }

  if (T_floor_cooling < sdT[0])
    ath_error("Cooling floor is smaller than first entry of cooling function (%e vs %e).",
              T_floor_cooling, sdT[0]);

  /* populate Yk following equation A6 in Townsend (2009) */
  for (i = 0; i < nfit_cool_d; i++)
  {
    Yk[i][n] = 0.0;
    for (k = n - 1; k >= 0; k--)
    {
      if (sdL[i][k] <= 0)
        ath_error("sdL[%d][%d]=%e. Has to be > 0.", i, k, sdL[i][k]);
      term = (sdL[i][n] / sdL[i][k]) * (sdT[k] / sdT[n]);

      if (sdexpt[i][k] == 1.0)
        term *= log(sdT[k] / sdT[k + 1]);
      else
        term *= ((1.0 - pow(sdT[k] / sdT[k + 1], sdexpt[i][k] - 1.0)) / (1.0 - sdexpt[i][k]));

      Yk[i][k] = Yk[i][k + 1] - term;

      if (isnan(Yk[i][k]))
        ath_error("Error initializing cooling. nan in Yk[%d][%d]", i, k);
    }
  }
  return;
}

/* piecewise power-law fit to the cooling curve with temperature in
   keV and L in 1e-23 erg cm^3 / s
   T can be handed over in code units, conversion has been done in init_cooling

   */
static Real sdLambda(const Real d0, const Real T)
{
  int iT, id; // bin indices for T,d
  Real L1, L2;
  const Real conv_fac = scaling_fac; // from units of 1e-23 erg cm^3 /s to code units.

  // const Real conv_fac = 1.311e-5; // from units of 1e-23 erg cm^3 /s to code units.
  const Real d = d0 * dens_conv;
  int interpolate = (nfit_cool_d > 1);

  /* first find the temperature bin */
  for (iT = nfit_cool_T - 1; iT >= 0; iT--)
  {
    if (T >= sdT[iT])
      break;
  }
  if (iT < 0)
    ath_error("[sdLambda] T %e d %e %e %d\n", T, d, sdT[0], iT);

  /* Find the density bin */
  if (d <= sdd[0])
  {
    interpolate = 0;
    id = 0;
  }
  else if (d >= sdd[nfit_cool_d - 1])
  {
    interpolate = 0;
    id = nfit_cool_d - 1;
  }
  else
  {
    for (id = nfit_cool_d - 1; id >= 0; id--)
    {
      if (d >= sdd[id])
        break;
    }
  }

  /* piecewise power-law; see equation A4 of Townsend (2009) */
  L1 = conv_fac * sdL[id][iT] * pow(T / sdT[iT], sdexpt[id][iT]);

  if (!interpolate)
    return L1;

  L2 = conv_fac * sdL[id + 1][iT] * pow(T / sdT[iT], sdexpt[id + 1][iT]);

  // Linear interpolation
  return L1 + (L2 - L1) / (sdd[id + 1] - sdd[id]) * d;
}

static Real tcool(const Real d, const Real T)
{
  const Real mu = 0.62, mu_e = 1.17;

  /* equation 13 of Townsend (2009) */
  return (SQR(mu_e) * T) / (Gamma_1 * d * sdLambda(d, T));
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
  for (iT = nT; iT >= 0; iT--)
  {
    if (T >= sdT[iT])
      break;
  }

  /* calculate Y using equation A5 in Townsend (2009) */
  term = (sdL[id][nT] / sdL[id][iT]) * (sdT[iT] / sdT[nT]);

  if (sdexpt[id][iT] == 1.0)
    term *= log(sdT[iT] / T);
  else
    term *= ((1.0 - pow(sdT[iT] / T, sdexpt[id][iT] - 1.0)) / (1.0 - sdexpt[id][iT]));

  return (Yk[id][iT] + term);
}

static Real Yinv(const Real Y1, const int id)
{
  // int iT,id;
  int nT = nfit_cool_T - 1;
  int nd = nfit_cool_d - 1;
  int iT;
  Real term;

  /* find the bin i in which the final temperature will be */
  for (iT = nT; iT > 0; iT--)
  {                           // use iT>0 instead of iT>=0 to force min(iT)=0
    if (Y(sdT[iT], id) >= Y1) // this means that newT<sdT[0] are wrong
      break;                  // but since we have T_floor>sdT[0] we're good.
  }

  /* calculate Yinv using equation A7 in Townsend (2009) */
  term = (sdL[id][iT] / sdL[id][nT]) * (sdT[nT] / sdT[iT]);
  term *= (Y1 - Yk[id][iT]);

  if (sdexpt[id][iT] == 1.0)
    term = exp(-1.0 * term);
  else
  {
    term = pow(1.0 - (1.0 - sdexpt[id][iT]) * term,
               1.0 / (1.0 - sdexpt[id][iT]));
  }

  if (!isfinite(term))
  {
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

  if (T <= T_floor_cooling)
    return T_floor_cooling;

  Tref = sdT[nfit_cool_T - 1];
  dref = sdd[nfit_cool_d - 1] / dens_conv;

  /* Find the density bin */
  if (d <= sdd[0])
  {
    interpolate = 0;
    id = 0;
  }
  else if (d >= sdd[nfit_cool_d - 1])
  {
    interpolate = 0;
    id = nfit_cool_d - 1;
  }
  else
  {
    for (id = nfit_cool_d - 1; id >= 0; id--)
    {
      if (d >= sdd[id])
        break;
    }
  }
  assert(id >= 0);
  assert(id < nfit_cool_d);

  /* Calculate new T */
  term1 = (T / Tref) * (sdLambda(d0, Tref) / sdLambda(d0, T)) * (dt_hydro / tcool(d0, T));
  T1 = Yinv(Y(T, id) + term1, id);
  if (isnan(T1))
  {
    printf("[newtemp] d=%.3e, id=%d, T=%.3e --> %.3e, dt_hydro=%e, tcool=%e, Tref=%e, sdLambda=%e, term1=%e, Y(T,id)=%e\n", d, id, T, T1, dt_hydro, tcool(d0, T), Tref, sdLambda(d0, T), term1, Y(T, id));
    assert(!isnan(T1));
  }
  if (!interpolate)
    return T1;

  T2 = Yinv(Y(T, id + 1) + term1, id + 1);
  assert(!isnan(T2));

  /* Linear interpolation along density */
  return T1 + (T2 - T1) / (sdd[id + 1] - sdd[id]) * d;
}

static void integrate_cooling(GridS *pG)
{
  int i, j, k, iprint = 0;
  int is, ie, js, je, ks, ke;

  PrimS W;
  ConsS U;
  // Changed this for T_floor!!
  Real temp, tempold, heat;

  /* ath_pout(0, "integrating cooling using Townsend (2009) algorithm.\n"); */

  is = pG->is;
  ie = pG->ie;
  js = pG->js;
  je = pG->je;
  ks = pG->ks;
  ke = pG->ke;

  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = is; i <= ie; i++)
      {

        W = Cons_to_Prim(&(pG->U[k][j][i]));
        pG->U[k][j][i].Erad = 0;

        temp = W.P / W.d;
        tempold = temp;

        /* do not cool above a certain threshold */
        if ((T_ceil_cool > 0) && (temp > T_ceil_cool))
          continue;

        temp = newtemp_townsend(W.d, temp, pG->dt);

        /* apply a temperature floor (nans tolerated) */
        if (isnan(temp) || temp < T_floor_cooling)
          temp = T_floor_cooling;

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
static void radiate_energy(MeshS *pM)
{
  ath_error("Does not work right now.");
  Real Erad_total = 0;
  GridS *pG;
  int i, j, k, is, ie, js, je, ks, ke;
  int nl, nd;
  Real dV;
  Real V = (pM->RootMaxX[2] - pM->RootMinX[2]) *
           (pM->RootMaxX[1] - pM->RootMinX[1]) *
           (pM->RootMaxX[0] - pM->RootMinX[0]);

#ifdef MPI_PARALLEL
  int ierr;
  Real my_Etot;
#endif

  // Calculate total energy
  for (nl = 0; nl <= (pM->NLevels) - 1; nl++)
  {
    for (nd = 0; nd <= (pM->DomainsPerLevel[nl]) - 1; nd++)
    {
      if (pM->Domain[nl][nd].Grid != NULL)
      {
        pG = pM->Domain[nl][nd].Grid;
        ath_poutfor
            is = pG->is;
        ie = pG->ie;
        js = pG->js;
        je = pG->je;
        ks = pG->ks;
        ke = pG->ke;

        for (k = ks; k <= ke; k++)
        {
          for (j = js; j <= je; j++)
          {
            for (i = is; i <= ie; i++)
            {
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
  for (nl = 0; nl <= (pM->NLevels) - 1; nl++)
  {
    for (nd = 0; nd <= (pM->DomainsPerLevel[nl]) - 1; nd++)
    {
      if (pM->Domain[nl][nd].Grid != NULL)
      {
        pG = pM->Domain[nl][nd].Grid;

        is = pG->is;
        ie = pG->ie;
        js = pG->js;
        je = pG->je;
        ks = pG->ks;
        ke = pG->ke;

        dV = pG->dx1 * pG->dx2 * pG->dx3;
        for (k = ks; k <= ke; k++)
        {
          for (j = js; j <= je; j++)
          {
            for (i = is; i <= ie; i++)
            {
              pG->U[k][j][i].E += Erad_total * dV / V;
            }
          }
        }
      }
    }
  }
  ath_pout(0, "[radiate_energy] Distributed a total energy of %e in fractional "
              "volumes of %e.\n",
           Erad_total, dV / V);
}
#endif /* ENERGY_HEATING */

static void test_cooling()
{
  int i, j, npts = 100;
  Real logt, temp, tc, logdt, dt;
  Real err;
  Real dens = 1.0;

  FILE *outfile;

  /* this sometimes crashes, so let's try it */
  newtemp_townsend(sdd[nfit_cool_d - 1] - 1e-4, 7.173e-04, 8.050312e-03);

  /* Now outputting some information from Townsend (2009) */
  outfile = fopen("lambda.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logt = log(1.0e-4) + (log(5.0) - log(1.0e-4)) * ((double)i / (npts - 1));
    temp = exp(logt);

    fprintf(outfile, "%e\t%e\t%e\n", temp, sdLambda(dens, temp), sdLambda(dens * drat, temp));
  }
  fclose(outfile);

  outfile = fopen("Yk.dat", "w");
  for (i = 0; i < nfit_cool_T; i++)
  {
    fprintf(outfile, "%d\t", i);
    for (j = 0; j < nfit_cool_d; j++)
      fprintf(outfile, "%e\t", Yk[j][i]);
    fprintf(outfile, "\n");
  }
  fclose(outfile);

  temp = 10.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-10kev.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logdt = log(0.1) + (log(2.0) - log(0.1)) * ((double)i / (npts - 1));
    dt = tc * exp(logdt);
    fprintf(outfile, "%e\t%e\n", dt / tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 3.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-3kev.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logdt = log(0.1) + (log(2.0) - log(0.1)) * ((double)i / (npts - 1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt / tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 1.0;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-1kev.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logdt = log(0.1) + (log(2.0) - log(0.1)) * ((double)i / (npts - 1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt / tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 0.3;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.3kev.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logdt = log(0.1) + (log(2.0) - log(0.1)) * ((double)i / (npts - 1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt / tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 0.1;
  tc = tcool(1.0, temp);

  outfile = fopen("townsend-fig1-0.1kev.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logdt = log(0.1) + (log(2.0) - log(0.1)) * ((double)i / (npts - 1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt / tc, newtemp_townsend(1.0, temp, dt));
  }
  fclose(outfile);

  temp = 0.1;
  tc = tcool(100.0, temp);

  outfile = fopen("townsend-fig1-0.1kev-100.dat", "w");
  for (i = 0; i < npts; i++)
  {
    logdt = log(0.1) + (log(2.0) - log(0.1)) * ((double)i / (npts - 1));
    dt = tc * exp(logdt);

    fprintf(outfile, "%e\t%e\n", dt / tc, newtemp_townsend(100.0, temp, dt));
  }
  fclose(outfile);

  ath_error("check cooling stuff done.\n");

  return;
}
/* end cooling routines */
/* ================================================================ */
#endif /* ENERGY_COOLING */

#ifdef INSTANTCOOL
static Real instant_cool(const Real rho, const Real P, const Real dt)
{
  /*returns cooling to decrease temperature to initial cloud temperature */
  Real temp = P / rho;
  Real tcloud = Gamma_1 / drat;
  Real Edot = 0.0;
  /* cool to tcloud, but no further */
  if ((temp > tcloud) && (temp < 10. * tcloud))
  {
    Edot = (temp - tcloud) * rho / dt / Gamma_1;
  }
  // Edot = (temp < tcloud) ? 0.0 : (temp-tcloud)/dt;

  return Edot;
}

static int after_cool(MeshS *pM, DomainS *pDomain, int fix)
{
  int i, j, k;
  int is, ie, js, je, ks, ke;
  Real x1, x2, x3;
  int V = 0; // verbose off
  int NO = 2;
  Real KE, rho, press, temp;

  GridS *pGrid = pDomain->Grid;
  Real tcloud = T_cloud;
  is = pGrid->is;
  ie = pGrid->ie;
  js = pGrid->js;
  je = pGrid->je;
  ks = pGrid->ks;
  ke = pGrid->ke;

  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = is; i <= ie; i++)
      {
        rho = pGrid->U[k][j][i].d;
        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
        KE = (SQR(pGrid->U[k][j][i].M1) +
              SQR(pGrid->U[k][j][i].M2) +
              SQR(pGrid->U[k][j][i].M3)) /
             (2.0 * rho);

        press = pGrid->U[k][j][i].E - KE;
        press *= Gamma_1;
        temp = press / rho;

        if ((temp > 1.1 * tcloud) && (temp < 10. * tcloud) && (pGrid->U[k][j][i].s[0] >= 0.1))
        {
          temp = 1.1 * tcloud;

          pGrid->U[k][j][i].E = temp * rho / Gamma_1 + KE;
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
  for (nd = 0; nd < (pM->DomainsPerLevel[nl]); nd++)
  {
    if (pM->Domain[nl][nd].Grid != NULL)
    {

      pG = pM->Domain[nl][nd].Grid;
      is = pG->is;
      ie = pG->ie;
      js = pG->js;
      je = pG->je;
      ks = pG->ks;
      ke = pG->ke;

      for (k = ks; k <= ke; k++)
      {
        for (j = js; j <= je; j++)
        {
          for (i = is; i <= ie; i++)
          {
            d = pG->U[k][j][i].d;
#if (NSCALARS > 0)
            s = pG->U[k][j][i].s[0];
#endif
            cc_pos(pG, i, j, k, &x1, &x2, &x3);
            tmp = s * pG->U[k][j][i].M1 / d;
            if (tmp == tmp && x1 < 0.0)
            {

              scal[0] += tmp * d;
              scal[1] += s * d;
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
#endif /* FOLLOW_CLOUD */

/*==============================================================================
 * HISTORY OUTPUTS:
 *

 *----------------------------------------------------------------------------*/

static Real _hst_mcut(const GridS *pG, const int i, const int j, const int k, const Real frac)
{
  Real s = 1.0;
  if (pG->U[k][j][i].d < frac * drat * s)
    return 0;
  return pG->U[k][j][i].d;
}

static Real hst_m13(const GridS *pG, const int i, const int j, const int k)
{
  return _hst_mcut(pG, i, j, k, 1 / 3.);
}

static Real hst_mT2(const GridS *pG, const int i, const int j, const int k)
{
  Real temp = get_pressure(&(pG->U[k][j][i])) / pG->U[k][j][i].d;
  const Real Tcl = T_cloud;
  if (temp > 2 * Tcl)
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
  if (pG->U[k][j][i].d < drat / 3. * s)
    return 0;
  return pG->U[k][j][i].M1;
}

static Real hst_Erad(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].Erad;
}

#ifdef FOLLOW_CLOUD
static Real hst_xshift(const GridS *pG, const int i, const int j, const int k)
{
  return x_shift;
}
#endif

#ifdef FOLLOW_CLOUD
static Real hst_v_wind(const GridS *pG, const int i, const int j, const int k)
{
  return v_wind;
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
  cc_pos(pG, i, j, k, &x1, &x2, &x3);
  return (pG->U[k][j][i].s[0] * x1);
}

static Real hst_cvx(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M1 / pG->U[k][j][i].d);
}

static Real hst_cvy(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M2 / pG->U[k][j][i].d);
}

static Real hst_cvz(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].s[0] * pG->U[k][j][i].M3 / pG->U[k][j][i].d);
}

static Real hst_cvx_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0] * pG->U[k][j][i].M1 / pG->U[k][j][i].d));
}

static Real hst_cvy_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0] * pG->U[k][j][i].M2 / pG->U[k][j][i].d));
}

static Real hst_cvz_sq(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].s[0] * pG->U[k][j][i].M3 / pG->U[k][j][i].d));
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
  // return 1.;
}
#endif

#endif /* NSCALARS */

/*==============================================================================
 * BOUNDARY CONDITIONS:
 *
 *----------------------------------------------------------------------------*/

static void bc_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i, j, k;
  Real presswind = T_cloud * drat;

  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = 1; i <= nghost; i++)
      {
        pGrid->U[k][j][is - i] = pGrid->U[k][j][is];

#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = 0.0;
#endif

        pGrid->U[k][j][is - i].d = 1.0;
        pGrid->U[k][j][is - i].M1 = 1.0 * v_wind;
        pGrid->U[k][j][is - i].M2 = 0.0;
        pGrid->U[k][j][is - i].M3 = 0.0;
        pGrid->U[k][j][is - i].E = presswind / Gamma_1 + 0.5 * SQR(v_wind);

        if ((pGrid->U[k][j][is - i].E < 0) || isnan(pGrid->U[k][j][is - i].E))
          ath_error("[bc_ix1] E %e %e %e %e %e %d %d %d\n",
                    pGrid->U[k][j][is - i].E, pGrid->U[k][j][is - i].M1, pGrid->U[k][j][is - i].M2,
                    pGrid->U[k][j][is - i].M3, pGrid->U[k][j][is - i].d, k, j, is - i);
      }
    }
  }

  return;
}

static void bc_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i, j, k;
  const int V = 0;
  int NO = 10;

  for (k = ks; k <= ke; k++)
  {
    for (j = js; j <= je; j++)
    {
      for (i = 1; i <= nghost; i++)
      {

        pGrid->U[k][j][ie + i] = pGrid->U[k][j][ie];

        if (pGrid->U[k][j][ie + i].M1 < 0.0)
        {
          if (V && (NO > 0))
          {
            printf("bc_ox1 %d %d %d %e\n",
                   i, j, k, pGrid->U[k][j][ie + i].M1);
            NO--;
          }
          pGrid->U[k][j][ie + i].E -= 0.5 * SQR(pGrid->U[k][j][ie + i].M1) / pGrid->U[k][j][ie + i].d;
          pGrid->U[k][j][ie + i].M1 = 0.0;
        }
      }
    }
  }

  return;
}