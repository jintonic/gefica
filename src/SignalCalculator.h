#ifndef GeFiCa_SignalCalculator
#define GeFiCa_SignalCalculator
/* Reference temperature for drift vel. corrections is 77K */
#define REF_TEMP 77.0
/* max, min temperatures for allowed range */
#define MIN_TEMP 40.0
#define MAX_TEMP 120.0
#define MAX_LINE 512
#define NET_SIGNAL_THRESH 0.55
#define WP_THRESH 0.55
#define WP_THRESH_ELECTRONS 1e-4 /*electrons are considered collected if
           they stop drifting where the wp is < this*/
// #include "PointContact.h"
namespace GeFiCa { class SignalCalculator;}

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include <vector>
#include "TMath.h"
using namespace TMath;
using namespace ROOT::Math;
// using namespace GeFiCa;
// class TTree;
#include "TTree.h"
#include "RhoZ.h"
#include "PointContact.h"

// from fields.c
struct velocity_lookup{
  float e;
  float e100;
  float e110;
  float e111;
  float h100;
  float h110;
  float h111;
  float ea; //coefficients for anisotropic drift 
  float eb;
  float ec;
  float ebp;
  float ecp;
  float ha;
  float hb;
  float hc;
  float hbp;
  float hcp;
  float hcorr;
  float ecorr;
};

class GeFiCa::SignalCalculator
{
public:
  SignalCalculator(const char*, const char*);
  // ~SignalCalculator();
  int MakeSignal(RhoZPhiPoint, float*, float);
  int GetSignal(RhoZPhiPoint, float*);
  int RcIntegrate(float*, float*, float, int);
  char *pt_to_str(char *str, int len, RhoZPhiPoint pt);
  void gears(const char *gears_file);
  int  PrintSignal(RhoZPhiPoint pt);

  int SetupVelo(const char*);
  int  efield_exists(RhoZPhiPoint pt, int* idxs);
  int  grid_weights(RhoZPhiPoint pt, int* idxs, double weights[4]);
  RhoZPhiVector efield(RhoZPhiPoint pt, int* idxs);
  int wpotential(RhoZPhiPoint pt, float *wp);
  int drift_velocity(RhoZPhiPoint pt, float q, RhoZPhiVector *velo);
  std::vector<double> cs1, cs2;

  RhoZ *grd;
  PointContact *det;
  ClassDef(SignalCalculator,1);

  // signal calculation 
  float xtal_temp;            // crystal temperature in Kelvin
  float preamp_tau;           // integration time constant for preamplifier, in ns
  int   time_steps_calc;      // number of time steps used in calculations
  float step_time_calc;       // length of time step used for calculation, in ns
  float step_time_out;        // length of time step for output signal, in ns
  //    nonzero values in the next few lines significantly slow down the code
  float charge_cloud_size;    // initial FWHM of charge cloud, in mm; set to zero for point charges
  int   use_diffusion;        // set to 0/1 for ignore/add diffusion as the charges drift
  float energy;               // set to energy > 0 to use charge cloud self-repulsion, in keV

protected:

  int N;

  char drift_name[256];       // drift velocity lookup table
  int  ntsteps_out = 100;

  // data for calc_signal.c
  RhoZPhiPoint *dpath_e, *dpath_h;        // electron and hole drift paths
  float surface_drift_vel_factor;  // ratio of velocity on passivated surface rather than in bulk
  float initial_vel, final_vel;    // initial and final drift velocities for charges collected to PC
  float dv_dE;     // derivative of drift velocity with field ((mm/ns) / (V/cm))
  float v_over_E;  // ratio of drift velocity to field ((mm/ns) / (V/cm))
  double final_charge_size;     // in mm
  double diffusion_coef;

  // to be added
  float impurity_z0;

  int  v_lookup_len = 0;
  std::vector<velocity_lookup> v_lookup;
  char driftName[256];       // drift velocity lookup table

};

#endif