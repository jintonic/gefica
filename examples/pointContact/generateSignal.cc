#include "Math/Point3D.h"
#include "Math/Vector3D.h"
using namespace ROOT::Math;

using namespace GeFiCa;
// run gears with pc.mac to generate pc.root
void generateSignal(const char *gears_file="pc.root")
{
  SignalCalculator SigCal("ICPC.root", "../../src/drift_vel_tcorr.tab");

  // configuration for signal calculation
  SigCal.xtal_temp         = 90;
  SigCal.preamp_tau        = 30;
  SigCal.time_steps_calc   = 80000;
  SigCal.step_time_calc    = 1.0;
  SigCal.step_time_out     = 10.0;
  SigCal.charge_cloud_size = 0;
  SigCal.use_diffusion     = 0;

  // RhoZPhiPoint pt(-1.8*cm, 1.4*cm, 0);
  // for (int i=0; i<10; i++){
  //   for (int j=0; j<10; j++){
  //     pt.SetRho(i*0.2*cm);
  //     pt.SetZ(j*0.2*cm);
  //     SigCal.PrintSignal(pt);
  //   }
  // }

  SigCal.gears(gears_file);
}