#include "Units.h"
#include "SignalCalculator.h"
using namespace GeFiCa; // ?????

#define HOLE_CHARGE 1.0
#define ELECTRON_CHARGE -1.0
#define REF_TEMP 77.0

// class TFile;
#include "TFile.h"
#include "Riostream.h"

SignalCalculator::SignalCalculator(const char *input="ICPC.root", const char *driftName="drift_vel_tcorr.tab")
{
  TFile *file = new TFile(input);
  grd = (RhoZ*) file->Get("rhoz");

  det = (PointContact*) file->Get("pc");
  N   = grd->GetN();

  SetupVelo(driftName);

  if (diffusion_coef == 0)
    diffusion_coef = v_over_E * 0.67;
  impurity_z0 = det->BottomImpurity;

}

char* SignalCalculator::pt_to_str(char *str, int len, RhoZPhiPoint pt){
  snprintf(str, len, "(%.1f %.1f %.1f)", pt.x(), pt.y(), pt.z());
  return str;
};

/* get_signal
   calculate signal for point pt. Result is placed in signal_out array
   returns -1 if outside crystal
   if signal_out == NULL => no signal is stored
*/
int SignalCalculator::GetSignal(RhoZPhiPoint pt, float *signal_out) {
  static float *signal, *sum, *tmp;
  static int tsteps = 0;
  float w, x, y;
  char  tmpstr[MAX_LINE];
  int   j, k, l, dt, err, comp_f;

  /* first time -- allocate signal and sum arrays */
  if (tsteps != time_steps_calc) {
    tsteps = time_steps_calc;
    if ((signal = (float *) malloc(tsteps*sizeof(*signal))) == NULL ||
  (tmp    = (float *) malloc(tsteps*sizeof(*tmp))) == NULL ||
  (sum    = (float *) malloc(tsteps*sizeof(*sum))) == NULL) {
      Error("GetSignal", "malloc failed in get_signal\n");
      return -1;
    }
  }

  for (j = 0; j < tsteps; j++) signal[j] = 0.0;

  if (det->OutsideDetector(pt)) {
    Error("GetSignal", "Point %s is outside detector!\n", pt_to_str(tmpstr, MAX_LINE, pt));
    return -1;
  }
  Info("GetSignal", "Calculating signal for %s...\n", pt_to_str(tmpstr, MAX_LINE, pt));

  memset(dpath_e, 0, tsteps*sizeof(RhoZPhiPoint));
  memset(dpath_h, 0, tsteps*sizeof(RhoZPhiPoint));

  err = MakeSignal(pt, signal, ELECTRON_CHARGE);
  err = MakeSignal(pt, signal, HOLE_CHARGE);
  /* make_signal returns 0 for success; require hole signal but not electron */

  /* change from current signal to charge signal, i.e.
     each time step contains the summed signals of all previous time steps */
  for (j = 1; j < tsteps; j++) signal[j] += signal[j-1];

  if (signal_out != NULL) {

    if (charge_cloud_size > 0.001 || use_diffusion) {
      /* convolute with a Gaussian to correct for charge cloud size
   and initial velocity
   charge_cloud_size = initial FWHM of charge cloud, in mm,
   NOTE this uses initial velocity of holes only;
   this may not be quite right if electron signal is strong */
      /* difference in time between center and edge of charge cloud */
      dt = (int) (1.5f + charge_cloud_size / (step_time_calc * initial_vel));
      
      if (initial_vel < 0.00001f) dt = 0;
      
      Info("GetSignal", "Initial vel, size, dt = %f mm/ns, %f mm, %d steps\n",
        initial_vel, charge_cloud_size, dt);
      
      if (use_diffusion) 
      {
        dt = (int) (1.5f + final_charge_size / (step_time_calc * final_vel));
        Info("GetSignal", "  Final vel, size, dt = %f mm/ns, %f mm, %d steps\n",
          final_vel, final_charge_size, dt);
      }
      
      if (dt > 1) {
  /* Gaussian */
  w = ((float) dt) / 2.355;
  l = dt/10;     // use l to speed up convolution of waveform with gaussian;
  if (l < 1) {   // instead of using every 1-ns step, use steps of FWHM/10
    l = 1;
  } else if (step_time_out > preamp_tau) {
    if (l > step_time_out/step_time_calc)
      l = step_time_out/step_time_calc;
  } else {
    if (l > preamp_tau/step_time_calc)
      l = preamp_tau/step_time_calc;
  }
  // TELL_CHATTY(">> l: %d\n", l);
  for (j = 0; j < tsteps; j++) {
    sum[j] = 1.0;
    tmp[j] = signal[j];
  }
  for (k = l; k < 2*dt; k+=l) {
    x = ((float) k)/w;
    y = exp(-x*x/2.0);
    for (j = 0; j < tsteps - k; j++){
      sum[j] += y;
      tmp[j] += signal[j+k] * y;
      sum[j+k] += y;
      tmp[j+k] += signal[j] * y;
    }
  }

  for (j = 0; j < tsteps; j++)
    signal[j] = tmp[j]/sum[j];
  
      }
    }

    /* now, compress the signal and place it in the signal_out array;
       truncate the signal if time_steps_calc % ntsteps_out != 0 */
    comp_f = time_steps_calc/ntsteps_out;
    for (j = 0; j < ntsteps_out; j++) signal_out[j] = 0;
    for (j = 0; j < ntsteps_out*comp_f; j++)
      signal_out[j/comp_f] += signal[j]/comp_f;

    /* do RC integration for preamp risetime */
    if (preamp_tau / step_time_out >= 0.1f)
      RcIntegrate(signal_out, signal_out,
      preamp_tau / step_time_out, ntsteps_out);
  }

  /* make_signal returns 0 for success; require hole signal but not electron */
  if (err) return -1;
  return 1;
}

/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
int SignalCalculator::MakeSignal(RhoZPhiPoint pt, float *signal, float q) 
{

  float  wpot, wpot2=0, dwpot=0;
  char   tmpstr[MAX_LINE];
  RhoZPhiPoint  new_pt;
  RhoZPhiVector v, dx;
  float  vel0, vel1 = 0, wpot_old=-1;
  // double diffusion_coeff;
  double repulsion_fact = 0.0, ds2, ds3, dv, ds_dt;
  int    ntsteps, i, t, n, collect2pc, low_field=0, surface_drift=0, stop_drifting = 0;

  new_pt = pt;
  collect2pc = ((q > 0 && impurity_z0 < 0) ||  // holes for p-type 
    (q < 0 && impurity_z0 > 0));   // electrons for n-type
  /*
  if (q > 0) {
    diffusion_coeff = TWO_TIMES_DIFFUSION_COEF_H;
  } else {
    diffusion_coeff = TWO_TIMES_DIFFUSION_COEF_E;
  }
  */
  ntsteps = time_steps_calc;
  for (t = 0; drift_velocity(new_pt, q, &v) >= 0 && !stop_drifting; t++) { 
    if (q > 0) {
      dpath_h[t] = new_pt;
    } else {
      dpath_e[t] = new_pt;
    }
    if (collect2pc) {
      if (t == 0) {
  vel1 = final_vel = initial_vel = sqrt(v.mag2());
  final_charge_size = charge_cloud_size;
  if (use_diffusion) {
    if (final_charge_size < 0.01) final_charge_size = 0.01;
    /* for a spherically symmetric charge cloud, the equivalent
       delta-E at a distance of 1 sigma from the cloud center is
       dE = Q/(4*pi*epsilon*sigma^2)  (Q is charge inside the 3D 1-sigma envelope)
       dE (V/cm) = Q (C) * 1/(4*pi*epsilon) (N m2 / C2) / sigma2 (mm2)
       1 V/m = 1 N/C
       dE (V/cm) = Q (C) * 1/(4*pi*epsilon) (V m / C) / sigma2 (mm2)
       dE (V/cm) = repulsion_fact * FWHM/sigma / (FWHM^2) (mm2), so
       repulsion_fact = (FWHM/sigma)^3 * Q (C) * 1/(4*pi*epsilon) (V m / C) * mm/m * mm/cm
    */
    if (energy > 0.1) {  // set up charge cloud self-repulsion
      repulsion_fact = energy * 0.67*0.67*0.67 / 0.003; // charge in 1 sigma (3D)
      repulsion_fact /= 6.241e18;        // convert to Coulombs
      repulsion_fact *= 9.0e13/16.0;     // 1/(4*pi*epsilon)  (N m2 / C2) * 1e4
      repulsion_fact *= 2.355*2.355*2.355;      // convert FWHM to sigma
    }
  }
  Info("MakeSignal", "initial v: %f (%e %e %e)\n",
        initial_vel, v.x(), v.y(), v.z());
      } else if (use_diffusion) {
  vel0 = vel1;
  vel1 = sqrt(v.mag2());
  final_charge_size *= vel1/vel0;  // effect of acceleration
  // include effects of acceleration and diffusion on cloud size
  dv = repulsion_fact * dv_dE /        // effect of repulsion
          (final_charge_size*final_charge_size);
  // FIXME? this next line could more more fine-grained
  if (dv > 0.05) dv = 0.05;  // on account of drift velocity saturation
  ds_dt = dv + diffusion_coef / final_charge_size;  // effect of diffusion
  if (ds_dt > 0.05 || ds_dt * step_time_calc > 0.1) {
    // nonlinear growth due to small size; need more careful calculation
    Info("MakeSignal", "ds_dt = %.2f; size = %.2f", ds_dt, final_charge_size);
    // ds_dt = 0.05;  // artificially limit nonlinear growth
    ds2 = 2.0 * diffusion_coef * step_time_calc; // increase^2 from diff.
    ds3 = (final_charge_size*final_charge_size *
     (final_charge_size +
      3.0 * dv * step_time_calc));         // FWHM^3 after repulsion
    final_charge_size = sqrt(ds2 + pow(ds3, 0.6667)); 
    Info("MakeSignal", " -> %.2f\n", final_charge_size);
  } else {
    final_charge_size +=  ds_dt * step_time_calc;  // effect of diff. + rep.
  }
      }
    }
    // uncomment me
    // Info("MakeSignal", "pt: (%.4f %.4f %.4f), v: (%e %e %e)",
    // new_pt.x(), new_pt.y(), new_pt.z(), v.x(), v.y(), v.z());
    if (0 && t >= ntsteps - 2) {   // DRC removed (if(0)) Oct 2019; t>ntsteps now dealt with below
      if (collect2pc || wpot > WP_THRESH_ELECTRONS) {
        /* for p-type, this is hole or electron+high wp */
        Info("MakeSignal", "\nExceeded maximum number of time steps (%d)\n", ntsteps);
        low_field = 1;
        // return -1;
      }
      break;
    }
    if (wpotential(new_pt, &wpot) != 0) {
      Info("MakeSignal", "\nCan calculate velocity but not WP at %s!\n",
      pt_to_str(tmpstr, MAX_LINE, new_pt));
      return -1;
    }
    if (wpot < 0.0) wpot = 0.0;

    /* ------------- DCR added Oct 2019: if WP is very small or large, then stop drifting */
    if (!collect2pc &&    wpot < 5.0e-5) stop_drifting = 2;    // drifting to outside
    if (collect2pc && 1.0-wpot < 5.0e-5) stop_drifting = 3;    // drifting to point contact
    if (t >= time_steps_calc - 2) stop_drifting = 1;    // have run out of time...

    if (t > 0) signal[t] += q*(wpot - wpot_old);
    // FIXME! Hack added by DCR to deal with undepleted point contact
    if (wpot >= 0.999 && (wpot - wpot_old) < 0.0002) {
      low_field = 1;
      break;
    }
    wpot_old = wpot;

    dx = v * step_time_calc;
    if (surface_drift && dx.z() < 0) { // currently SetX and SetY are brocken.
      // Hmmm... should the default be zero or one?
      dx.SetXYZ(dx.x() * surface_drift_vel_factor, 
        dx.y() * surface_drift_vel_factor, dx.z());
    }
    new_pt = new_pt + dx;
    // q = charge_trapping(dx, q); //FIXME

    // look for charges on passivated surface of a PPC detector
    if (new_pt.z() < 0 &&                    // at or below surface, and
        (det->WrapAroundR <= det->PointContactR ||  // this is a PPC detector
         new_pt.x()*new_pt.x() + new_pt.y()*new_pt.y() < // or point is inside wrap-around
         det->WrapAroundR * det->WrapAroundR)) {
      Info("MakeSignal", " -> Passivated surface! q = %.2f  r = %.2f\n",
                  q, sqrt(new_pt.x()*new_pt.x() + new_pt.y()*new_pt.y()));
      //break;
      surface_drift = 1;
      new_pt.SetZ(0);
    }

  }

  if (t == 0) {
    Info("MakeSignal", "The starting point %s is outside the WP or field.\n",
    pt_to_str(tmpstr, MAX_LINE, pt));
    return -1;
  }

  if (low_field) {
    Info("MakeSignal", "Low field near point contact; this may or may not be a problem.\n");
  } else {
    Info("MakeSignal", "Drifted to edge of WP or field grid, point: %s q: %.2f\n", 
    pt_to_str(tmpstr, MAX_LINE, new_pt), q);
  }
  if (!low_field && stop_drifting<2) {
    /* figure out how much we must drift to get to the crystal boundary */
    for (n = 0; n+t < ntsteps && !det->OutsideDetector(new_pt); n++){
      new_pt = new_pt + dx;
      if (q > 0) dpath_h[t+n] = new_pt;
      else dpath_e[t+n] = new_pt;
    }
    if (n == 0) n = 1; /* always drift at least one more step */
    // TELL_CHATTY(
    Info("MakeSignal", "q: %.1f t: %d n: %d ((%.2f %.2f %.2f)=>(%.2f %.2f %.2f))\n", 
    q, t, n, pt.x(), pt.y(), pt.z(), new_pt.x(), new_pt.y(), new_pt.z());

    if (n + t >= ntsteps){
      n = ntsteps - t;
      if (q > 0 || wpot > WP_THRESH_ELECTRONS) { /* hole or electron+high wp */
        Info("MakeSignal", "Exceeded maximum number of time steps (%d)\n", ntsteps);
        /* check WP to see if we have produced most of the signal */
        if ((wpot < 0.95 || wpot > 0.05) &&
            wpotential(new_pt, &wpot2) != 0) {
          Info("MakeSignal", "Cannot finish drifting to make at least 95 percent of signal.\n");
          return -1;  /* FIXME: could this be improved? */
        }
        /* drift to new_pt and wpot2 */
        dwpot = (wpot2 - wpot)/n;
      }
    } else {
      /* make WP go gradually to 1 or 0 */
      if (wpot > 0.3) {
        dwpot = (1.0 - wpot)/n;
      } else {
        dwpot = - wpot/n;
      }
    }

    /* now drift the final n steps */
    dx = v * step_time_calc;
    if (new_pt.z() > 0) {               // charges NOT on passivated surface
      for (i = 0; i < n; i++){
        signal[i+t] += q * dwpot;
        // q = charge_trapping(dx, q); //FIXME
      }
    }
  }
  // Info("MakeSignal", "q:%.2f pt: %s\n", q, pt_to_str(tmpstr, MAX_LINE, pt));
  if (q > 0) final_vel = sqrt(v.mag2());

  return 0;
}

//FIXME -- placeholder function. Even parameter list is dodgy
/*
static float charge_trapping(vector dx, float q){
  return q;
}
*/

int SignalCalculator::RcIntegrate(float *s_in, float *s_out, float tau, int time_steps){
  int   j;
  float s_in_old, s;  /* DCR: added so that it's okay to
       call this function with s_out == s_in */
  
  if (tau < 1.0f) {
    for (j = time_steps-1; j > 0; j--) s_out[j] = s_in[j-1];
    s_out[0] = 0.0;
  } else {
    s_in_old = s_in[0];
    s_out[0] = 0.0;
    for (j = 1; j < time_steps; j++) {
      s = s_out[j-1] + (s_in_old - s_out[j-1])/tau;
      s_in_old = s_in[j];
      s_out[j] = s;
    }
  }
  return 0;
}

/* wpotential
   gives (interpolated) weighting potential at point pt, stored in wp
   returns 0 for success, 1 on failure
*/
int SignalCalculator::wpotential(RhoZPhiPoint pt, float *wp)
{
   double weights[4];
   int    idxs[4];

   if (efield_exists(pt, idxs) < 0) return 1;
   grid_weights(pt, idxs, weights);
   *wp = 0.0;
   for (int i = 0; i < 4; i++){
      *wp += weights[i] * (grd->Wp[idxs[i]]);
   }

  return 0;

}


int SignalCalculator::SetupVelo(const char* drift_name)
{
  int vlook_sz = v_lookup_len;

  char  line[MAX_LINE], *c;
  FILE  *fp;
  int   i; //, v_lookup_len;
  struct velocity_lookup tmp, v, v0;
  float sumb_e, sumc_e, sumb_h, sumc_h;

  double be=1.3e7, bh=1.2e7, thetae=200.0, thetah=200.0;  // parameters for temperature correction
  double pwre=-1.680, pwrh=-2.398, mue=5.66e7, muh=1.63e9; //     adopted for Ge   DCR Feb 2015
  double mu_0_1, mu_0_2, v_s_1, v_s_2, E_c_1, E_c_2, e, f;

  if (vlook_sz == 0) {
    vlook_sz = 10; }

  if ((fp = fopen(drift_name, "r")) == NULL){
    Error("SetupVelo", "failed to open velocity lookup table file: '%s'\n", drift_name);
    return -1;
  }
  line[0] = '#';
  c = line;
  while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
  if (c == NULL) {
    Error("SetupVelo", "Failed to read velocity lookup table from file: %s\n", drift_name);
    fclose(fp);
    return -1;
  }
  Info("SetupVelo", "Drift velocity table:\n"
         "  e          e100    e110    e111    h100    h110    h111\n");   
  for (v_lookup_len = 0; ;v_lookup_len++){
    if (v_lookup_len == vlook_sz - 1){
      vlook_sz += 10;

    }
    if (sscanf(&line[0], "%f %f %f %f %f %f %f", 
      &tmp.e, &tmp.e100, &tmp.e110, &tmp.e111, 
      &tmp.h100, &tmp.h110, &tmp.h111) != 7){
      break;
    }
    v_lookup.push_back(tmp);

    Info("SetupVelo", "print v_lookup: %10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
      tmp.e, tmp.e100, tmp.e110, tmp.e111, tmp.h100, tmp.h110,tmp.h111);
    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0' ||
       line[0] == '\n' || line[0] == '\r') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
    if (line[0] == 'e' || line[0] == 'h') break; /* no more velocities data;
                      now reading temp correction data */
  }

  /* check for and decode temperature correction parameters */
  while (line[0] == 'e' || line[0] == 'h') {
    if (line[0] == 'e' &&
   sscanf(line+2, "%lf %lf %lf %lf", 
          &mue, &pwre, &be, &thetae) != 4) break;//asume EOF
    if (line[0] == 'h' &&
   sscanf(line+2, "%lf %lf %lf %lf", 
          &muh, &pwrh, &bh, &thetah) != 4) break;//asume EOF
    if (line[0] == 'e')
      Info("SetupVelo", "electrons: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
        mue, pwre, be, thetae);
    if (line[0] == 'h')
      Info("SetupVelo", "    holes: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
        muh, pwrh, bh, thetah);

    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
  }

  if (v_lookup_len == 0){
    Error("Failed to read velocity lookup table from file: %s\n", drift_name);
    return -1;
  }  
  v_lookup_len++;
  if (vlook_sz != v_lookup_len){
    vlook_sz = v_lookup_len;
  }
  Info("SetupVelo", "Drift velocity table has %d rows of data\n", v_lookup_len);
  fclose(fp);

  /*
    apply temperature dependence to mobilities;
    see drift_velocities.doc and tempdep.c
    The drift velocity reduces at higher temperature due to the increasing of
    scattering with the lattice vibration. We used a model by M. Ali Omar and
    L. Reggiani (Solid-State Electronics Vol. 30, No. 12 (1987) 1351) to
    calculate the temperature dependence.
  */
  /* electrons */
  Info("SetupVelo", "Adjusting mobilities for temperature, from %.1f to %.1f\n", REF_TEMP, xtal_temp);
  Info("SetupVelo", "Index  field  vel_factor\n");
  mu_0_1 = mue * pow(REF_TEMP, pwre);
  v_s_1 = be * sqrt(tanh(0.5 * thetae / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = mue * pow(xtal_temp, pwre);
  v_s_2 = be * sqrt(tanh(0.5 * thetae / xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].e100 *= f;
    v_lookup[i].e110 *= f;
    v_lookup[i].e111 *= f;
    Info("SetupVelo", "%2d %5.0f %f\n", i, e, f);
  }

  /* holes */
  mu_0_1 = muh * pow(REF_TEMP, pwrh);
  v_s_1 = bh * sqrt(tanh(0.5 * thetah / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = muh * pow(xtal_temp, pwrh);
  v_s_2 = bh * sqrt(tanh(0.5 * thetah / xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].h100 *= f;
    v_lookup[i].h110 *= f;
    v_lookup[i].h111 *= f;
    Info("SetupVelo", "%2d %5.0f %f\n", i, e, f);
  }
  /* end of temperature correction */

  for (i = 0; i < vlook_sz; i++){
    v = v_lookup[i];
    v_lookup[i].ea =  0.5 * v.e100 -  4 * v.e110 +  4.5 * v.e111;
    v_lookup[i].eb = -2.5 * v.e100 + 16 * v.e110 - 13.5 * v.e111;
    v_lookup[i].ec =  3.0 * v.e100 - 12 * v.e110 +  9.0 * v.e111;
    v_lookup[i].ha =  0.5 * v.h100 -  4 * v.h110 +  4.5 * v.h111;
    v_lookup[i].hb = -2.5 * v.h100 + 16 * v.h110 - 13.5 * v.h111;
    v_lookup[i].hc =  3.0 * v.h100 - 12 * v.h110 +  9.0 * v.h111;
  }

  v_lookup[0].ebp = v_lookup[0].ecp = v_lookup[0].hbp = v_lookup[0].hcp = 0.0;
  sumb_e = sumc_e = sumb_h = sumc_h = 0.0;

  for (i = 1; i < vlook_sz; i++){
    v0 = v_lookup[i-1];
    v = v_lookup[i];
    sumb_e += (v.e - v0.e)*(v0.eb+v.eb)/2;
    sumc_e += (v.e - v0.e)*(v0.ec+v.ec)/2;
    sumb_h += (v.e - v0.e)*(v0.hb+v.hb)/2;
    sumc_h += (v.e - v0.e)*(v0.hc+v.hc)/2;
    v_lookup[i].ebp = sumb_e/v.e;
    v_lookup[i].ecp = sumc_e/v.e;
    v_lookup[i].hbp = sumb_h/v.e;
    v_lookup[i].hcp = sumc_h/v.e;
  }

  return 0;

}

// efield_exist + nearest_field_grid_index
int SignalCalculator::efield_exists(RhoZPhiPoint pt, int *idxs){
  // cyl_int_pt ipt;
  char ptstr[MAX_LINE];
  sprintf(ptstr, "(r,z) = (%.1f, %.1f)", pt.rho(), pt.z());

  PointContact& pc = (PointContact&) *det;
   if(pc.OutsideDetector(pt))
      Error("efield_exists", "point %s is outside crystal\n", ptstr);

   double min_c1 = *std::min_element(grd->C1.begin(), grd->C1.end());
   int idx1 = (pt.rho() - min_c1) / grd->dC1p[1];
   int idx2 = pt.z()              / grd->dC2p[1];

   //              N2
   //       ---------------
   //      |* * * * * * * *
   //      |* * * * * * * *
   // idx1 |* * * * * * * *
   //      |* * * * * * * *
   //      |* * * * * 
   //       ---------
   //         idx2

   // int *idxs;
   idxs[0] = (idx2 + 1) * grd->N1 + idx1;     // Lower left
   idxs[1] = (idx2 + 1) * grd->N1 + idx1 + 1; // Upper left
   idxs[2] = (idx2)     * grd->N1 + idx1;     // Lower right
   idxs[3] = (idx2)     * grd->N1 + idx1 + 1; // Upper right

   if (idx1 +1 >= grd->N1 || idx2 +1 >= grd->N2){
      Error("efield_exists", "point %s is outside wp table\n", ptstr);
      return -1;
   }

   for (int i = 0; i < 4 ; i++){
      if (grd->Et[idxs[i]] == 0.0) {
         Error("efield_exists", "point %s has no efield\n", ptstr);
         return 0;
      }
   }

  return 1;
}

int SignalCalculator::grid_weights(RhoZPhiPoint pt, int* idxs, double* weights)
{
   RhoZPhiPoint ll(grd->C1[idxs[0]], grd->C2[idxs[0]], 0); // grd->C3[idxs[0]]);
   RhoZPhiPoint ul(grd->C1[idxs[1]], grd->C2[idxs[1]], 0); // grd->C3[idxs[1]]);
   RhoZPhiPoint lr(grd->C1[idxs[2]], grd->C2[idxs[2]], 0); // grd->C3[idxs[2]]);
   RhoZPhiPoint ur(grd->C1[idxs[3]], grd->C2[idxs[3]], 0); // grd->C3[idxs[3]]);

   RhoZPhiVector llv = ll - pt;
   RhoZPhiVector ulv = ul - pt;
   RhoZPhiVector lrv = lr - pt;
   RhoZPhiVector urv = ur - pt;

   weights[0] = llv.rho() * llv.z();
   weights[1] = ulv.rho() * ulv.z();
   weights[2] = lrv.rho() * lrv.z();
   weights[3] = urv.rho() * urv.z();

   return 0;

}

/* Find (interpolated or extrapolated) electric field for this point */
RhoZPhiVector SignalCalculator::efield(RhoZPhiPoint pt, int* idxs){
  double w[4];
  float er, ez, ephi;
  int    i; //, j;

  er = ez = ephi = 0;

  grid_weights(pt, idxs, w);
  for (i = 0; i < 4; i++){
      er += grd->E1[idxs[i]] * w[i];
      ez += grd->E2[idxs[i]] * w[i];
  }

  RhoZPhiVector e(er, ez, ephi);

  return e;
}

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, 1 on success but extrapolation was necessary,
   and -1 for failure
   anisotropic drift: crystal axes are assumed to be (x,y,z)
*/
int SignalCalculator::drift_velocity(RhoZPhiPoint pt, float q, RhoZPhiVector *velo){
  int   i, sign;
  float f, a, b, c;
  float bp, cp, en4, en6;
  struct velocity_lookup v_lookup1, v_lookup2;

   int idxs[4];
   float abs_e;
   RhoZPhiVector e_normal;

   if (efield_exists(pt, idxs) < 0) return -1;

   RhoZPhiVector e = efield(pt, idxs);

   abs_e = sqrt(e.mag2());
   e_normal.SetRho(e.rho() / abs_e);
   e_normal.SetZ  (e.z()   / abs_e);

   if (pt.rho() > 0.001) { // currently SetX and SetY are brocken.
      e_normal.SetXYZ(e_normal.rho() * pt.x() / pt.rho(),
        e_normal.rho() * pt.y() / pt.rho(), e_normal.z());
   } else {
      e_normal.SetXYZ(0, 0, e_normal.z());

   }

  /* find location in table to interpolate from */
   for (i = 0; i < v_lookup_len - 2 && abs_e > v_lookup[i+1].e; i++){};
      v_lookup1 = v_lookup[i];
      v_lookup2 = v_lookup[i+1];

      f = (abs_e - v_lookup1.e)/(v_lookup2.e - v_lookup1.e);
      if (q > 0){
         a = (v_lookup2.ha - v_lookup1.ha)*f+v_lookup1.ha;
         b = (v_lookup2.hb - v_lookup1.hb)*f+v_lookup1.hb;
         c = (v_lookup2.hc - v_lookup1.hc)*f+v_lookup1.hc;
         bp = (v_lookup2.hbp - v_lookup1.hbp)*f+v_lookup1.hbp;
         cp = (v_lookup2.hcp - v_lookup1.hcp)*f+v_lookup1.hcp;
         dv_dE = (v_lookup2.h100 - v_lookup1.h100)/(v_lookup2.e - v_lookup1.e);
      } else {
         a = (v_lookup2.ea - v_lookup1.ea)*f+v_lookup1.ea;
         b = (v_lookup2.eb - v_lookup1.eb)*f+v_lookup1.eb;
         c = (v_lookup2.ec - v_lookup1.ec)*f+v_lookup1.ec;
         bp = (v_lookup2.ebp - v_lookup1.ebp)*f+v_lookup1.ebp;
         cp = (v_lookup2.ecp - v_lookup1.ecp)*f+v_lookup1.ecp;
         dv_dE = (v_lookup2.e100 - v_lookup1.e100)/(v_lookup2.e - v_lookup1.e);
   }

  /* velocity can vary from the direction of the el. field
     due to effect of crystal axes */
#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  en4 = POW4(e_normal.x()) + POW4(e_normal.y()) + POW4(e_normal.z());
  en6 = POW6(e_normal.x()) + POW6(e_normal.y()) + POW6(e_normal.z());
  float abs_v = a + b*en4 + c*en6;
  sign = (q < 0 ? -1 : 1);
  v_over_E = abs_v / abs_e;
  double xv = sign*e_normal.x()*(abs_v+bp*4*(e_normal.x()*e_normal.x() - en4)
               + cp*6*(POW4(e_normal.x()) - en6));
  double yv = sign*e_normal.y()*(abs_v+bp*4*(e_normal.y()*e_normal.y() - en4)
               + cp*6*(POW4(e_normal.y()) - en6));
  double zv = sign*e_normal.z()*(abs_v+bp*4*(e_normal.z()*e_normal.z() - en4)
               + cp*6*(POW4(e_normal.z()) - en6));
  velo->SetXYZ(xv, yv, zv);
#undef POW4
#undef POW6
  return 0;
}

// attaches signal tree to existing gears output file.
void SignalCalculator::gears(const char *gears_file="pc.root")
{
  TFile *f     = new TFile(gears_file, "UPDATE");
  TTree *tg    = (TTree*) f->Get("t");
  TTree *ts    = new TTree("s", "Signal tree");
  int nentries = (Int_t) tg->GetEntries();

  std::vector<double> *x  = new std::vector<double>();
  std::vector<double> *y  = new std::vector<double>();
  std::vector<double> *z  = new std::vector<double>();
  std::vector<double> *de = new std::vector<double>();

  tg->SetBranchAddress("x",  &x);
  tg->SetBranchAddress("y",  &y);
  tg->SetBranchAddress("z",  &z);
  tg->SetBranchAddress("de", &de);

  int evt, stp;
  std::vector<double> sgl;

  ts->Branch("evt", &evt, "evt/I");
  ts->Branch("stp", &stp, "stp/I");
  ts->Branch("sgl","std::vector<double>", &sgl);

  static float *signal;
  RhoZPhiPoint pt;
  // init
  if ((dpath_e = (RhoZPhiPoint*) malloc(time_steps_calc*sizeof(RhoZPhiPoint))) 
    == NULL ||
      (dpath_h = (RhoZPhiPoint*) malloc(time_steps_calc*sizeof(RhoZPhiPoint))) 
    == NULL ||
      (signal = (float*) malloc(ntsteps_out*sizeof(float))) 
    == NULL ) {
    Error("PrintSignal", "Path malloc failed\n");
    return;
  }

  for (int evt_no=0; evt_no<nentries; evt_no++){
    tg->GetEntry(evt_no);
    if(x->size() > 0){
      for (int stp_no = 0; stp_no < x->size(); stp_no++){
        sgl.clear();
        for (int k=0; k < ntsteps_out; k++) signal[k] = 0;
        pt.SetXYZ((*x)[stp_no], (*y)[stp_no], (*z)[stp_no]);
        energy = (*de)[stp_no];
        if (energy == 0) continue; // nothing happens
        if (GetSignal(pt, signal) < 0) {
          printf("point not in crystal or has no field:"
            " (x = %.1f, y = %.1f, z = %.1f)\n",
            pt.x(), pt.y(), pt.z());
          continue;
        }
        
        evt = evt_no;
        stp = stp_no;
        for (int k=0; k < ntsteps_out; k++) sgl.push_back(signal[k]);
        ts->Fill();
        printf("Event/Step number %d/%d filled.\n", evt_no, stp_no);
      }
    }
    
  }

  f->Write();
  f->Close();
}

int SignalCalculator::PrintSignal(RhoZPhiPoint pt){

  static float *signal;

  // init 
  if ((dpath_e = (RhoZPhiPoint*) malloc(time_steps_calc*sizeof(RhoZPhiPoint))) 
    == NULL ||
      (dpath_h = (RhoZPhiPoint*) malloc(time_steps_calc*sizeof(RhoZPhiPoint))) 
    == NULL ||
      (signal = (float*) malloc(ntsteps_out*sizeof(float))) 
    == NULL ) {
    Error("PrintSignal", "Path malloc failed\n");
    return -1;
  }


  if (GetSignal(pt, signal) < 0) {
    printf("point not in crystal or has no field: "
      "(x = %.1f, y = %.1f, z = %.1f)\n",
      pt.x(), pt.y(), pt.z());
    return 1;
  }
  printf("signal: \n");
  for (int i = 0; i < ntsteps_out; i++){
    printf("%f ", signal[i]);
    if (i%10 == 9) printf("\n");
  }

  printf("\n");
  return 0;
}