#ifndef GEFICA_HEMISPHERIC_H
#define GEFICA_HEMISPHERIC_H

#include "Spherical.h"

namespace GEFICA { class HemiSpheric; }

class GEFICA::HemiSpheric : public GEFICA::Spherical
{
   public :
      double Radius, PointContactR;

      HemiSpheric(int n_rho,int n_theta,int n_phi) : Spherical(n_rho, n_theta, n_phi) {};
      void SetVoltage(double voltage,double r); 

      ClassDef(HemiSpheric,1);
};

#endif
