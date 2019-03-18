#ifndef GeFiCa_RHOPHIZ_H
#define GeFiCa_RHOPHIZ_H

#include "XYZ.h"

namespace GeFiCa { class RhoPhiZ; }

/**
 * 3D cylindrical coordinates.
 */
class GeFiCa::RhoPhiZ : public GeFiCa::XYZ
{
   public:
      RhoPhiZ(int n_rho, int n_phi, int n_z): XYZ(n_rho,n_phi,n_z) {};

      ClassDef(RhoPhiZ,1);

   protected:
      virtual void OverRelaxAt(int idx); 

      virtual double GetData(double x, double y, double z, double *data);
};

#endif
