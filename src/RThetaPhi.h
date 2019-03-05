#ifndef GeFiCa_RTHETAPHI_H
#define GeFiCa_RTHETAPHI_H

#include "XYZ.h"

namespace GeFiCa { class RThetaPhi; }

/**
 * 3D spherical coordinates.
 */
class GeFiCa::RThetaPhi : public GeFiCa::XYZ
{
   public:
      /**
       * Default constructor.
       */
      RThetaPhi(int n_r=101, int n_theta=181, int n_phi=10,
            const char *name="rtp",
            const char *title="3D spherical coordinates")
         : XYZ(n_r, n_theta, n_phi*2, name, title) {};

      ClassDef(RThetaPhi,1);

   protected:
      virtual void DoSOR2(int idx);

      virtual double GetData(double x, double y, double z, double *data);
};
#endif

