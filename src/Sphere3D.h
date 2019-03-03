#ifndef GeFiCa_SPHERE3D_H
#define GeFiCa_SPHERE3D_H

#include "RThetaPhi.h"

namespace GeFiCa { class Sphere3D;}

/**
 * Grid setup for 3D spherical detectors.
 */
class GeFiCa::Sphere3D: public GeFiCa::RThetaPhi
{
   public :
      double OuterR; ///< outer radius
      double InnerR; ///< inner radius

      Sphere3D(int nr=101, ///< number of grid points along radius
            int nt=181, ///< number of grid points along theta
            int np=10, ///< number of grid points along phi
            const char *name="s3d",
            const char *title="3D sphere detector"); ///< default constructor

      ClassDef(Sphere3D,1);

   protected:
      virtual void Initialize();
};
#endif

