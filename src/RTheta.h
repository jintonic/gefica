#ifndef GeFiCa_RTheta
#define GeFiCa_RTheta

#include "R.h"

namespace GeFiCa { class RTheta; }

/**
 * 2D spherical coordinates.
 */
class GeFiCa::RTheta : public GeFiCa::R
{
   public:
      RTheta(size_t n_r, size_t n_theta);

      ClassDef(RTheta,1);

   protected:
      void OverRelaxAt(size_t idx) {};
};
#endif
