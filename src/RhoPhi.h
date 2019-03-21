#ifndef GeFiCa_RHOPHI_H
#define GeFiCa_RHOPHI_H

#include "XY.h"

namespace GeFiCa { class RhoPhi; }

/**
 * 2D cylindrical coordinates.
 */
class GeFiCa::RhoPhi : public GeFiCa::XY
{
   public:
      /**
       * Default constructor
       */
      RhoPhi(int N1=101, int N2=101,
            const char *name="rp",
            const char *title="2D cylindrical coordinates")
         : XY(N1, N2, name, title) {}; 

      ClassDef(RhoPhi,1);

   protected:
      void OverRelaxAt(int idx);
      virtual double GetData(double x, double y, double z, double *data);
};
#endif

