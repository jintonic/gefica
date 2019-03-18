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
      RhoPhi(int fN1=101, int fN2=101,
            const char *name="rp",
            const char *title="2D cylindrical coordinates")
         : XY(fN1, fN2, name, title) {}; 

      ClassDef(RhoPhi,1);

   protected:
      void OverRelaxAt(int idx);
      virtual double GetData(double x, double y, double z, double *data);
};
#endif

