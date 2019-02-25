#ifndef GeFiCa_RhoZ_H
#define GeFiCa_RhoZ_H

#include "XY.h"
class TF2;

namespace GeFiCa { class RhoZ; }

/**
 * 2D cylindrical coordinates.
 */
class GeFiCa::RhoZ : public GeFiCa::XY
{
   public:
      RhoZ(int n1=101, int n2=101, const char *name="Rhoz",
            const char *title="2D cylindrical coordinates")
         : XY(n1, n2, name, title) {}; ///< Default constructor

      double GetCapacitance();

      ClassDef(RhoZ,1);

   protected:
      virtual void DoSOR2(int idx);
};
#endif 
