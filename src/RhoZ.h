#ifndef GeFiCa_RhoZ_H
#define GeFiCa_RhoZ_H

#include "XY.h"

namespace GeFiCa { class RhoZ; }

/**
 * 2D cylindrical coordinates.
 */
class GeFiCa::RhoZ : public GeFiCa::XY
{
   public:
      /**
       * Default Constructor.
       */
      RhoZ(int N1=101, int N2=101, const char *name="rhoz",
            const char *title="2D cylindrical coordinates")
         : XY(N1, N2, name, title) {};

      double GetC();

      ClassDef(RhoZ,1);

   protected:
      virtual void OverRelaxAt(int idx);
};
#endif 
