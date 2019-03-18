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
      RhoZ(int fN1=101, int fN2=101, const char *name="rhoz",
            const char *title="2D cylindrical coordinates")
         : XY(fN1, fN2, name, title) {};

      double GetC();

      ClassDef(RhoZ,1);

   protected:
      virtual void OverRelaxAt(int idx);
};
#endif 
