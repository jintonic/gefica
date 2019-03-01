#ifndef GeFiCa_R_H
#define GeFiCa_R_H

#include "X.h"

namespace GeFiCa { class R; }

/**
 * 1D spherical coordinate.
 */
class GeFiCa::R : public X 
{
   public:
      /**
       * Default constructor
       */
      R(int n=101, const char *name="r",
            const char *title="1D spherical coordinate")
         : X(n, name, title) {};

      ClassDef(R,1);

   protected:
      virtual void DoSOR2(int idx);
      virtual void DoSOR4(int idx);
};
#endif

