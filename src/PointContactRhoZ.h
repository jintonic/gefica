#ifndef GEFICA_POINTCONTACTRHOZ_H
#define GEFICA_POINTCONTACTRHOZ_H

#include "PointContactDZ.h"

namespace GeFiCa { class PointContactRhoZ; }

/**
 * Grid setup for 2D point contact detectors.
 * The grid is setup in [0, Radius] and [0, Height].
 */
class GeFiCa::PointContactRhoZ : public GeFiCa::PointContactDZ
{
   public:
      /**
       * Default constructor.
       */
      PointContactRhoZ(int nr=101, int nz=101,
            const char *name="pcrz",
            const char *title="2D point contact detector")
         : PointContactDZ(nr, nz, name, title) {};

      ClassDef(PointContactRhoZ,1);

   protected:
      void Initialize();
};

#endif
