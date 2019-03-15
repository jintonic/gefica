#ifndef GeFiCa_XYZ_H
#define GeFiCa_XYZ_H

#include "XY.h"

namespace GeFiCa { class XYZ; }

/**
 * 3D coordinates.
 */
class GeFiCa::XYZ : public GeFiCa::XY
{
   public:
      /**
       * Default constructor.
       */
      XYZ(int fN1=101, int fN2=101, int fN3=11, const char *name="xyz",
            const char *title="3D coordinates");
      virtual ~XYZ();

      ClassDef(XYZ,1);

   protected:
      virtual double GetData(double x, double y, double z, double *data);
      virtual void SetStepLength
         (double steplength1,double steplength2,double steplength3);
      virtual void DoSOR2(int idx); 
      virtual bool CalculateField(int idx);
};
#endif

