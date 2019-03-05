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
      int n3; ///< number of grid points along the 3nd coordinate

   public:
      /**
       * Default constructor.
       */
      XYZ(int n1=101, int n2=101, int n3=11, const char *name="xyz",
            const char *title="3D coordinates");
      virtual ~XYZ();

      double GetPotential(double x,double y,double z)
      {return GetData(x,y,z,kPotential);};
      virtual double GetE1(double x,double y,double z)
      {return GetData(x,y,z,kE1);};
      virtual double GetE2(double x,double y,double z)
      {return GetData(x,y,z,kE2);};
      virtual double GetE3(double x,double y,double z)
      {return GetData(x,y,z,kE3);};
      virtual double GetImpurity(double x,double y,double z)
      {return GetData(x,y,z,kImpurity);};

      ClassDef(XYZ,1);

   protected:
      virtual double GetData
         (double tarx,double tary,double tarz, EOutput output);
      virtual void SetStepLength
         (double steplength1,double steplength2,double steplength3);
      virtual int FindIdx(double tarx,double tary,
            double tarz,int begin,int end);
      virtual void DoSOR2(int idx); 
      virtual bool CalculateField(int idx);
};
#endif

