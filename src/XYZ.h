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
      XYZ(int n1=101, int n2=101, int n3=101, const char *name="xyz",
            const char *title="3D coordinates"); ///< Default constructor
      virtual ~XYZ();

      virtual void SetImpurity(TF3 *Im); ///< Set impurity distribution

      virtual void SaveField(const char *fout);
      virtual void LoadField(const char *fin);

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
      double *fE3; ///< [n] electric field along the 3rd coordinate
      double *fC3; ///< [n] the 3rd coordinate
      double *fdC3p; ///< [n] step length to next grid point alone C3
      double *fdC3m; ///< [n] step length to previous grid point alone C3

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

