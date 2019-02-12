//////////////////////
//R
//1D Field under Polar co

#ifndef GeFiCa_R_H
#define GeFiCa_R_H

#include "X.h"

namespace GeFiCa { 
   class R;
}

class GeFiCa::R : public X 
{
   public:
      R(int nx=101): X(nx) { SetName("R"); SetTitle("R"); }

      ClassDef(R,1);

      virtual double GetPotential(double r)
      {return GetData(r,kPotential);};
      virtual double GetE1(double r,double theta,double phi)
      {return GetData(r,kE1);};
      virtual double GetImpurity(double r)
      {return GetData(r,kImpurity);};

   protected:
      virtual void DoSOR2(int idx);
      virtual void DoSOR4(int idx);
};
#endif

