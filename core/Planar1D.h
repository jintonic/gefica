#ifndef GeFiCa_PLANAR1D_H
#define GeFiCa_PLANAR1D_H

#include "X.h"

namespace GeFiCa { class Planar1D; }

class GeFiCa::Planar1D : public GeFiCa::X
{
   public :
      double UpperBound,LowerBound; // boundary of the planar detector
      double cathode_voltage,annode_voltage;

   public :
      Planar1D(int nx=101) : X(nx),
	UpperBound(10),LowerBound(1),cathode_voltage(2000),annode_voltage(0) {};
      
      void initialize(); 
      bool Analytic();
      bool CalculateField(EMethod method=kSOR2);
      ClassDef(Planar1D, 1);
};

#endif
