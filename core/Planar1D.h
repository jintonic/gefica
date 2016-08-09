#ifndef GEFICA_PLANAR1D_H
#define GEFICA_PLANAR1D_H

#include "X.h"

namespace GEFICA { class Planar1D; }

class GEFICA::Planar1D : public GEFICA::X
{
   public :
      double UpperBound,LowerBound; // boundary of the planar detector
      double cathode_voltage,annode_voltage;

   public :
      Planar1D(int nx=101) : X(nx),
	UpperBound(10),LowerBound(1),cathode_voltage(2000),annode_voltage(0) {};
      
      void initialize(); 
      bool Analyic();
      bool CalculateField(EMethod method=kRK2);
      ClassDef(Planar1D, 1);
};

#endif
