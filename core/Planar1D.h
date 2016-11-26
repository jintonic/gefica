#ifndef GeFiCa_PLANAR1D_H
#define GeFiCa_PLANAR1D_H

#include "X.h"

namespace GeFiCa { class Planar1D; }
/**
 * 1 D planar detector geometry for which the static electronic potential and
 * field are calculated.
 */
class GeFiCa::Planar1D : public GeFiCa::X
{
   public :
      double UpperBound;///< upper boundary of the planar detector
      double LowerBound;///< lower boundary of the planar detector
      double Vpos;///< Volage of the cathode
      double Vneg;///< Voltage of the anode

   public :
      /**
       * The constructor for Planar1D, of given a value for nx, no input is
       * needed.
       */
      Planar1D(int nx=101) : X(nx), UpperBound(10), LowerBound(0), 
      Vpos(0), Vneg(2000){}; 

      /**
       *Calculates the step length for the detector
       */
      void Initialize(); 
      /**
       * potential(x) = a x^2 + b x + c with boundary conditions:
       * p(LowerBound)=fPotential[0]; p(UpperBound)=fPotential[n-1];
       * d^2 p / dx^2 = rho/epsilon; so
       * a = 
       * b = 
       * c = 
       */
      bool Analytic();
      bool CalculateField(EMethod method=kSOR2);
      /**
       *This defines the class Planar1D for the cint dictionary.
       */
      ClassDef(Planar1D, 1);
};
#endif

