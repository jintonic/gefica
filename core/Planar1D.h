#ifndef GeFiCa_PLANAR1D_H
#define GeFiCa_PLANAR1D_H

#include "X.h"

namespace GeFiCa { class Planar1D; }
/**
 * Field for 1D planar detector.
 *
 * The default number of grids is 101. The default location of the lower
 * boundary of the detector is at 0, the upper boundary is at 1*cm. The default
 * voltage is 2000*volt with the potential at lower boundary to be 0 and upper
 * boundary 2000*volt. If the calculated potential is bent more than 2000*volt,
 * it is likely that the thickness specified is too large that the detector is
 * not depleted with 2000*volt.
 */
class GeFiCa::Planar1D : public GeFiCa::X
{
   public :
      double UpperBound;///< upper boundary of the planar detector
      double LowerBound;///< lower boundary of the planar detector

   public :
      Planar1D(int nx=101) : X(nx), UpperBound(1), LowerBound(0) 
      {}; 

      /**
       * Calculate the step length of the grid.
       *
       * double stepLength=(UpperBound-LowerBound)/(n-1);
       */
      void Initialize(); 

      /**
       * Analytic calculation of 1D field with fixed impurity concentration.
       *
       * In case of fixed impurity, potential(x) = a x^2 + b x + c with
       * boundary conditions:
       *
       * - potential(0) = Vneg,
       * - potential(d) = Vpos,
       *
       * where d = UpperBound - LowerBound. It also obeys Gauss's Law:
       *
       * - d^2 p / dx^2 = -rho/epsilon
       *
       * So, 
       *
       * - a = - rho/2/epsilon
       * - b = (Vpos-Vneg - ad^2)/d
       * - c = Vneg
       */
      bool Analytic();
      bool CalculateField(EMethod method=kSOR2);
      /**
       *This defines the class Planar1D for the cint dictionary.
       */
      ClassDef(Planar1D, 1);
};
#endif

