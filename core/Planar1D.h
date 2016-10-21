#ifndef GeFiCa_PLANAR1D_H
#define GeFiCa_PLANAR1D_H

#include "X.h"

namespace GeFiCa { class Planar1D; }
/**
*1 D planar detector geometry for which the static electronic potential and field are calculated.
*/
class GeFiCa::Planar1D : public GeFiCa::X
{
   public :
      double UpperBound,/**< upper boundary of the planar detector.*/LowerBound;/**< lower boundary of the planar detector.*/ 
      double cathode_voltage/**< Volage of the cathode. */,annode_voltage;/**< Voltage of the anode. */

   public :
   /**
   *The constructor for Planar1D, of given a value for nx, no input is needed.
   */
      Planar1D(int nx=101) : X(nx),
	UpperBound(10),/**< Input upper boundary of the planar detector. */
	LowerBound(1),/**< Input lower boundary of the planar detector. */
	cathode_voltage(2000),/**< Volage of the cathode. */
	annode_voltage(0)/**< Voltage of the anode. */ {}; 
	
	  /**
	  *Calculates the step length for the detector
	  */
      void initialize(); 
	  bool Analytic();
      bool CalculateField(EMethod method=kSOR2);
	  /**
	  *This defines the class Planar1D for the cint dictionary.
	  */
      ClassDef(Planar1D, 1);
};

#endif
