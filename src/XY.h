#ifndef GeFiCa_XY_H
#define GeFiCa_XY_H

#include "X.h"
class TF2;

namespace GeFiCa { class XY; }

/**
 * 2D grid for field calculation in Cartesian coordinates.
 */
class GeFiCa::XY : public GeFiCa::X
{
   public:
      int n2; ///< number of grid points along the 2nd coordinate

      XY(int nx=101, int ny=101);
      
      virtual ~XY();

      virtual void SaveField(const char *fout=NULL);
      virtual void LoadField(const char *fin=NULL);
		/**
		* Returns the potential in two dimensions
		*/
      double GetPotential(double x,double y){return GetData(x,y,kPotential);};
      /**
      * Returns the value for the electric field under the first direction
      */
      double GetE1(double x,double y){return GetData(x,y,kE1);};
      /**
      * Returns the value for the electic field under the second direction
      */
      double GetE2(double x,double y){return GetData(x,y,kE2);};
      /**
      * Returns the two dimensional impurity
      */
      double GetImpurity(double x,double y){return GetData(x,y,kImpurity);};
   	/**
   	* Method involved in setting the impurity. 
   	* This is used for a variable impurity level that changes with x
   	*/
      void SetImpurity(TF2 *Im);
		/**
		* This defines the class for the CINT library
		*/
      ClassDef(XY,1);

   protected:
      virtual void SOR2(int idx,bool elec);
      double *fE2; /**< Electric field under the second coordinate (x, r, or rho) direction */
      double *fC2; /**< The location under the second coordinate (x, r, or rho) direction*/
      double *fdC2p; ///< distance between this and next grid points alone C2
      double *fdC2m; ///< distance between this and previous grid points alone C2
      /**
      * Uses a binary search to return the index in two dimensions
      */
      int FindIdx(double tarx,double tary,int ybegin,int yend);
      /**
       * Returns data for various variables. 
       */
      double GetData(double tarx,double tary,EOutput output); 
      void SetStepLength(double steplength1,double steplength2); 
      virtual bool CalculateField(int idx);
};
#endif 

