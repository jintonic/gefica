#ifndef GeFiCa_XY_H
#define GeFiCa_XY_H

#include "X.h"

namespace GeFiCa { class XY; }

/**
 * 2D coordinates.
 */
class GeFiCa::XY : public GeFiCa::X
{
   public:
      int n2; ///< number of grid points along the 2nd coordinate

      /**
       * Default constructor.
       */
      XY(int nx=101, int ny=101, const char *name="xy",
            const char *title="2D coordinates");
      virtual ~XY();
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

      ClassDef(XY,1);

   protected:
      void SetStepLength(double steplength1,double steplength2); 
      /**
       * Uses a binary search to return the index in two dimensions
       */
      int FindIdx(double tarx,double tary,int ybegin,int yend);
      /**
       * Returns data for various variables. 
       */
      double GetData(double tarx,double tary,EOutput output); 
      virtual void DoSOR2(int idx);
      virtual bool CalculateField(int idx);
};
#endif 

