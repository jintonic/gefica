#ifndef GeFiCa_XY_H
#define GeFiCa_XY_H

class TGraph; class TMultiGraph;

#include "X.h"

namespace GeFiCa { class XY; }

/**
 * 2D coordinates.
 */
class GeFiCa::XY : public GeFiCa::X
{
   public:
      /**
       * Default constructor.
       */
      XY(int nx=101, int ny=101, const char *name="xy",
            const char *title="2D coordinates");
      virtual ~XY();

      /**
       * Get an electric field line originated from (\param x, \param y).
       */
      TGraph* GetFieldLineFrom(double x, double y);

      ClassDef(XY,1);

   protected:
      TMultiGraph *fEgraphs; ///< graphs of electric field lines

      void SetStepLength(double steplength1,double steplength2); 
      int FindIdx(double tarx,double tary,int ybegin,int yend);
      virtual double GetData(double x,double y, double z, double *data); 
      virtual void DoSOR2(int idx);
      virtual bool CalculateField(int idx);
};
#endif 

