#ifndef GeFiCa_XY
#define GeFiCa_XY

#include "Grid.h"

namespace GeFiCa { class XY; }
/**
 * 2D Cartesian coordinates.
 */
class GeFiCa::XY : public GeFiCa::Grid
{
   public:
      XY(size_t nx=101, size_t ny=101) :
         Grid(nx,ny) { SetName("xy"); SetTitle("2D Cartesian coordinates"); }
      virtual ~XY();
      /**
       * Get an electric field line originated from (\param x, \param y).
       * If \param positive, propagate along E direction;
       * else propagate against E direction.
       */
      FieldLine* GetFieldLineFrom(double x, double y, bool positive=true);

      ClassDef(XY,1);

   protected:
      double GetData(const std::vector<double> &data,
            double x, double y, double z) const; 
      virtual void OverRelaxAt(size_t idx);
      virtual bool CalculateField(size_t idx);
};
#endif 

