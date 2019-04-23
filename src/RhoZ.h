#ifndef GeFiCa_RhoZ
#define GeFiCa_RhoZ

#include <TNamed.h>
#include "Grid.h"
namespace GeFiCa { class RhoZ; class PointContact; class Segmented; }
/**
 * 2D cylindrical coordinates in rho-z plane.
 */
class GeFiCa::RhoZ : public Grid, public TNamed
{
   public:
      RhoZ(size_t n1=200, size_t n2=201) : Grid(n1, n2),
      TNamed("rhoz", "2D cylindrical coordinates in rho-z") {};

      void GetBoundaryConditionFrom(Detector &detector);
      double GetC();

      ClassDef(RhoZ,1);
   protected:
      void OverRelaxAt(size_t idx); 
      void GetInfoFrom(Segmented &detector) {};
      void GetInfoFrom(PointContact &detector);
      void ReallocateGridPointsNearBoundaries(PointContact &detector);
      void CalculateE();
};
#endif 

