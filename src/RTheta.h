#ifndef GeFiCa_RTheta
#define GeFiCa_RTheta

#include "Grid.h"
namespace GeFiCa { class RTheta; }
/**
 * 2D spherical coordinates.
 */
class GeFiCa::RTheta : public Grid
{
   public:
      RTheta(size_t n1=101, size_t n2=181) : Grid(n1, n2) {
         fName="rt"; fTitle="2D spherical coordinate"; }
      void SetupWith(Detector &detector);
   protected:
      virtual void OverRelaxAt(size_t idx);
      void CalculateE();
      ClassDef(RTheta,1);
};
#endif

