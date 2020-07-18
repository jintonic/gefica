#ifndef GeFiCa_XYZ
#define GeFiCa_XYZ

#include "Grid.h"
namespace GeFiCa { class XYZ; class SquarePointContact; }

/**
 * 3D Cartesian coordinates.
 */
class GeFiCa::XYZ : public GeFiCa::Grid
{
   public:
      XYZ(size_t n1=50, size_t n2=50, size_t n3=50) : Grid(n1, n2, n3) {
         fName="xyz", fTitle="3D Cartecian coordinates"; }
      void SetupWith(Detector &detector);
      double GetC();
   protected:
      void OverRelaxAt(size_t idx); 
      void GetInfoFrom(SquarePointContact &detector);
      void CalculateE();
      ClassDef(XYZ,1);
};
#endif
