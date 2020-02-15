#ifndef GeFiCa_XYZ
#define GeFiCa_XYZ

#include "Grid.h"
#include "SquarePointContact.h"

namespace GeFiCa { class XYZ; class PointContact; }

/**
 ** 3D coordinates.
 **/
class GeFiCa::XYZ : public GeFiCa::Grid
{
   public:
   /**
   ** Default constructor.
   **/
      XYZ(int N1=50, int N2=50, int N3=50) :Grid(N1,N2,N3) {fName="xyz",
         fTitle="3D coordinates";}
      void SetupWith(Detector &detector);
      double GetC();
   protected:
      void OverRelaxAt(size_t idx); 
      void GetInfoFrom(SquarePointContact &detector);
      void CalculateE();
      void GeneralSetup(SquarePointContact &idetector);
      ClassDef(XYZ,1);
};
#endif
