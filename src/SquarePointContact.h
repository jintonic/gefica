#ifndef GeFiCa_SquarePointContact
#define GeFiCa_SquarePointContact

#include "Detector.h"
namespace GeFiCa { class SquarePointContact; }
/**
 * Configuration of Squre point contact detectors.
 */
class GeFiCa::SquarePointContact : public Detector
{
   public:
      double Width; ///< Widthi(x) of crystal
      double Length; ///< Length(y) of crystal

      double PointContactW; ///< Width of point contact
      double PointContactL; ///< Length of point contact
      double PointContactH; ///< Height of point contact

      double WrapAroundW; ///< Inner radius of outer contact 

      SquarePointContact(const char *name="spc",
            const char *title="squre point-contact detector");
      void CheckConfigurations();
      void Draw(Option_t* option="side");
      ClassDef(SquarePointContact,1);
};
#endif

