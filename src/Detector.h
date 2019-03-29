#ifndef GeFiCa_Detector
#define GeFiCa_Detector
/**
 * The only namespace in GeFiCa.
 */
namespace GeFiCa { class Crystal; class Detector; }
/**
 * Crystal properties.
 */
class GeFiCa::Crystal
{
   public:
      double Height; ///< height of crystal
      double TopImpurity; ///< net impurity concentration at top of crystal
      double BottomImpurity; ///<net impurity concentration at bottom of crystal
      Crystal(); ///< Default constructor.
      /**
       * Return net impurity concentration at \param height.
       */
      double GetImpurity(double height)
      { return (TopImpurity-BottomImpurity)*height/Height+BottomImpurity; }
      void SetAverageImpurity(double impurity)
      { TopImpurity=impurity; BottomImpurity=impurity; }
};
#include <vector>
#include <TNamed.h>
/**
 * Detector & crystal properties.
 * It provides TClass::ClassName() through ClassDef for
 * Grid::SetBoundaryCondition to know the specific type of detector.
 */
class GeFiCa::Detector : public Crystal, public TNamed
{
   public:
      std::vector<double> Bias; ///< bias on electrodes
      Detector(const char *name="detector", const char *title="detector");
      /**
       * Check if detector configurations make sense.
       */
      void CheckConfigurations();
      ClassDef(Detector, 1);
};
#endif
