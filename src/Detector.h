#ifndef GeFiCa_Detector
#define GeFiCa_Detector

namespace GeFiCa { class Crystal; class Detector; }

/**
 * Crystal properties.
 */
class GeFiCa::Crystal
{
   public:
      /**
       * Height of crystal.
       * Since the derived classes can be used to configure not only 1D but
       * also 2D and 3D grids, this property exists in all of them.
       */
      double Height;
      double TopImpurity; ///< net impurity concentration at top of crystal
      double BottomImpurity; ///< net impurity concentration at bottom of crystal
      /**
       * \return net impurity concentration at \param height.
       */
      double GetImpurity(double height)
      { return (TopImpurity-BottomImpurity)*height/Height+BottomImpurity; }
      void SetAverageImpurity(double impurity)
      { TopImpurity=impurity; BottomImpurity=impurity; }
};

#include "Grid.h"

/**
 * Detector & crystal properties.
 */
class GeFiCa::Detector : public GeFiCa::Crystal
{
   public:
      std::vector<double> Bias; ///< bias on electrodes
      virtual ~Detector() {};

      virtual void Configure(Grid& grid)=0; //< Configure \param grid.
};
#endif
