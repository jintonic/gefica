#ifndef GeFiCa_Crystal
#define GeFiCa_Crystal

namespace GeFiCa { class Crystal; }

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
};
#endif
