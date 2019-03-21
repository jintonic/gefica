#ifndef GeFiCa_Detector
#define GeFiCa_Detector

#include "Grid.h"

#include "Crystal.h"

namespace GeFiCa { class Detector; class Grid; }

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
