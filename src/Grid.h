#ifndef GeFiCa_Grid
#define GeFiCa_Grid

#include <vector>

#include <TTree.h>

namespace GeFiCa { class Grid; }

/**
 * Define data structure of a grid.
 */
class GeFiCa::Grid
{
   public:
      size_t N1; ///< number of points along 1st coordinate
      size_t N2; ///< number of points along 2nd coordinate
      size_t N3; ///< number of points along 3rd coordinate

      std::vector<double> V; ///< potentials of points in grid
      std::vector<double> E; ///< total E of points in grid
      std::vector<double> E1; ///< E1 of points in grid
      std::vector<double> E2; ///< E2 of points in grid
      std::vector<double> E3; ///< E3 of points in grid
      std::vector<double> C1; ///< 1st coordinates of points in grid
      std::vector<double> C2; ///< 2nd coordinates of points in grid
      std::vector<double> C3; ///< 3rd coordinates of points in grid
      std::vector<double> Src; ///< -impurity*|Qe|/epsilon
      std::vector<double> dC1p; ///< step length to next grid point alone C1
      std::vector<double> dC1m; ///< step length to previous grid point alone C1
      std::vector<double> dC2p; ///< step length to next grid point along C2
      std::vector<double> dC2m; ///< step length to previous grid point along C2
      std::vector<double> dC3p; ///< step length to next grid point alone C3
      std::vector<double> dC3m; ///< step length to previous grid point alone C3

      double GetV(double c1, double c2=0, double c3=0) const
      { return GetData(c1,c2,c3,V); }
      double GetE(double c1, double c2=0, double c3=0) const
      { return GetData(c1,c2,c3,E); }
      double GetE1(double c1, double c2=0, double c3=0) const
      { return GetData(c1,c2,c3,E1); }
      double GetE2(double c1, double c2=0, double c3=0) const
      { return GetData(c1,c2,c3,E2); }
      double GetE3(double c1, double c2=0, double c3=0) const
      { return GetData(c1,c2,c3,E3); }

      double RelaxationFactor; ///< within (0,2)

      Grid(size_t n1=0, size_t n2=0, size_t n3=0) :
         N1(n1), N2(n2), N3(n3), RelaxationFactor(1.95) {};
      virtual ~Grid() {};

      void SuccessiveOverRelax()
      { for (size_t i=0; i<V.size(); i++) OverRelaxAt(i); }

      size_t GetN() { return V.size(); } ///< total number of grid points
      TTree* GetTree();

   protected:
      TTree* fTree; ///<! ROOT tree to visualize fields

      /**
       * Over relax potential V at \param idx.
       */
      virtual void OverRelaxAt(size_t idx)=0;

      size_t GetIdxOfPointNear(double c1, size_t begin, size_t end=0);
      size_t GetIdxOfPointNear(double c1, double c2, size_t begin, size_t end=0);
      size_t GetIdxOfPointNear(double c1, double c2, double c3, size_t begin, size_t end=0);
      virtual double GetData(double c1, double c2, double c3,
            const std::vector<double> &data) const = 0;
};
#endif
