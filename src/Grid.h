#ifndef GeFiCa_Grid
#define GeFiCa_Grid
namespace GeFiCa { class Points; class FieldLine; class Grid; class Detector; }
#include <vector>
/**
 * A group of discrete points.
 */
class GeFiCa::Points
{
   public:
      std::vector<double> C1; ///< the 1st coordinates of the points
      std::vector<double> C2; ///< the 2nd coordinates of the points
      std::vector<double> C3; ///< the 3rd coordinates of the points
      std::vector<double> Vp; ///< potential at each point
      std::vector<double> Et; ///< total electric field strength
      std::vector<double> E1; ///< projection of Et on C1
      std::vector<double> E2; ///< projection of Et on C2
      std::vector<double> E3; ///< projection of Et on C3
      std::vector<double> dC1p; ///< step length to next point alone C1
      std::vector<double> dC1m; ///< step length to previous point alone C1
      std::vector<double> dC2p; ///< step length to next point along C2
      std::vector<double> dC2m; ///< step length to previous point along C2
      std::vector<double> dC3p; ///< step length to next point alone C3
      std::vector<double> dC3m; ///< step length to previous point alone C3
      size_t GetN() { return C1.size(); } ///< total number of points
};
#include <TNamed.h>
class TGraph;
/**
 * Electric field line data.
 * It inherits from TNamed the ability to be retreived from a ROOT file by its
 * name. One can get either a TGraph or C1, C2... from it for plotting.
 */
class GeFiCa::FieldLine : public Points, public TNamed
{
   public:
      FieldLine() : Points(), TNamed("fl", "electric field line"), fGl(0) {};
      TGraph* GetGraph();
   private:
      TGraph* fGl; ///< a TGraph object to draw the field line
      ClassDef(FieldLine,1);
};
class TTree;
/**
 * Data structure of a electric field grid.
 * It inherits the flat data structure from Points instead of using the class
 * to save typing. Compare `C1[i]` with `point.C1[i]`.
 */
class GeFiCa::Grid : public Points
{
   public:
      std::vector<double> Src; ///< -(net impurity concentration)x|Qe|/epsilon
      size_t N1; ///< number of points along the 1st coordinate
      size_t N2; ///< number of points along the 2nd coordinate
      size_t N3; ///< number of points along the 3rd coordinate
      size_t MaxIterations; ///< maximal iterations of SOR to be performed
      double RelaxationFactor; ///< within (0,2), used to speed up convergence
      double Precision; ///< difference between two consecutive SOR iterations
      /**
       * Default constructor.
       * It also defines ROOT drawing style.
       */
      Grid(size_t n1=0, size_t n2=0, size_t n3=0);
      virtual ~Grid() {};
      /**
       * Fix potentials on boundaries based on \param detector geometry.
       * It fills Points data based on \param detector geometry and N1,
       * N2 and/or N3, and raises the flag fIsFixed for points on/outside
       * electrodes. It has to be called before SuccessiveOverRelax().
       */
      virtual void GetBoundaryConditionFrom(Detector &detector);
      /**
       * Successively over-relax potentials on grid points.
       */
      void SuccessiveOverRelax();
      /**
       * Get number of iterations for SOR to converge.
       */
      size_t GetIterations() { return fIterations; }
      /**
       * Solve Poisson's Equation analytically.
       * It only accepts a constant impurity throughout the grid.
       */
      void SolveAnalytically();
      /**
       * Get potential at (c1,c2,c3) by interpolation.
       */
      double GetV(double c1, double c2=0, double c3=0) const
      {return GetData(Vp,c1,c2,c3); }
      double GetE(double c1, double c2=0, double c3=0) const
      { return GetData(Et,c1,c2,c3); }
      double GetE1(double c1, double c2=0, double c3=0) const
      { return GetData(E1,c1,c2,c3); }
      double GetE2(double c1, double c2=0, double c3=0) const
      { return GetData(E2,c1,c2,c3); }
      double GetE3(double c1, double c2=0, double c3=0) const
      { return GetData(E3,c1,c2,c3); }
      /**
       * Get detector capacitance.
       * Calculate C based on \f$CV^2/2 = \epsilon \int E^2/2 dx^3\f$.
       * http://hyperphysics.phy-astr.gsu.edu/hbase/electric/capeng.html
       */
      double GetC();
      /**
       * Create &/or return a TTree with field data.
       * \param [in] createNew is a flag
       * - if false (default), the function returns the point of an existing
       *   tree, or create one if there is none.
       * - if true, the function always create a new tree and delete the old
       *   one if it exists.
       */
      TTree* GetTree(bool createNew=false);
      /**
       * Check if every grid point is depleted.
       */
      bool IsDepleted();
      /**
       * Potentials at all points are multiplied by \param scale.
       */
      Grid& operator*=(double scale);
      /**
       * Potentials of this grid are summed with those of \param other grid.
       */
      Grid& operator+=(Grid& other);
      /**
       * Propogate a field line from (c1,c2,c3).
       */
      virtual FieldLine* GetFieldLineFrom(double c1, double c2, double c3=0)
      { return 0; }
   protected:
      std::vector<bool> fIsFixed; ///< true if field values are fixed
      std::vector<bool> fIsDepleted; ///< true if a grid point is depleted
      TTree* fTree; ///<! ROOT tree to visualize fields
      Detector* fDetector; ///<! Pointer to associated detector object
      size_t fIterations; ///< number of iterations of SOR performed
      /**
       * Over relax potential Vp[\param idx].
       */
      virtual void OverRelaxAt(size_t idx) {};
      /**
       * Get index of point near \param c1 in between \param begin & \param end.
       */
      size_t GetIdxOfPointToTheRightOf(double c1,
            size_t begin, size_t end) const;
      size_t GetIdxOfPointToTheRightOf(double c1, double c2,
            size_t begin, size_t end) const;
      size_t GetIdxOfPointToTheRightOf(double c1, double c2, double c3,
            size_t begin, size_t end) const;
      /**
       * Interpolate grid data at (c1,c2,c3).
       */
      virtual double GetData(const std::vector<double> &data,
            double c1, double c2, double c3) const;
      /**
       * Get index of the grid point with max potential.
       */
      size_t GetIdxOfMaxV();
      /**
       * Get index of the grid point with min potential.
       */
      size_t GetIdxOfMinV();
      /**
       * Calculate Et, E1, E2, E3 from Vp.
       */
      virtual void CalculateE();
      double twopoint(double dataset[2],double tarlocationset,double pointxset[2])const;
      double threepoint(double dataset[3],double tarlocationset[2],double pointxset[3],double pointyset[3])const;
      double fourpoint(double dataset[4],double tarlocationset[2],double pointxset[4],double pointyset[4])const;
};
#endif
