#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TNamed.h>

class TF3; class TTree; class TGraph;

/**
 * The only namespace in GeFiCa.
 */
namespace GeFiCa {
   /**
    * Switches to use different calculation methods.
    */
   enum EMethod {
      kAnalytic, ///< Analytic calculation
      kSOR2, ///< Successove over-relaxation method to the 2nd order
      kSOR4, ///< Successove over-relaxation method to the 4th order
      kFEniCS ///< Not implemented yet
   };
   class X;
}

/**
 * 1D coordinate.
 */
class GeFiCa::X : public TNamed 
{
   public:
      double V0;///< voltage of inner/lower electrode
      double V1;///< voltage of outer/higher electrode

      int MaxIterations; ///< maximal iteration to be performed
      double Csor; ///< 1<=Csor<2, used to boost iteration speed
      double Precision; ///< difference between two consecutive iterations
      TGraph *Gsor; ///< graph of current precision VS # of iterations

      /**
       * Default constructor for GeFiCa::X.
       * \param [in] nx number of grid points
       * \param [in] name name of an object of this class saved in a ROOT file
       * \param [in] title description of this class
       */
      X(int nx=101, const char *name="x", const char *title="1D coordinate");
      virtual ~X();
      /**
       * Calculate potential using various methods.
       * Current available methods are kAnalytic and kSOR2.
       */
      bool CalculatePotential(EMethod method=kSOR2);
      /**
       * Check if detector is depleted.
       */
      bool IsDepleted();

      X& operator*=(double);
      X& operator+=(X*);
      /**
       * Set average impurity of the crystal as a single number.
       */
      void SetAverageImpurity(double density)
      { for (int i=0; i<fN; i++) fImpurity[i]=density; }
      /**
       * Set impurity distribution.
       * \param [in] fi is a 3D function, but works for lower dimensions as well.
       */
      void SetImpurity(TF3 *fi) { if (fi) fImpDist = fi; }
 
      double GetV(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,fV); }
      virtual double GetE1(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,fE1); }
      virtual double GetE2(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,fE2); }
      virtual double GetE3(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,fE3); }
      virtual double GetImpurity(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,fImpurity); }
      /**
       * Get detector capacitance.
       * calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
       */
      double GetC();
      /**
       * Create &/or return a TTree with field data.
       * \param [in] createNew is a flag
       * - if false (default), the function returns the point of an existing
       *   tree, or create one if there is none.
       * - if true, the function always create a new tree and delete the old
       *   one if there is one.
       */
      TTree* GetTree(bool createNew=false);
      /**
       * Get number of iterations for SOR to converge.
       */
      int GetNsor();

      ClassDef(X,1);

   protected:
      int fN; ///< total number of grid points (fN = fN1 in 1D case)
      int fN1; ///< number of grid points along the 1st coordinate
      int fN2; ///< number of grid points along the 2nd coordinate
      int fN3; ///< number of grid points along the 3nd coordinate

      double *fV; ///< [fN] electric potential
      double *fE1; ///< [fN] electric field along the 1st coordinate
      double *fE2; ///< [fN] electric field along the 2nd coordinate
      double *fE3; ///< [fN] electric field along the 3rd coordinate
      double *fC1; ///< [fN] the 1st coordinate
      double *fC2; ///< [fN] the 2nd coordinate
      double *fC3; ///< [fN] the 3rd coordinate

      double *fdC1p; ///< [fN] step length to next grid point alone C1
      double *fdC1m; ///< [fN] step length to previous grid point alone C1
      double *fdC2p; ///< [fN] step length to next grid point along C2
      double *fdC2m; ///< [fN] step length to previous grid point along C2
      double *fdC3p; ///< [fN] step length to next grid point alone C3
      double *fdC3m; ///< [fN] step length to previous grid point alone C3

      double *fImpurity; ///< [fN] net impurity concentration (Nacceptor-Ndonor)
      bool *fIsFixed; ///< [fN] true if field values are fixed
      bool *fIsDepleted; ///< [fN] true if a grid point is depleted
      TTree *fTree; ///<! ROOT tree to assist field visualization
      TF3 *fImpDist; ///<! Impurity distribution
 
      /**
       * Setup and initialize grid.
       */
      virtual void Initialize() { InitializeGrid(); SetGridImpurity(); }
      virtual void InitializeGrid() {};
      virtual void SetGridImpurity();
      /**
       * find surrounding index and return in int array
       */
      virtual int* FindSurroundingMatrix(int idx);
      /**
       * Initialize fC1, fdC1p, fdC1m, fIsFixed
       */
      void SetStepLength(double stepLength);
      /**
       * Uses a binary search to return the index .
       */
      int FindIdx(double tarx,int begin,int end);
      /**
       * Calculates the field using an analytical method.
       */
      virtual bool Analytic();
      /**
       * Interpolate grid data at a given location.
       */
      double GetData(double x, double y, double z, double *data); 
      /**
       * Calculate fields at @param idx using SOR2.
       */
      virtual void DoSOR2(int idx); 
      /**
       * Calculate electric field after CalculatePotential.
       */
      virtual bool CalculateField(int idx);

      int GetIdxOfMaxV(); ///< Get index of the grid point with max potential
      int GetIdxOfMinV(); ///< Get index of the grid point with min potential
};
#endif

