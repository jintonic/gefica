#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TNamed.h>

class TF3; class TTree; class TGraph;

/**
 * The only namespace in GeFiCa.
 */
namespace GeFiCa { class Grid; class X; }

/**
 * 1D Cartesian coordinate.
 */
class GeFiCa::X : public GeFiCa::Grid, public TNamed 
{
   public:
      /**
       * Default constructor.
       */
      X(size_t nx=101) : Grid(nx), TName("x", "1D coordinate") {};
      virtual ~X();
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
      { return GetData(x,y,z,V); }
      virtual double GetE1(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,E1); }
      virtual double GetE2(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,E2); }
      virtual double GetE3(double x=0, double y=0, double z=0)
      { return GetData(x,y,z,E3); }
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
      int GetN1() { return N1; }
      int GetN2() { return N2; }
      int GetN3() { return N3; }

      double* GetVs() { return V; }

      double* GetE1s() { return E1; }
      double* GetE2s() { return E2; }
      double* GetE3s() { return E3; }

      double* GetC1s() { return C1; }
      double* GetC2s() { return C2; }
      double* GetC3s() { return C3; }

      ClassDef(X,1);

   protected:
      int fN; ///< total number of grid points (fN = N1 in 1D case)
      int N1; ///< number of grid points along the 1st coordinate
      int N2; ///< number of grid points along the 2nd coordinate
      int N3; ///< number of grid points along the 3nd coordinate

      double *V; ///< [fN] electric potential
      double *E1; ///< [fN] electric field along the 1st coordinate
      double *E2; ///< [fN] electric field along the 2nd coordinate
      double *E3; ///< [fN] electric field along the 3rd coordinate
      double *C1; ///< [fN] the 1st coordinate
      double *C2; ///< [fN] the 2nd coordinate
      double *C3; ///< [fN] the 3rd coordinate

      double *dC1p; ///< [fN] step length to next grid point alone C1
      double *dC1m; ///< [fN] step length to previous grid point alone C1
      double *dC2p; ///< [fN] step length to next grid point along C2
      double *dC2m; ///< [fN] step length to previous grid point along C2
      double *dC3p; ///< [fN] step length to next grid point alone C3
      double *dC3m; ///< [fN] step length to previous grid point alone C3

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
       * Initialize C1, dC1p, dC1m, fIsFixed
       */
      void SetStepLength(double stepLength);
      /**
       * Uses a binary search to return the index .
       */
      int FindIdx(double tarx,int begin=0,int end=-1);
      int FindIdx(double tarx,double tary,int begin=0,int end=-1);
      int FindIdx(double tarx,double tary, double tarz,int begin=0,int end=-1);
      /**
       * Interpolate grid data at a given location.
       */
      virtual double GetData(double x, double y, double z, double *data); 
      /**
       * Calculate fields at @param idx using SOR2.
       */
      virtual void OverRelaxAt(int idx); 
      /**
       * Calculate electric field after SuccessiveOverRelax.
       */
      virtual bool CalculateField(int idx);

      int GetIdxOfMaxV(); ///< Get index of the grid point with max potential
      int GetIdxOfMinV(); ///< Get index of the grid point with min potential
};
#endif

