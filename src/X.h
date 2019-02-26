/** @file X.h
 * Header defining base class GeFiCa::X.
 */
#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TNamed.h>

class TF3;
class TTree;

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
   };
   /**
    * Switches to different output.
    */
   enum EOutput {
      kImpurity, ///< crystal impurity concentration
      kPotential, ///< electric potential
      kE1, ///< Electric field conponent in the 1st coordinate
      kE2, ///< Electric field conponent in the 2nd coordinate
      kE3, ///< Electric field conponent in the 3rd coordinate
   };

   class X;
}

/**
 * 1D coordinate.
 */
class GeFiCa::X : public TNamed 
{
   public:
      int n1; ///< number of grid points along the 1st coordinate
      int n; ///< total number of grid points (n = n1 in 1D case)
      double Csor; ///< 1<=Csor<2, used to boost iteration speed
      double Precision; ///< difference between two consecutive iterations
      int MaxIterations; ///< maximal iteration to be performed

      double V0;///< voltage of one electrode
      double V1;///< voltage of the other electrode

      /**
       * Default constructor for GeFiCa::X.
       * \param nx number of grid points
       * \param name name of an object of this class saved in a ROOT file
       * \param title description of this class
       */
      X(int nx=101, const char *name="x", const char *title="1D coordinate");
      virtual ~X();
      /**
       * Calculate potential using various methods.
       * Current available methods are kAnalytic and kSOR2.
       */
      bool CalculatePotential(EMethod method=kSOR2);
      
      bool IsDepleted(); ///< check if the detector is depleted

      X& operator*=(double);
      X& operator+=(X*);
      
      /**
       * Setup and initialize grid.
       */
      virtual void Initialize() {};
      /**
       * find surrounding index and return in int array
       */
      virtual int* FindSurroundingMatrix(int idx);
      /**
       * Set average impurity of the crystal as a single number.
       */
      void SetAverageImpurity(double density)
      { for (int i=0; i<n; i++) fImpurity[i]=density; }
      /**
       * Set impurity that changes with position.
       */
      virtual void SetImpurity(TF3 *fi);
 
      double GetE1(double x){return GetData(x,kE1);};
      double GetE2(double x){return 0;};
      double GetE3(double x){return 0;};
      double GetImpurity(double x){return GetData(x,kImpurity);};
      double GetPotential(double x){return GetData(x,kPotential);};
      /**
       * Get Capacitance, Cdet.
       * calculate Cdet based on CV^2/2 = epsilon int E^2 dx^3 / 2
       */
      double GetCapacitance();
      /**
       * create &/or return a TTree with field data
       */
      virtual TTree* GetTree();

      ClassDef(X,1);

   protected:
      double *fV; ///< [n] electric potential
      double *fE1; ///< [n] electric field along the 1st coordinate
      double *fC1; ///< [n] the 1st coordinate
      double *fdC1p; ///< [n] step length to next grid point alone C1
      double *fdC1m; ///< [n] step length to previous grid point alone C1
      double *fImpurity; ///< [n] net impurity concentration (Nacceptor-Ndonor)
      bool *fIsFixed; ///< [n] true if field values are fixed
      bool *fIsDepleted; ///< [n] true if a grid point is depleted
      TTree *fTree; ///<! ROOT tree to assist field visualization
 
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
       * Returns data of certain EOutput.
       */
      double GetData(double tarx, EOutput output); 
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

