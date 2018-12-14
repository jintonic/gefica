#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TObject.h>
#include <string>

class TF1;

namespace GeFiCa { 
   enum EMethod ///< Different methods to calculate fields
   {
      kAnalytic,
      kSOR2, ///< Successove over-relaxation method to the 2nd order
      kSOR4, ///< Successove over-relaxation method to the 4th order
   };
   enum EOutput ///< Different output of the calculate field
   {
      kImpurity,
      kPotential,
      kE1,
      kE2,
      kE3,
   };

   class X;
}

/**
 * 1D grid for field calculation.
 *
 * Successive Over-Relaxation (SOR) method is used to calculate the field in a
 * grid. Please refer to https://mediatum.ub.tum.de/node?id=969435 for detailed
 * description of SOR method. It will have some error compared with actual data
 * but should be close.  Analytic solution is also provided for comparison.
 */
class GeFiCa::X : public TObject 
{
   public:
      int n1; ///< number of grid points along the 1st coordinate
      int n; ///< total number of grid points (n = n1 in 1D case)
      double Csor; ///< 1<=Csor<2, used to boost iteration speed
      double Precision; ///< difference between two consecutive iterations
      int MaxIterations; ///< maximal iteration to be performed

      double V0;///< voltage of one electrode
      double V1;///< voltage of the other electrode

      bool *DepletedData;///<is a grid depleted

   public:
      X(int nx=101 /**< [in] Number of grid points */);
      virtual ~X();
      /**
       * Calculate potential using various methods.
       * Current available methods are kAnalytic and kSOR2.
       */
      bool CalculatePotential(EMethod method=kSOR2);
      
      bool IsDepleted(); ///< check if the detector is depleted
      int Findmax();
      int Findmin();

      X& operator*=(double p);
      X& operator+=(X *anotherfield);
      
      virtual void Initialize() {};
      /**
       * find surrounding index and return in int array
       */
      virtual int* FindSurroundingMatrix(int idx);
      /**
       * This function creates a new TFile and TTree and fills it from data
       * created by X::CalculatePotential.    
       */
      virtual void SaveField(const char *fout);
      /**
       * calculate electric field after load
       */
      virtual void LoadField(const char *fin);
      /**
       * Set average impurity of the crystal as a single number.
       */
      void SetImpurity(double density)
      { for (int i=0; i<n; i++) fImpurity[i]=density; }
      /**
       * Set impurity that changes with x.
       */
      void SetImpurity(TF1 *fi1);
      /**
       * Returns the value for E under the first direction.
       */
      double GetE1(double x){return GetData(x,kE1);};
      double GetE2(double x){return 0;};
      double GetE3(double x){return 0;};
      /**
       * Returns the impurity level.
       */ 
      double GetImpurity(double x){return GetData(x,kImpurity);};
      /**
       * Returns the potential.
       */ 
      double GetPotential(double x){return GetData(x,kPotential);};
      void CopyFrom(const X&);
      double  CalculateCapacitance();
      /**
       *This defines the class for the cint dictionary.
       */
      ClassDef(X,1);

   protected:
      bool * fIsFixed; ///< Is used to check if a value is fixed or if it can be modified. It is usually used to check boundary conditions and find out if you are on the edge.
      bool fIsLoaded; ///< fIsLoaded is used to check if points in the grid have a value or not. If fIsLoaded returns true, the points if the grid have value and you do not need to initialize, if it returns false you do.
      double *fE1; ///< Electric field under the first coordinate (x, r, or rho) direction 
      double  *fPotential; ///< Potential in the referenced grid
      double *fC1; ///< the location under the first coordinate (x, r, or rho) direction
      double *fdC1p; ///< distance between this and next grid points alone C1
      double *fdC1m; ///< distance between this and previous grid points alone C1
      double *fImpurity; ///< Value of the impurity level at a point on the grid
      /**
       * Sets the field step length.
       */
      void SetStepLength(double steplength);
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
       * Uses Second order Successive Over-Relaxation method to calculate the field.
       */
      virtual void SOR2(int idx,bool calculateElectricField); 
      /**
       * Calculate electric field after CalculatePotential.
       */
      virtual bool CalculateField(int idx);
};
#endif

