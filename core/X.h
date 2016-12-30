#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TObject.h>
class TF1;

namespace GeFiCa { 
   enum EMethod ///< Different methods to calculate fields
   {
      kAnalytic,
      kSOR2, ///< Successove over-relaxation method to the 2nd order
      kSOR4, ///< Successove over-relaxation method to the 4th order
   };

   class X;
}

/**
 * 1 D grid for field calculation.
 *
 * Successive Over-Relaxation (SOR) method is used to calculate the field in a
 * grid. Please refer to https://mediatum.ub.tum.de/node?id=969435 for detailed
 * description of SOR method. It will have some error compared with actual data
 * but should be close.  Analytic solution is also provided for comparison.
 */
class GeFiCa::X : public TObject 
{
   public:
      int n1; ///< number of grid along the 1st axis
      int MaxIterations; ///< max one turn Iteration number
      int n; ///< n = n1 total number of grid
      double Csor; ///< boost Iteration speed
      double Precision; ///< X limit
      int t,d;

      double V1;///< Volage of the cathode
      double V0;///< Voltage of the anode

   public:

      /**
       * X is a constructor, if given a number, no input is needed
       */
      X(int nx=101 /**< Number of grid lines. */);

      virtual ~X();
      /**
       * Function for deciding which method to use for the field calculation. 
       * Current methods are kAnalytic and kSOR2.
       */
      virtual bool CalculateField(EMethod method=kSOR2);
      /**
       * This function creates a new TFile and TTree and fills it from data created by X::CalculateField.    
       */
      virtual void SaveField(const char *fout);
      /**
       * calculate electric field after load
       */
      virtual void LoadField(const char *fin);
      /*! \brief Ionizing impurity level method
       * 
       * This function takes an argument for the variable density. It can be used if you consider impurity to be constant.
       * It is usualy in the form of 1e10/cm3 or some variation where cm3 is defined in the namespace GeFiCa.
       */
      virtual void SetImpurity(double density); // If you can consider impurity to e constant
      /**
       * Another Important method involved in setting the impurity. This is used for a variable impurity level that changes with x.
       */
      virtual void SetImpurity(TF1 * Im); // Used for the real change in impurity over x 
      /**
       * Returns the value for E under the first direction.
       */
      virtual double GetE1(double x,double y,double z){return GetData(x,2);};
      virtual double GetE2(double x,double y,double z){return 0;};
      virtual double GetE3(double x,double y,double z){return 0;};
      /**
       * Returns the impurity level.
       */ 
      virtual double GetImpurity(double x){return GetData(x,0);};
      /**
       * Returns the potential.
       */ 
      virtual double GetPotential(double x){return GetData(x,1);};
      /**
       * Returns the location under the first directions.
       */ 
      virtual double GetXEdge(bool beginorend); 
      /**
       *This defines the class for the cint dictionary.
       */
      ClassDef(X,1);

   protected:
      bool * fIsFixed; ///< Is used to check if a value is fixed or if it can be modified. It is usally used to check boundary conditions and find out if you are on the edge.
      bool fIsLoaded; ///< fIsLoaded is used to check if points in the grid have a value or not. If fIsLoaded returns true, the the points if the grid have value and you do not need to initialize, if it returns false you do.
      double *fE1; ///< Electric field under the first coordinate (x, r, or rho) direction 
      double  *fPotential; ///< Potential in the referenced grid
      double *fC1; ///< the location under the first coordinate (x, r, or rho) direction
      double *fDistanceToNext;///< Distance to the next point under the current coordinate direction
      double *fDistanceToPrevious;///< Distance to the previous point under the current coordinate direction
      double *fImpurity; ///< Value of the impurity level at a point on the grid

      /**
       * Sets the field step length.
       */
      virtual void SetStepLength(double steplength);
      /**
       * Uses a binary search to return the index .
       */
      virtual int FindIdx(double tarx,int begin,int end);
      /**
       * Calculates the field using an analytical method.
       */
      virtual bool Analytic();
      /**
       * Returns data for various variables. 
       * Case 0: Returns impurity
       * Case 1: Returns the potential
       * Case 2: Returns E1
       */
      virtual double GetData(double tarx,int thing); 
      /**
       * Uses Second order Successive Over-Relaxation method to calculate the field.
       */
      virtual void SOR2(int idx,bool calculateElectricField); 

};
#endif

