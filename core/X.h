#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TObject.h>
class TF1;

namespace GeFiCa { 

	/**
	* EMethod is a namespace with three different methods of field calculation.
	*/
   enum EMethod {
      kAnalytic,	/**< Analytical method*/
      kSOR2,		/**< Successove over-relaxation method to the 2nd order*/
      kSOR4,		/**< Successove over-relaxation method to the 4th order*/
   };

   /**
   * X is a member of GeFiCa.
   */
   class X;

   static const double epsilon = 16*8.854187817e-12; 
   static const double cm =1;
   static const double volt=1;
   static const double cm3=cm*cm*cm; 
}

/**
 * 1 D grid for field calculation.
 */
class GeFiCa::X : public TObject 
{
   public:
      int n1; ///< number of grid along the 1st axis
      int MaxIterations; ///< max one turn Iteration number
      int n; ///< n = n1 total number of grid
      double Csor; ///< boost Iteration speed
      double Precision; ///< X limit

   public:
      X(int nx=101);

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
	   * This function takes an argument for the variable density.
	   * It is usualy in the form of 1e10/cm3 or some variation.
	   */
      virtual void SetImpurity(double density);
      /**
      * Another Important method involved in setting the impurity.
      */
      virtual void SetImpurity(TF1 * Im);
      /**
      * Returns the value for E under the first direction.
      */
      virtual double GetE1(double x){return GetData(x,2);};
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

      ClassDef(X,1);

   protected:
      bool * fIsFixed; // will this grid calculate
      bool fIsLoaded; // is this grid calculated before
      double *fE1; // Electric field under first direction 
      double  *fPotential; // Potential in this grid
      double *fC1; // the location under first direction
      double *fDistanceToNext,*fDistanceToPrevious,*fImpurity;
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
