#ifndef GeFiCa_X_H
#define GeFiCa_X_H

#include <TObject.h>
class TF1;

namespace GeFiCa { 
   enum EMethod {
      kAnalytic,
      kSOR2,
      kSOR4,
   };

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


      virtual bool CalculateField(EMethod method=kSOR2);
      
      virtual void SaveField(const char *fout);
      virtual void LoadField(const char *fin);
      virtual void SetImpurity(double density);
      virtual void SetImpurity(TF1 * Im);

      virtual double GetE1(double x){return GetData(x,2);}; 
      virtual double GetImpurity(double x){return GetData(x,0);}; 
      virtual double GetPotential(double x){return GetData(x,1);}; 
      virtual double GetXEdge(bool beginorend); 

      ClassDef(X,1);

   protected:
      bool * fIsFixed; // will this grid calculate
      bool fIsLoaded; // is this grid calculated before
      double *fE1; // Electric field under first direction 
      double  *fPotential; // Potential in this grid
      double *fC1; // the location under first direction
      double *fDistanceToNext,*fDistanceToPrevious,*fImpurity;

      virtual void SetStepLength(double steplength);
      virtual int FindIdx(double tarx,int begin,int end);

      virtual bool Analytic();
      
      virtual double GetData(double tarx,int thing); 
      virtual void SOR2(int idx,bool calculateElectricField); 
      
};

#endif
