// definition of necessary units
static const double cm=1;
static const double cm3=cm*cm*cm;
static const double volt=1;
static const double C=1; // Coulomb
static const double e=1.6e-19*C; // elementary charge
static const double epsilon0=8.854187817e-14*C/volt/cm; // vacuum permittivity
// https://link.springer.com/chapter/10.1007/10832182_519
static const double epsilon=15.8; // Ge dielectric constant
//______________________________________________________________________________
// V"(x)=a, https://www.wolframalpha.com/input/?i=V%27%27(x)%3Da
double V(double *coordinates, double *parameters)
{
   double x = coordinates[0];// there is no y and z dependence
   double x0= parameters[0]; // lower electrode
   double x1= parameters[1]; // upper electrode
   double v0= parameters[2]; // lower voltage
   double v1= parameters[3]; // upper voltage
   double rho=parameters[4]; // space charge density [C/cm3]

   double a =-rho/epsilon0/epsilon;
   double c2= (v1-v0)/(x1-x0) - a/2*(x1+x0);
   double c1= (v0*x1-v1*x0)/(x1-x0) + a/2*x0*x1;
   return a*x*x/2 + c2*x + c1;
}
//______________________________________________________________________________
// E=-V'
double E(double *coordinates, double *parameters)
{
   double x = coordinates[0];
   double x0= parameters[0];
   double x1= parameters[1];
   double v0= parameters[2];
   double v1= parameters[3];
   double rho=parameters[4];

   double a =-rho/epsilon0/epsilon;
   double c2= (v1-v0)/(x1-x0) - a/2*(x1+x0);
   return -a*x - c2;
}
//______________________________________________________________________________
//
const int n=5; // number of curves
double rho[n]={-3.5e10*e/cm3, -1.5e10*e/cm3, 0, 1.5e10*e/cm3, 3.5e10*e/cm3};

// search for depletion voltage of a planar detector with given impurity

void voltage(double impurity=5e10/*1/cm^3*/, double thickness=1/*cm*/)
{
   // potential due to impurity alone
   GeFiCa::Planar1D impurityPotential;
   impurityPotential.SetImpurity(impurity/GeFiCa::cm3);
   impurityPotential.V1=0;
   impurityPotential.UpperBound=thickness*GeFiCa::cm;
   impurityPotential.CalculatePotential(GeFiCa::kSOR2);

   // weighting potential
   GeFiCa::Planar1D weightingPotential;
   weightingPotential.SetImpurity(0);
   weightingPotential.V1=1*GeFiCa::volt;
   weightingPotential.UpperBound=thickness*GeFiCa::cm;
   weightingPotential.CalculatePotential(GeFiCa::kSOR2);

   // total potential = impurity potential + V * weighting potential
   GeFiCa::Planar1D totalPotential;
   double lower=0, upper=20000*GeFiCa::volt; // limits of trials
   double vdep;
   while (lower<upper) {
      vdep=(upper+lower)/2;
      totalPotential.Copy(weightingPotential);
      totalPotential*=vdep;
      totalPotential+=&impurityPotential;
      if (totalPotential.IsDepleted()) upper=vdep-1e-5;
      else lower=vdep+1e-5;
   }
   cout<<"depletion voltage is: "<<vdep/GeFiCa::volt<<" V"<<endl;
}

//______________________________________________________________________________
//
bool isDepleted(double we0,double we1,double im0, double im1)
{
   return ((we0+im0<0)==(we1+im1<0))||(we0+im0==0)||(we1+im1==0);
}
//______________________________________________________________________________
//
double *findDepletedVoltage()
{
   TF1 *fwe=new TF1("fwe", E, 0, 1*cm,5);
   fwe->SetParameters(0,1*cm,0,1*cm,0);

   TF1 *fim[n]={0};
   double x0[n]={0}, x1[n], v0[n]={0}, v1[n];
   double *dv=new double[n];
   for (int i=0; i<n; i++) 
   {
      x1[i] = 1*cm;
      v1[i] = 0*volt;
      fim[i] = new TF1(Form("fim%d",i), E, x0[i], x1[i], 5);
      fim[i]->SetParameters(x0[i],x1[i],v0[i],v1[i],rho[i]);
      double upperdepletedvoltage=1e10; 
      double lowerdepletedvoltage=0; 
      while(upperdepletedvoltage>=lowerdepletedvoltage)
      {
         double mid=(upperdepletedvoltage+lowerdepletedvoltage)/2;
         if(isDepleted(mid* fwe->Eval(x0[i]),mid*fwe->Eval(x1[i]),fim[i]->Eval(x0[i]),fim[i]->Eval(x1[i])))
         {
            upperdepletedvoltage=mid-1e-5;
         }
         else lowerdepletedvoltage=mid+1e-5;
      }
      dv[i]=upperdepletedvoltage;
   }
   return dv; 
}
Vdep()
{
   double *result=findDepletedVoltage();
   for(int i=0;i<5;i++)
   {
      cout<<rho[i]<<" "<<result[i]<<endl;
   }
}
