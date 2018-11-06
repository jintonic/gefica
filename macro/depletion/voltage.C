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
