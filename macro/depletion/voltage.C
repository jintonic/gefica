// search for depletion voltage of a planar detector with given impurity
{
   // potential due to impurity alone
   GeFiCa::Planar1D impurityPotential = GeFiCa::Planar1D();
   impurityPotential.SetImpurity(2e10/GeFiCa::cm3);
   impurityPotential.V0=0;
   impurityPotential.V1=0;
   impurityPotential.Csor=1.994;
   impurityPotential.CalculatePotential(GeFiCa::kSOR2);

   // weighting potential
   GeFiCa::Planar1D weightingPotential = GeFiCa::Planar1D();
   weightingPotential.SetImpurity(0);
   weightingPotential.V0=0;
   weightingPotential.V1=1*GeFiCa::volt;
   weightingPotential.Csor=1.994;
   weightingPotential.CalculatePotential(GeFiCa::kSOR2);

   // total potential = impurity potential + V * weighting potential
   GeFiCa::Planar1D totalPotential = GeFiCa::Planar1D();
   double upper=10000*GeFiCa::volt; // upper limit of trials
   double lower=0, mid;
   while (lower<=upper) {
      mid=(upper+lower)/2;
      totalPotential.CopyField(&weightingPotential);
      totalPotential*=mid;
      totalPotential+=&impurityPotential;
      if (totalPotential.IsDepleted()) upper=mid-1e-5;
      else lower=mid+1e-5;
   }
   cout<<"depletion voltage is: "<<mid<<endl;
   totalPotential.SaveField("depleted.root");
}
