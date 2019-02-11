// binary search for the depletion voltage of a point contact detector
using namespace GeFiCa;
void voltage()
{
   // potentials from both impurity and bias
   PointContactDZ *totalPotential = new PointContactDZ(690,505);

   totalPotential->V0=2500*volt;
   totalPotential->V1=0*volt;

   totalPotential->Z=5.05*cm;
   totalPotential->Radius=3.45*cm;
   totalPotential->PointContactZ=0.21*cm;
   totalPotential->PointContactR=0.14*cm;

   totalPotential->Csor=1.994;
   totalPotential->MaxIterations=1e5;
   totalPotential->Precision=1e-7*volt;

   TF3 *impurityDistribution = new TF3("fi","-0.318e10+0.025e10*y");
   totalPotential->SetImpurity(impurityDistribution);

   // potential due to space charge alone
   PointContactDZ *impurityPotential = new PointContactDZ(690,505);

   impurityPotential->V0=0; // no bias
   impurityPotential->V1=0; // no bias

   impurityPotential->Z=5.05*cm;
   impurityPotential->Radius=3.45*cm;
   impurityPotential->PointContactZ=0.21*cm;
   impurityPotential->PointContactR=0.14*cm;

   impurityPotential->Csor=1.994;
   impurityPotential->MaxIterations=1e5;
   impurityPotential->Precision=1e-7*volt;

   impurityPotential->SetImpurity(impurityDistribution);

   impurityPotential->CalculatePotential(kSOR2);
   //impurityPotential->SaveField("im");
   //impurityPotential->LoadField("im");

   // weighting potential
   PointContactDZ *weightingPotential = new PointContactDZ(690,505);

   weightingPotential->V0=1*volt;
   weightingPotential->V1=0*volt;

   weightingPotential->Z=5.05*cm;
   weightingPotential->Radius=3.45*cm;
   weightingPotential->PointContactZ=0.21*cm;
   weightingPotential->PointContactR=0.14*cm;

   weightingPotential->Csor=1.994;
   weightingPotential->MaxIterations=1e5;
   weightingPotential->Precision=1e-7*volt;

   weightingPotential->SetImpurity(new TF3("wpi","0"));

   weightingPotential->CalculatePotential(kSOR2);
   //weightingPotential->SaveField("wp");
   //weightingPotential->LoadField("wp");

   totalPotential = (PointContactDZ*) weightingPotential->Clone("total");
   (*totalPotential) *= 2.;

   PointContactDZ *multiweightPotential = new PointContactDZ(100,100);
   multiweightPotential = (PointContactDZ*) weightingPotential->Clone("total");

   (*totalPotential) += impurityPotential;
   totalPotential->SaveField("beforejump");
   int stepsize=2;
   while(!totalPotential->IsDepleted()) {
      (*multiweightPotential)*=(double)2;
      (*totalPotential)+=multiweightPotential;
      stepsize=stepsize*2;
   }
   double upper=stepsize;
   std::cout<<stepsize;
   double lower=0;
   double mid=upper/2;
   while(lower<=upper) {
      mid=(upper+lower)/2;
      totalPotential = (PointContactDZ*) weightingPotential->Clone("total");
      //totalPotential->LoadField("wp");
      (*totalPotential)*=mid;
      (*totalPotential)+=impurityPotential;
      if(totalPotential->IsDepleted()) {
         upper=mid-1e-5;
      } else {
         lower=mid+1e-5;
      }
   }
   totalPotential->SaveField("result.root");
}

