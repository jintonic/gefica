// binary search for the depletion voltage of a point contact detector
using namespace GeFiCa;
void voltage()
{
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

   TF3 *impurityDistribution = new TF3("fi","-0.318e10+0.025e10*y");
   impurityPotential->SetImpurity(impurityDistribution);

   if (FILE *input = fopen("impurity.root","r")) { // already calculated
      fclose(input);
      impurityPotential->LoadField("impurity.root");
   } else { // calculate the fields if not yet calculated
      impurityPotential->CalculatePotential(kSOR2);
      impurityPotential->SaveField("impurity.root");
   }

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

   weightingPotential->SetImpurity(new TF3("wpi","0")); // no impurity

   if (FILE *input = fopen("weighting.root","r")) { // already calculated
      fclose(input);
      weightingPotential->LoadField("weighting.root");
   } else { // calculate the fields if not yet calculated
      weightingPotential->CalculatePotential(kSOR2);
      weightingPotential->SaveField("weighting.root");
   }

   double bias, vlower=0*volt, vupper=2e4*volt; // range of search
   while (vupper-vlower>1e-3*volt) { // binary search
      bias=(vupper+vlower)/2; // bias voltage
      PointContactDZ totalPotential = *(PointContactDZ*) weightingPotential->Clone();
      totalPotential*=bias;
      totalPotential+=impurityPotential;
      if(totalPotential.IsDepleted()) vlower=bias;
      else vupper=bias;
   }
   cout<<"depletion voltage: "<<bias<<endl;
}

