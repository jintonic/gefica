// binary search for the depletion voltage of a point contact detector
using namespace GeFiCa;
void voltage()
{
   TFile *fi, *fw;
   PointContactDZ *vi, *vw, *vt;

   // potential due to space charge alone
   if (FILE *input = fopen("impurity.root","r")) { // if already calculated
      fclose(input);
      fi = new TFile("impurity.root"); // read only
      vi = (PointContactDZ*) fi->Get("vi");
   } else { // if not calculated before
      fi = new TFile("impurity.root", "recreate");

      vi = new PointContactDZ(690,505);
      vi->SetName("vi");

      vi->V0=0; // no bias
      vi->V1=0; // no bias

      vi->Z=5.05*cm;
      vi->Radius=3.45*cm;
      vi->PointContactZ=0.21*cm;
      vi->PointContactR=0.14*cm;

      vi->Csor=1.994;
      vi->MaxIterations=1e5;
      vi->Precision=1e-7*volt;

      TF3 *impurityDistribution = new TF3("fi","-0.318e10+0.025e10*y");
      vi->SetImpurity(impurityDistribution);

      vi->CalculatePotential(kSOR2);

      vi->Write();
   }
   vi->Dump();

   // weighting potential
   if (FILE *input = fopen("weighting.root","r")) { // if already calculated
      fclose(input);
      fw = new TFile("weighting.root"); // read only
      vw = (PointContactDZ*) fw->Get("vw");
   } else { // if not calculated before
      fw = new TFile("weighting.root", "recreate");

      vw = new PointContactDZ(690,505);
      vw->SetName("vw");

      vw->V0=1*volt;
      vw->V1=0*volt;

      vw->Z=5.05*cm;
      vw->Radius=3.45*cm;
      vw->PointContactZ=0.21*cm;
      vw->PointContactR=0.14*cm;

      vw->Csor=1.994;
      vw->MaxIterations=1e5;
      vw->Precision=1e-7*volt;

      vw->SetImpurity(new TF3("wpi","0")); // no impurity

      vw->CalculatePotential(kSOR2);

      vw->Write();
   }
   vw->Dump();

   double bias, vlower=0*volt, vupper=2e4*volt; // range of search
   while (vupper-vlower>1e-3*volt) { // binary search
      bias=(vupper+vlower)/2; // bias voltage
      vt = (PointContactDZ*) vw->Clone("vt");
      (*vt)*=bias;
      (*vt)+=vi;
      if (vt->IsDepleted()) vlower=bias;
      else vupper=bias;
      delete vt;
      cout<<"bias: "<<bias<<", u: "<<vupper<<", l: "<<vlower<<endl;
   }
   cout<<"depletion voltage: "<<bias<<endl;

   fw->Close();
   fi->Close();
   delete fi;
   delete fw;
}

