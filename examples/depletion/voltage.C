// binary search for the depletion voltage of a point contact detector
using namespace GeFiCa;
void voltage()
{
   PointContactDZ *vi; // potential (v) due to impurity (i) alone
   PointContactDZ *vw; // weighting (w) potential (v)
   PointContactDZ *vt; // totoal (t) potential (v) = vw * bias + vi
   TFile *fi; // ROOT file used to save vi
   TFile *fw; // ROOT file used to save vw

   cout<<"\ncalculate or load potential due to space charge alone:"<<endl;
   if (FILE *input = fopen("impurity.root","r")) { // load from fi
      fclose(input);
      fi = new TFile("impurity.root"); // read only
      vi = (PointContactDZ*) fi->Get("vi");
   } else { // a fresh calculation
      fi = new TFile("impurity.root", "recreate");

      vi = new PointContactDZ(690,505); vi->SetName("vi");

      vi->V0=0; vi->V1=0; // no bias

      vi->Z=5.05*cm; vi->Radius=3.45*cm;
      vi->PointContactZ=0.21*cm; vi->PointContactR=0.14*cm;

      vi->Csor=1.994;
      vi->MaxIterations=1e5;
      vi->Precision=1e-7*volt;

      // x in TF3 -> r in PointContactDZ, y in TF3 -> z in PointContactDZ
      TF3 *impurityDistribution = new TF3("fi","-0.318e10+0.025e10*y");
      vi->SetImpurity(impurityDistribution);

      vi->CalculatePotential(kSOR2);

      vi->Write(); // save itself to fi
   }
   vi->Dump(); // print configurations

   cout<<"\ncalculate or load weighting potential:"<<endl;
   if (FILE *input = fopen("weighting.root","r")) { // load from fw
      fclose(input);
      fw = new TFile("weighting.root"); // read only
      vw = (PointContactDZ*) fw->Get("vw");
   } else { // a fresh calculation
      fw = new TFile("weighting.root", "recreate");

      vw = new PointContactDZ(690,505); vw->SetName("vw");

      vw->V0=1*volt; vw->V1=0*volt; // weighting potential

      vw->Z=5.05*cm; vw->Radius=3.45*cm;
      vw->PointContactZ=0.21*cm; vw->PointContactR=0.14*cm;

      vw->Csor=1.994;
      vw->MaxIterations=1e5;
      vw->Precision=1e-7*volt;

      vw->SetImpurity(new TF3("wpi","0")); // no impurity

      vw->CalculatePotential(kSOR2);

      vw->Write(); // save itself to fw
   }
   vw->Dump(); // print configurations

   double bias, vLower=0*volt, vUpper=2e4*volt; // range of search
   cout<<"\nStart binary search in ["<<vLower<<", "<<vUpper<<"] V"<<endl;
   while (vUpper-vLower>0.1*volt) { // binary search
      bias=(vUpper+vLower)/2; // a new guess
      // vt = vw*bias + vi
      vt = (PointContactDZ*) vw->Clone("vt");
      (*vt)*=bias; (*vt)+=vi;
      // updating range
      if (vt->IsDepleted()) vUpper=bias;
      else vLower=bias;

      delete vt;
      cout<<"Current guess: "<<bias<<" V, ";
      cout<<"new search range: ["<<vLower<<", "<<vUpper<<"]"<<endl;
   }
   cout<<"The depletion voltage is found to be: "<<bias<<" V"<<endl;

   fw->Close();
   fi->Close();
   delete fi;
   delete fw;
}

