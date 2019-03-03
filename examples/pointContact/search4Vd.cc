// binary search for the depletion voltage of a point contact detector
using namespace GeFiCa;
void search4Vd()
{
   PointContactDZ *vi; // potential (v) due to impurity (i) alone
   PointContactDZ *vu; // potential (v) due to unit (u) bias alone
   PointContactDZ *vt; // totoal (t) potential (v) = vu * bias + vi
   TFile *fi; // ROOT file used to save vi
   TFile *fu; // ROOT file used to save vu

   cout<<"\ncalculate or load potential due to space charge alone:"<<endl;
   if (FILE *input = fopen("impurity.root","r")) { // load from fi
      fclose(input);
      fi = new TFile("impurity.root"); // read only
      vi = (PointContactDZ*) fi->Get("vi");
   } else { // a fresh calculation
      fi = new TFile("impurity.root", "recreate");

      vi = new PointContactDZ(690,505); vi->SetName("vi");
      vi->V0=0; vi->V1=0; // no bias
      vi->Height=5.05*cm; vi->Radius=3.45*cm;
      vi->PointContactH=0.21*cm; vi->PointContactR=0.14*cm;

      // x in TF3 -> r in PointContactDZ, y in TF3 -> z in PointContactDZ
      TF3 *fid = new TF3("fImpDistr","-0.318e10+0.025e10*y");
      vi->SetImpurity(fid);

      vi->Csor=1.994; // speed up SOR
      vi->CalculatePotential(kSOR2);

      vi->Write(); // save itself to fi
   }
   vi->Dump(); // print configurations

   cout<<"\ncalculate or load potential due to unit bias:"<<endl;
   if (FILE *input = fopen("oneVolt.root","r")) { // load from fu
      fclose(input);
      fu = new TFile("oneVolt.root"); // read only
      vu = (PointContactDZ*) fu->Get("vu");
   } else { // a fresh calculation
      fu = new TFile("oneVolt.root", "recreate");

      vu = new PointContactDZ(690,505); vu->SetName("vu");
      vu->V0=1*volt; vu->V1=0*volt;
      vu->Height=5.05*cm; vu->Radius=3.45*cm;
      vu->PointContactH=0.21*cm; vu->PointContactR=0.14*cm;
      vu->SetAverageImpurity(0); // no impurity

      vu->Csor=1.994;
      vu->CalculatePotential(kSOR2);

      vu->Write(); // save itself to fu
   }
   vu->Dump(); // print configurations

   double bias, vLower=0*volt, vUpper=2e4*volt; // range of search
   cout<<"\nStart binary search in ["<<vLower<<", "<<vUpper<<"] V"<<endl;
   while (vUpper-vLower>0.1*volt) { // binary search
      bias=(vUpper+vLower)/2; // a new guess
      // vt = vu*bias + vi
      vt = (PointContactDZ*) vu->Clone("vt");
      (*vt)*=bias; (*vt)+=vi;
      // updating range
      if (vt->IsDepleted()) vUpper=bias;
      else vLower=bias;

      delete vt;
      cout<<"Current guess: "<<bias<<" V, ";
      cout<<"new search range: ["<<vLower<<", "<<vUpper<<"]"<<endl;
   }
   cout<<"The depletion voltage is found to be: "<<bias<<" V"<<endl;

   fu->Close();
   fi->Close();
   delete fi;
   delete fu;
}
