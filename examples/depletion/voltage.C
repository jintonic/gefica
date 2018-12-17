namespace GeFiCa;
voltage(){
   GeFiCa::PointContactRZ *wholedetector = new GeFiCa::PointContactRZ(690,505);
   wholedetector->Radius=3.45*GeFiCa::cm;
   wholedetector->ZUpperBound=5.05*GeFiCa::cm;
   wholedetector->PointR=0.14*GeFiCa::cm;
   wholedetector->PointDepth=0*GeFiCa::cm;
   detector2->ContactInnerR=1.450;

   wholedetector->MaxIterations=1e5;
   wholedetector->Precision=1e-7*GeFiCa::volt;
   wholedetector->Csor=1.994;
   wholedetector->V0=2500*GeFiCa::volt;
   wholedetector->V1=0*GeFiCa::volt;

   TF2 *im=new TF2("f","-0.318e10+0.025e10*y");
   wholedetector->SetImpurity(im);
   //wholedetector->Impurity="-0.318e10+0.025e10*y";//-0.01e10/GeFiCa::cm3);


   GeFiCa::PointContactRZ *justimpurity = new GeFiCa::PointContactRZ(690,505);
   justimpurity->Radius=3.45;
   justimpurity->ZUpperBound=5.05;
   justimpurity->PointR=0.14;
   justimpurity->PointDepth=0.;
   justimpurity->ContactInnerR=1.450;

   justimpurity->MaxIterations=1e5;
   justimpurity->Precision=1e-7;
   justimpurity->Csor=1.994;

   TF2 *im2=new TF2("f","-0.318e10+0.025e10*y");
   justimpurity->SetImpurity(im2);
   //justimpurity->Impurity="-0.318e10+0.025e10*y";//-0.01e10/GeFiCa::cm3);


   justimpurity->V1=0;
   justimpurity->V0=0;
   justimpurity->CalculatePotential(GeFiCa::kSOR2);
   justimpurity->SaveField("im");
   justimpurity->LoadField("im");

   GeFiCa::PointContactRZ *weightingPotential = new GeFiCa::PointContactRZ(690,505);
   weightingPotential->Radius=3.45*GeFiCa::cm;
   weightingPotential->ZUpperBound=5.05*GeFiCa::cm;
   weightingPotential->PointR=0.14*GeFiCa::cm;
   weightingPotential->PointDepth=0.*GeFiCa::cm;
   weightingPotential->ContactInnerR=1.450;

   weightingPotential->MaxIterations=1e5;
   weightingPotential->Precision=1e-7*GeFiCa::volt;
   weightingPotential->Csor=1.994;
   weightingPotential->V0=1*GeFiCa::volt;
   weightingPotential->V1=0*GeFiCa::volt;


   TF2 *im3=new TF2("f","0.0+0.*y");
   weightingPotential->SetImpurity(im3);
   weightingPotential->CalculatePotential(GeFiCa::kSOR2);
   weightingPotential->SaveField("wp");
   weightingPotential->LoadField("wp");
   weightingPotential->CalculateCapacitance();
   //wholedetector->Copy(*weightingPotential);
   wholedetector->LoadField("wp");
   (*wholedetector) *= (double)2;
   GeFiCa::PointContactRZ * multiweightPotential=new GeFiCa::PointContactRZ(100,100);
    //multiweightPotential  ->Copy(weightingPotential);
    multiweightPotential  ->LoadField("wp");
   (*wholedetector) += justimpurity;
   wholedetector->SaveField("beforejump");
   int stepsize=2;
   while(!wholedetector->IsDepleted())
   {
      (*multiweightPotential)*=(double)2;
      (*wholedetector)+=multiweightPotential;
      stepsize=stepsize*2;
   }
   double upper=stepsize;
   std::cout<<stepsize;
   double lower=0;
   double mid=upper/2;
   while(lower<=upper)
   {
      mid=(upper+lower)/2;
      //wholedetector->Copy(weightingPotential);
      wholedetector->LoadField("wp");
      (*wholedetector)*=mid;
      (*wholedetector)+=justimpurity;
      if(wholedetector->IsDepleted())
      {
         upper=mid-1e-5;
      }
      else
      {
         lower=mid+1e-5;
      }
   }
   wholedetector->SaveField("result.root");
}

