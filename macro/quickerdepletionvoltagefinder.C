GeFiCa::PointContactRZ *im()
{
}

quickerdepletionvoltagefinder()
{
   GeFiCa::PointContactRZ *wholedetector = new GeFiCa::PointContactRZ(690,505);
   wholedetector->Radius=3.45*GeFiCa::cm;
   wholedetector->ZUpperBound=5.05*GeFiCa::cm;
   wholedetector->PointR=0.14*GeFiCa::cm;
   wholedetector->PointDepth=0.21*GeFiCa::cm;

   wholedetector->MaxIterations=1e5;
   wholedetector->Precision=1e-7*GeFiCa::volt;
   wholedetector->Csor=1.994;
   wholedetector->V0=2500*GeFiCa::volt;
   wholedetector->V1=0*GeFiCa::volt;

   wholedetector->Impurity="-0.318e10+0.025e10*y";//-0.01e10/GeFiCa::cm3);


   GeFiCa::PointContactRZ *justimpurity = new GeFiCa::PointContactRZ(690,505);
   justimpurity->Radius=3.45;
   justimpurity->ZUpperBound=5.05;
   justimpurity->PointR=0.14;
   justimpurity->PointDepth=0.21;

   justimpurity->MaxIterations=1e5;
   justimpurity->Precision=1e-7;
   justimpurity->Csor=1.994;

   justimpurity->Impurity="-0.318e10+0.025e10*y";//-0.01e10/GeFiCa::cm3);


   justimpurity->V1=0;
   justimpurity->V0=0;
   //justimpurity->CalculatePotential(GeFiCa::kSOR2);
   //justimpurity->SaveField("im");
   justimpurity->LoadField("im");

   GeFiCa::PointContactRZ *weightingPotential = new GeFiCa::PointContactRZ(690,505);
   weightingPotential->Radius=3.45*GeFiCa::cm;
   weightingPotential->ZUpperBound=5.05*GeFiCa::cm;
   weightingPotential->PointR=0.14*GeFiCa::cm;
   weightingPotential->PointDepth=0.21*GeFiCa::cm;

   weightingPotential->MaxIterations=1e5;
   weightingPotential->Precision=1e-7*GeFiCa::volt;
   weightingPotential->Csor=1.994;
   weightingPotential->V0=1*GeFiCa::volt;
   weightingPotential->V1=0*GeFiCa::volt;


   weightingPotential->Impurity="0+0*y";
   //weightingPotential->CalculatePotential(GeFiCa::kSOR2);
   //weightingPotential->SaveField("wp");
   weightingPotential->LoadField("wp.root");
   weightingPotential->CalculateCapacitance();
   int stepsize=1;
   wholedetector->CopyField(weightingPotential);
   wholedetector=wholedetector*2;
   GeFiCa::PointContactRZ * multiweightPotential=new GeFiCa::PointContactRZ(100,100);
    multiweightPotential  ->CopyField(weightingPotential);
   wholedetector=wholedetector+justimpurity;
   wholedetector->SaveField("beforejump");
   int stepsize=2;
   while(!wholedetector->Depleattest())
   {
      multiweightPotential=multiweightPotential*2;
      wholedetector=wholedetector+multiweightPotential;
      stepsize=stepsize*2;
   }
   double upper=stepsize;
   double lower=0;
   double mid=upper/2;
   while(lower<=upper)
   {
      mid=(upper+lower)/2;
      wholedetector->CopyField(weightingPotential);
      wholedetector=wholedetector*mid;
      wholedetector=wholedetector+justimpurity;
      if(wholedetector->Depleattest())
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

