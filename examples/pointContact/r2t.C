{
   GeFiCa::PointContactRZ *detector2 = new GeFiCa::PointContactRZ(1036,506);
   detector2->LoadField("point2dSOR2.root");

   ifstream infile("ev.new");
   ofstream outfile("result.txt");

   double x,y,v,er,ez,e,anotherV,E1,E2,E,sizeofr,sizeofz;
   while (infile>>x>>y>>v>>e>>er>>ez) {
      sizeofr=x;
      sizeofz=y;
   }
   outfile<<"# radius of detector: "<<sizeofr<<"\n";
   outfile<<"# height of detector: "<<sizeofz<<"\n";
   outfile<<"# number of grid on r: "<<n1<<"\n";
   outfile<<"# number of grid on z: "<<n2<<"\n";
   outfile<<"# voltage of one end: "<<detector2->V0<<"\n";
   outfile<<"# voltage of other end: "<<detector2->V1<<"\n";
   outfile<<"# r z V(gefica) dV Er(V/cm) Ez(V/cm) E(V/cm) dEr dEz dE "<<sizeofz<<"\n";
   outfile.seekg(0, std::ios::beg); 
   while (infile>>x>>y>>v>>e>>er>>ez) {
      sizeofr=x;
      sizeofz=y;
      anotherV=detector2->GetPotential(x/10,y/10);
      E1=detector2->GetE1(x/10,y/10);
      E2=detector2->GetE2(x/10,y/10);
      E=sqrt(E1*E1+E2*E2);
      outfile<<x<<"  "<<y<<"  "<<anotherV<<"  "<<v-anotherV<<"  "<<E1<<"  "<<E2<<"  "<<E<<"  "<<er-E1<<"  "<<ez-E2<<"  "<<e-E<<endl;
   }
   infile.close();
   outfile.close();
}
