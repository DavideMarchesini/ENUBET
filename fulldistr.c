double CustomEnergySpectrum(double E) {
    double k = 1.0;
    double alpha = 2;
    return k * TMath::Power(E, -alpha);
}

void fulldistr(int N=1000){
	TFile *file=new TFile("muons.root","RECREATE");
	TNtupleD* muons=new TNtupleD("fMuons","fMuons","fmuX:fmuY:fmuZ:fmuPX:fmuPY:fmuPZ:fmuE");
	double R=80.;
	double X0= 47.;
	double Y0= 0.;
	double Z0= 100.;
	double ThSky,PhSky,xbox,ybox,zbox,theta,pt,phi,px,py,pz, energy;
	TRandom3 *rand = new TRandom3();
	for (int i = 0; i < N; i++) {
          ThSky = TMath::ACos(TMath::Sqrt(rand->Rndm()));
          PhSky=2.*TMath::Pi()*gRandom->Rndm();
	  xbox=R*sin(ThSky)*cos(PhSky)+X0;			
	  ybox=R*sin(ThSky)*sin(PhSky)+Y0;
	  zbox=R*cos(ThSky)+Z0;
          double E = rand->Rndm(); 
          energy = CustomEnergySpectrum(E+0.1); 
          theta=TMath::ACos(cbrt(1.-2*gRandom->Rndm()));
          pt = TMath::Sin(theta);
          phi = 2*TMath::Pi()*gRandom->Rndm();
          px = pt * TMath::Cos(phi);
          py = pt * TMath::Sin(phi);
          pz = TMath::Sqrt(1 - px * px - py * py);
          if (pz < 0.) muons->Fill(xbox, zbox, ybox, px, pz, py, energy); 
          if (pz > 0.) muons->Fill(xbox, zbox, ybox, -px, -pz, -py, energy); 
        }
	muons->Print();
	muons->Write();
	file->Close();
}
