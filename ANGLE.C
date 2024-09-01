#include <cmath>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include <cassert>
using namespace ROOT::Math;

void line(double t, const double *p, double &x, double &y, double &z) {
   x = 150  + t;
   y = p[2] + p[3]*t;
   z = p[0] + p[1]*t;
}
struct SumDistance2 {
   TGraph2DErrors *fGraph;
   SumDistance2(TGraph2DErrors *g) : fGraph(g) {}
   double distance2(double x, double y, double z, double ex, double ey, double ez, const double *p) {
      XYZVector xp(x,y,z);
      XYZVector x0(150, p[2], p[0]);
      XYZVector x1(151, p[2] + p[3], p[0] + p[1]);
      XYZVector u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      //double weight = 1 / (ex * ex + ey * ey + ez * ez);
      double weight = 1 / (ex * ex);
      return d2 * weight;
   }
   double operator() (const double *par) {
      assert(fGraph != nullptr);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      
      double *ex = fGraph->GetEX();
      double *ey = fGraph->GetEY();
      double *ez = fGraph->GetEZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i  = 0; i < npoints; ++i) {
         double d = distance2(x[i],y[i],z[i], ex[i], ey[i], ez[i], par);
         sum += d;
      }
      return sum;
   }
};
void cylindrical_to_cartesian(int r, double z, double phi, double &x, double &y, double &z_out) {
    phi += 8;
    double PH0 = 0.9;
    double PH_I = 66.6;
    if (r == 0) {
        double RRR = 98.5;
        double PH = 1.8;
        double ZZ = 6.95;
        double DT0 = 11;
        double angle = (PH * phi + PH_I + PH0) * M_PI / 180;
        x = RRR * std::sin(angle);
        y = RRR * std::cos(angle);
        z_out = DT0 * z + ZZ;
    } else if (r == 1) {
        double RRR = 98.5;
        double PH = 1.8;
        double ZZ = 9.15;
        double DT0 = 11;
        double angle = (PH * phi + PH_I + PH0) * M_PI / 180;
        x = RRR * std::sin(angle);
        y = RRR * std::cos(angle);
        z_out = DT0 * z + ZZ;
    } else if (r >= 2 && r <= 4) {
        double RRR = 95.5 + r * 3;
        double PH = 1.8;
        double ZZ = 5.5;
        double DT = 11;
        double angle = (PH * phi + PH_I + PH0) * M_PI / 180;
        x = RRR * std::sin(angle);
        y = RRR * std::cos(angle);
        z_out = DT * z + ZZ;
    }
    }
void ANGLE(){
    TFile *fc = TFile::Open("output.root");
    if (!fc || fc->IsZombie()) {
        std::cerr << "Error: cannot open input file 'outputG.root'" << std::endl;
        return;
    }
    TTree *fChainc = (TTree*)fc->Get("tree");
    if (!fChainc) {
        std::cerr << "Error: cannot find tree 't' in the file" << std::endl;
        return;
    }
    double PrimaryDirX, PrimaryDirY, PrimaryDirZ, DirX, DirY, DirZ;
    std::vector<int> *r= nullptr;
    std::vector<int> *phi= nullptr;
    std::vector<int> *zp= nullptr;
    std::vector<float> *LG= nullptr;
    fChainc->SetBranchAddress("r", &r);
    fChainc->SetBranchAddress("phi", &phi);
    fChainc->SetBranchAddress("zp", &zp);
    fChainc->SetBranchAddress("LG", &LG);
    TFile* output = TFile::Open("ANGLEoutput.root", "RECREATE");
    TTree* tree = new TTree("tree", "Fitted Data");
    double param1, param2, param3, param4, param5, err1, errdirZ, err3, errdirY, chi2, NDF, MinFCN, energyTot, dirX, dirY, dirZ, azimuth, elevation, errazimuth, errelevation;
    tree->Branch("dirX", &dirX);
    tree->Branch("dirY", &dirY);
    tree->Branch("dirZ", &dirZ);
    tree->Branch("errdirY", &errdirY);
    tree->Branch("errdirZ", &errdirZ);
    tree->Branch("azimuth", &azimuth);
    tree->Branch("elevation", &elevation);
    tree->Branch("errazimuth", &errazimuth);
    tree->Branch("errelevation", &errelevation);
    tree->Branch("energyTot", &energyTot);
    tree->Branch("chi2", &chi2);
    tree->Branch("NDF", &NDF);
    tree->Branch("MinFCN", &MinFCN);
    int nbytes = 0, nb = 0;
    int nentries=fChainc->GetEntries();  
    double sum = 0.0;
    for (Long64_t jentry=0; jentry<nentries; jentry++) { 
       Long64_t ientry = fChainc->LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChainc->GetEntry(jentry);   
       nbytes += nb;
       cout<<"== EVENT "<<jentry<<" =="<<endl;
       std::vector<float> x, y, z, energyF;
       double A, B, C;
       energyTot=0; 
       gStyle->SetOptStat(0);
       gStyle->SetOptFit();
       TGraph2DErrors * gr = new TGraph2DErrors();
       for(int i = 0; i < r->size(); ++i){
          cylindrical_to_cartesian(r->at(i), zp->at(i), phi->at(i), A, B, C);
          x.push_back(A);
          y.push_back(B);
          z.push_back(C);
          energyF.push_back(LG->at(i));
          energyTot=energyTot+LG->at(i);
          gr->SetPoint(i,A,B,C);
          gr->SetPointError(i,1/LG->at(i),1/LG->at(i),1/LG->at(i));
       }
       ROOT::Fit::Fitter  fitter;
       SumDistance2 sdist(gr);
       ROOT::Math::Functor fcn(sdist,4);
       double P_i_0 = 0.0, P_i_1 = 0.0, P_i_2 = 0.0, P_i_3 = 0.0;
       double sum_x = 0.0, sum_y = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
       for (int i = 0; i < r->size(); ++i) {
         sum_x += x[i];
         sum_y += y[i];
       }
       P_i_0 = sum_x / x.size();
       P_i_2 = sum_y / y.size();
       for (int i = 0; i < r->size(); ++i) {
         sum_x2 += TMath::Power(x[i] - P_i_0, 2);
         sum_y2 += TMath::Power(y[i] - P_i_2, 2);
       }
       P_i_1 = TMath::Sqrt(sum_x2 / x.size());
       P_i_3 = TMath::Sqrt(sum_y2 / y.size());
       
       double pStart[4] = {P_i_0, P_i_1, P_i_2, P_i_3};
       
       fitter.SetFCN(fcn,pStart);
       for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.000001);
       bool ok = fitter.FitFCN();
       if (!ok) {
         Error("line3Dfit","Line3D Fit failed");
         continue;
       }
       const ROOT::Fit::FitResult & result = fitter.Result();
       double n=sqrt(result.Parameter(1)*result.Parameter(1)+result.Parameter(3)*result.Parameter(3)+1);
       param1= result.Parameter(0);
       dirZ= result.Parameter(1)/n;
       param3= result.Parameter(2);
       dirY= result.Parameter(3)/n;
       dirX= 1/n;
       err1= result.ParError(0);
       errdirZ= result.ParError(1)/n;
       err3= result.ParError(2);
       errdirY= result.ParError(3)/n;
       azimuth = (TMath::ATan(dirY/dirZ))* 180 / TMath::Pi();
       elevation = abs(TMath::ATan(dirX/sqrt(dirY*dirY+dirZ*dirZ)) * 180 / TMath::Pi());
       
       double dy_1=1/(1+pow(dirY/dirZ, 2))/dirZ, dz_1=1/(1+pow(dirY/dirZ, 2))*dirY/pow(dirZ, 2);
       double dy_2=1/(1+pow(dirX/sqrt(dirY*dirY+dirZ*dirZ), 2))*pow(dirX/(dirY*dirY+dirZ*dirZ), 2)*2*dirY, dz_2=1/(1+pow(dirX/sqrt(dirY*dirY+dirZ*dirZ), 2))*pow(dirX/(dirY*dirY+dirZ*dirZ), 2)*2*dirZ;
       errazimuth=180/TMath::Pi()*sqrt(pow(dy_1*errdirY,2)+pow(dz_1*errdirZ,2)+dy_1*dz_1*errdirY*errdirZ); 
       errelevation=180/TMath::Pi()*sqrt(pow(dy_2*errdirY,2)+pow(dz_2*errdirZ,2)+dy_1*dz_1*errdirY*errdirZ); 
       NDF= x.size();
       MinFCN=result.MinFcnValue();
       const double * parFit = result.GetParams();
       chi2=0;
       double dx, dy, dz, xi, yi, zi;
       for (int i = 0; i < r->size(); ++i) {
        line(z[i], parFit, xi, yi, zi);
        dx = (x[i] - xi) ;
        dy = (y[i] - yi) ;
        dz = (z[i] - zi) ;
        chi2 += dx*dx + dy*dy + dz*dz;
       }
       tree->Fill();
  }
  tree->Write();
  output->Close();
  cout<<"out of events loop"<<endl;

}
