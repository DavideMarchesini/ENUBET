#include "MappingScript.h"
void fG(){
    const int NNR=5;
    const int NNZ=15;
    const int NNPH=25;
    string filename="templateAnalysis.root";
    string fileInputCal="mappingMarchesiniActive.txt";
    string fileInputMap="mappingMarchesini.txt";
    std::map<std::vector<int>, std::vector<int>> dMap1=MappingAnodeCord(fileInputMap.c_str());
    ifstream in (fileInputCal.c_str());
    if (!in.is_open()) {
        std::cerr << "Error opening file " << fileInputCal << std::endl;
        return;
    }
    int numberOfAnodes = 64;
    int numberOfBoards = 13;
    int mRR, mZZ, mPP;
    int b1, a1, coeff;
    double corr[NNR][NNZ][NNPH]={0};
   cout<<"DONE"<<endl;
   int k=0;
    while (k<746) {
        in >> b1 >> a1 >> coeff;
        mZZ=dMap1[{b1,a1}][2];
        mRR=dMap1[{b1,a1}][0];
        mPP=dMap1[{b1,a1}][1];
        if (in.eof() == true) break;
        corr[mRR][mZZ][mPP+8]=coeff;
    }
    in.close();
    cout<<"DONE"<<endl;
  cout<<"================= RUN is ON ===================== "<<endl;
  TFile *fc = new TFile("templateAnalysis.root");
  fc->cd();
  TTree *fChainc = (TTree*)fc->Get("Demo");
  std::vector<int> *caloR= nullptr;
  std::vector<int> *caloZ= nullptr;
  std::vector<int> *caloPhi= nullptr;
  std::vector<int> *t0Z= nullptr;
  std::vector<int> *t0Phi= nullptr;
  std::vector<double> *t0UpEDep= nullptr;
  std::vector<double> *t0DwEDep= nullptr;
  std::vector<float> *ucmEDep= nullptr;
  double PrimaryDirX, PrimaryDirY, PrimaryDirZ, DirX, DirY, DirZ;
  fChainc->SetBranchAddress("caloR", &caloR);
  fChainc->SetBranchAddress("caloZ", &caloZ);
  fChainc->SetBranchAddress("caloPhi", &caloPhi);
  fChainc->SetBranchAddress("t0Z", &t0Z);
  fChainc->SetBranchAddress("t0Phi", &t0Phi);
  fChainc->SetBranchAddress("t0UpEDep", &t0UpEDep);
  fChainc->SetBranchAddress("t0DwEDep", &t0DwEDep);
  fChainc->SetBranchAddress("ucmEDep", &ucmEDep);
  fChainc->SetBranchAddress("PrimaryDirX", &PrimaryDirX);
  fChainc->SetBranchAddress("PrimaryDirY", &PrimaryDirY);
  fChainc->SetBranchAddress("PrimaryDirZ", &PrimaryDirZ);
  TFile *file = TFile::Open("outputG.root", "RECREATE");
  TTree *tree = new TTree("tree", "A tree with a vector");
  std::vector<int> *zp= nullptr;
  std::vector<int> *r= nullptr;
  std::vector<int> *phi= nullptr;
  std::vector<float> *energy= nullptr;
  tree->Branch("r", &r);
  tree->Branch("zp", &zp);
  tree->Branch("phi", &phi);
  tree->Branch("energy", &energy);
  tree->Branch("DirX", &DirX);
  tree->Branch("DirY", &DirY);
  tree->Branch("DirZ", &DirZ);

  int nbytes = 0, nb = 0;
  int nentries=fChainc->GetEntries();  
  for (Long64_t jentry=0; jentry<nentries; jentry++) { 
    Long64_t ientry = fChainc->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChainc->GetEntry(jentry);   nbytes += nb;
    cout<<"== EVENT "<<jentry<<" =="<<endl;
    DirX=PrimaryDirY;
    DirY=PrimaryDirX;
    DirZ=PrimaryDirZ;
    r = new std::vector<int>();
    zp = new std::vector<int>();
    phi = new std::vector<int>();
    energy = new std::vector<float>();
    
    for(int i = 0; i < ucmEDep->size(); ++i){
        if(corr[caloR->at(i)+2][caloZ->at(i)][caloPhi->at(i)+8]==1){
                r->push_back(caloR->at(i)+2);
	        zp->push_back(caloZ->at(i));
	        phi->push_back(caloPhi->at(i));
	        energy->push_back(ucmEDep->at(i));
        }
    }
    
    for(int i = 0; i < t0DwEDep->size(); ++i){
        if(corr[0][t0Z->at(i)][t0Phi->at(i)+8]==1 && t0DwEDep->at(i)>0){
            r->push_back(0);
	        zp->push_back(t0Z->at(i));
	        phi->push_back(t0Phi->at(i));
	        energy->push_back(t0DwEDep->at(i));
        }
        if(corr[1][t0Z->at(i)][t0Phi->at(i)+8]==1 && t0UpEDep->at(i)>0){
            r->push_back(1);
	        zp->push_back(t0Z->at(i));
	        phi->push_back(t0Phi->at(i));
	        energy->push_back(t0UpEDep->at(i));
        }
    }
    if(energy->size()>0)  tree->Fill();
      
  }
  tree->Write();
  file->Close();
  cout<<"out of events loop"<<endl;

}
