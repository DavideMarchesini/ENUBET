void merge(){


  string filename0="output_105.root";
  string filename1="output_104.root";
  string filename2="output_15.root";
  string output="output.root";
  cout<<"================= RUN is ON ===================== "<<endl;
  TFile *fc0 = new TFile(filename0.c_str());
  fc0->cd();
  TTree *fChainc0 = (TTree*)fc0->Get("tree");
  std::vector<int> *zp0 = nullptr;
  std::vector<int> *r0 = nullptr;
  std::vector<int> *phi0 = nullptr;
  std::vector<double> *LG0 = nullptr;
  std::vector<double> *HG0 = nullptr;
  fChainc0->SetBranchAddress("r", &r0);
  fChainc0->SetBranchAddress("zp", &zp0);
  fChainc0->SetBranchAddress("phi", &phi0);
  fChainc0->SetBranchAddress("LG", &LG0);
  fChainc0->SetBranchAddress("HG", &HG0);
  int n0=fChainc0->GetEntries();
  
  TFile *fc1 = new TFile(filename1.c_str());
  fc1->cd();
  TTree *fChainc1 = (TTree*)fc1->Get("tree");
  std::vector<int> *zp1 = nullptr;
  std::vector<int> *r1 = nullptr;
  std::vector<int> *phi1 = nullptr;
  std::vector<double> *LG1 = nullptr;
  std::vector<double> *HG1 = nullptr;
  fChainc1->SetBranchAddress("r", &r1);
  fChainc1->SetBranchAddress("zp", &zp1);
  fChainc1->SetBranchAddress("phi", &phi1);
  fChainc1->SetBranchAddress("LG", &LG1);
  fChainc1->SetBranchAddress("HG", &HG1);
  int n1=fChainc1->GetEntries();
  
  TFile *fc2 = new TFile(filename2.c_str());
  fc2->cd();
  TTree *fChainc2 = (TTree*)fc2->Get("tree");
  std::vector<int> *zp2 = nullptr;
  std::vector<int> *r2 = nullptr;
  std::vector<int> *phi2 = nullptr;
  std::vector<double> *LG2 = nullptr;
  std::vector<double> *HG2 = nullptr;
  fChainc2->SetBranchAddress("r", &r2);
  fChainc2->SetBranchAddress("zp", &zp2);
  fChainc2->SetBranchAddress("phi", &phi2);
  fChainc2->SetBranchAddress("LG", &LG2);
  fChainc2->SetBranchAddress("HG", &HG2);
  int n2=fChainc2->GetEntries();
  
  TFile *file = TFile::Open(output.c_str(), "RECREATE");
  TTree *tree = new TTree("tree", "A tree with a vector");
  std::vector<int> *zM = nullptr;
  std::vector<int> *rM = nullptr;
  std::vector<int> *phiM = nullptr;
  std::vector<double> *LGM = nullptr;
  std::vector<double> *HGM = nullptr;
  tree->Branch("r", &rM);
  tree->Branch("zp", &zM);
  tree->Branch("phi", &phiM);
  tree->Branch("LG", &LGM);
  tree->Branch("HG", &HGM);
  
  
  
  int k=0, triggerTS=0, nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<n0; jentry++) {
    Long64_t ientry = fChainc0->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChainc0->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000==0) cout<<"== EVENT "<<jentry<<" =="<<endl;
    rM = new std::vector<int>();
    zM = new std::vector<int>();
    phiM = new std::vector<int>();
    LGM = new std::vector<double>();
    HGM = new std::vector<double>();
    for(int i=0; i<r0->size(); i++){
      rM->push_back(r0->at(i));
      zM->push_back(zp0->at(i));
      phiM->push_back(phi0->at(i));
      LGM->push_back(LG0->at(i));
      HGM->push_back(HG0->at(i));
    }
    tree->Fill();
    delete rM;
    delete zM;
    delete phiM;
    delete LGM;
    delete HGM;
  }
  
  for (Long64_t jentry=0; jentry<n1; jentry++) {
    Long64_t ientry = fChainc1->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChainc1->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000==0) cout<<"== EVENT "<<jentry<<" =="<<endl;
    rM = new std::vector<int>();
    zM = new std::vector<int>();
    phiM = new std::vector<int>();
    LGM = new std::vector<double>();
    HGM = new std::vector<double>();
    for(int i=0; i<r1->size(); i++){
      rM->push_back(r1->at(i));
      zM->push_back(zp1->at(i));
      phiM->push_back(phi1->at(i));
      LGM->push_back(LG1->at(i));
      HGM->push_back(HG1->at(i));
    }
    tree->Fill();
    delete rM;
    delete zM;
    delete phiM;
    delete LGM;
    delete HGM;
  }
  
  for (Long64_t jentry=0; jentry<n2; jentry++) {
    Long64_t ientry = fChainc2->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChainc2->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000==0) cout<<"== EVENT "<<jentry<<" =="<<endl;
    rM = new std::vector<int>();
    zM = new std::vector<int>();
    phiM = new std::vector<int>();
    LGM = new std::vector<double>();
    HGM = new std::vector<double>();
    for(int i=0; i<r2->size(); i++){
      rM->push_back(r2->at(i));
      zM->push_back(zp2->at(i));
      phiM->push_back(phi2->at(i));
      LGM->push_back(LG2->at(i));
      HGM->push_back(HG2->at(i));
    }
    tree->Fill();
    delete rM;
    delete zM;
    delete phiM;
    delete LGM;
    delete HGM;
  }
  tree->Write();
  file->Close();
  fc0->Close(); 
  fc1->Close(); 
  fc2->Close(); 
}
