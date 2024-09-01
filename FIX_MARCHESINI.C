// D. Marchesini July 2024
#include "MappingScript.h"
void FIX_MARCHESINI(){
  /*Only sector to modify by defining the parameters
  - p is the file generated as mean and std of the pedestal
  - filename is the input file given by Janus
  - output is the name of the final output file
  - fileInput is the calibration of channels given by board/channel/weight
  - numberOfBoards is the nunmber of boards
  - thresholdTrigger is the time coincidence in mu s to impose
  - sigma is the number of std of the pedestal to which it eliminates the value
  - numberOfAnodes is the number of channels connected to the board
  */
  string p="pedestal_data.txt";
  string filename="Run105_list.root";
  string output="output_105.root";
  string fileInput = "calibrationMarchesini.txt";
  string fileInputMap="mappingMarchesini.txt";
  int numberOfBoards = 13;
  double thresholdTrigger=0.010;
  double sigma=4.0;  
  int numberOfAnodes = 64;
  int majority=3;
  
  //Mapping Part and Equalization
  const int NNR=5;
  const int NNZ=15;
  const int NNPH=25;
  std::map<std::vector<int>, std::vector<int>> dMap1=MappingAnodeCord(fileInputMap.c_str());
  double mRR[numberOfAnodes * numberOfBoards];
  double mZZ[numberOfAnodes * numberOfBoards];
  double mPP[numberOfAnodes * numberOfBoards];
  ifstream in (fileInput.c_str());
  int b1,a1;
  double mpv1, x1, y1, mpv_err1, chi21, coeff1;
  vector <double> weight1;
  std::map<int,double> map;
  in.ignore(1000, '\n'); 
  while (true) {
    in >> b1 >> a1 >> coeff1;
    int f1 = b1 * numberOfAnodes + a1;
    map[f1] = coeff1;
    mZZ[f1]=dMap1[{b1,a1}][2];
    mRR[f1]=dMap1[{b1,a1}][0];
    mPP[f1]=dMap1[{b1,a1}][1];
    if (in.eof() == true) break;
    weight1.push_back( coeff1 );
  }
  in.close();
  
  //Begin work on the root file, so decalaration of all the varaibles in input and output
  cout<<"================= RUN is ON ===================== "<<endl;
  TFile *fc = new TFile(filename.c_str());
  if (!fc || fc->IsZombie()) {
        std::cerr << "Error opening file " << std::endl;
        return;
    }
  fc->cd();
  TTree *fChainc = (TTree*)fc->Get("t");
  if (!fChainc) {
        std::cerr << "Error: Tree 't' not found in the file" << std::endl;
        fc->Close();
        return;
    }
  Short_t lg[numberOfAnodes];
  Short_t hg[numberOfAnodes];
  UChar_t boardID;
  Double_t triggerTimeStamp;
  fChainc->SetBranchAddress("lg", lg);
  fChainc->SetBranchAddress("hg", hg);
  fChainc->SetBranchAddress("boardID", &boardID);
  fChainc->SetBranchAddress("triggerTimeStamp", &triggerTimeStamp); 
  TFile *file = TFile::Open(output.c_str(), "RECREATE");
  TTree *tree = new TTree("tree", "A tree with a vector");
  std::vector<int> *zp = nullptr;
  std::vector<int> *r = nullptr;
  std::vector<int> *phi = nullptr;
  std::vector<double> *LG = nullptr;
  std::vector<double> *HG = nullptr;
  double Timestamp;
  tree->Branch("r", &r);
  tree->Branch("zp", &zp);
  tree->Branch("phi", &phi);
  tree->Branch("LG", &LG);
  tree->Branch("HG", &HG);
  tree->Branch("Timestamp", &Timestamp);
  std::vector<double> fersHG(numberOfAnodes * numberOfBoards, 0);
  std::vector<double> fersLG(numberOfAnodes * numberOfBoards, 0);
  std::vector<double> thresholdLG(numberOfAnodes * numberOfBoards, 0);
  std::vector<double> thresholdHG(numberOfAnodes * numberOfBoards, 0);
  int k=0, triggerTS=0, nbytes = 0, nb = 0, nentries=fChainc->GetEntries();
  
  //Calculation of the pedestal value for all channels for the low gain and high gain
  std::ofstream df(p);
  for (int i = 0; i < numberOfBoards; i++) {
      for (int j = 0; j < numberOfAnodes; j++) {
          TH1F *hist = new TH1F("hist", "hist", 100, 10, 150);
          TF1 *fitFunc = new TF1("fitFunc", "gaus", 10, 150);
          for (Long64_t jentry = 0; jentry < nentries; jentry++) {
            Long64_t ientry = fChainc->LoadTree(jentry);
            nb = fChainc->GetEntry(jentry);   
            if (boardID == i)  
              hist->Fill(lg[j]);          
          }
          fitFunc->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
          hist->Fit(fitFunc, "Q");
          if (fitFunc->GetParameter(0) != 0) 
            df << i << "  " << j << "  " << fitFunc->GetParameter(1) << "  " << fitFunc->GetParameter(2) << "\n";   
    }
  }
  for (int i = 0; i < numberOfBoards; i++) {
      for (int j = 0; j < numberOfAnodes; j++) {
          TH1F *hist = new TH1F("hist", "hist", 100, 10, 150);
          TF1 *fitFunc = new TF1("fitFunc", "gaus", 10, 150);
          for (Long64_t jentry = 0; jentry < nentries; jentry++) {
            Long64_t ientry = fChainc->LoadTree(jentry);
            nb = fChainc->GetEntry(jentry);   
            if (boardID == i)  
              hist->Fill(hg[j]);          
          }
          fitFunc->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
          hist->Fit(fitFunc, "Q");
          if (fitFunc->GetParameter(0) != 0) 
            df << i << "  " << j << "  " << fitFunc->GetParameter(1) << "  " << fitFunc->GetParameter(2) << "\n";   
    }
  }
  df.close();
  
  //Reading the written values of the txt file and put them in a vector
  std::vector<double> dataLG(numberOfBoards * numberOfAnodes, 0.0);
  std::vector<double> data1LG(numberOfBoards * numberOfAnodes, 0.0);
  std::vector<double> dataHG(numberOfBoards * numberOfAnodes, 0.0);
  std::vector<double> data1HG(numberOfBoards * numberOfAnodes, 0.0);
  std::ifstream inputFile(p);
  if (!inputFile.is_open()) {
      std::cerr << "Error opening input file: " << p << std::endl;
      return 1;
  }
  std::string line;
  int v=0;
  while (std::getline(inputFile, line)) {
      std::istringstream iss(line);
      double col1, col2, col3, col4;
      if (!(iss >> col1 >> col2 >> col3 >> col4)) {
          std::cerr << "Error reading line: " << line << std::endl;
          continue; 
      }
      int board = static_cast<int>(col1); 
      int anode = static_cast<int>(col2); 
      if (board < numberOfBoards && anode < numberOfAnodes) {
          int index = board * numberOfAnodes + anode;
          if(v<numberOfBoards*numberOfAnodes){
            dataLG[index] = col3;
            data1LG[index]= sigma*col4;
          }else{
             dataHG[index] = col3;
             data1HG[index]= sigma*col4;
          }
        }else{
           std::cerr << "Index out of bounds: board=" << board << ", anode=" << anode << std::endl;
        }
      v++;
  }
  inputFile.close();
  
  //Cycle over all the events and disposition in the new root file
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = fChainc->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChainc->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000==0) cout<<"== EVENT "<<jentry<<" =="<<endl;
    r = new std::vector<int>();
    zp = new std::vector<int>();
    phi = new std::vector<int>();
    LG = new std::vector<double>();
    HG = new std::vector<double>();
    if(k==0){ 
      std::fill(fersLG.begin(), fersLG.end(), 0);
      std::fill(fersHG.begin(), fersHG.end(), 0);
    }
    for(int i=0; i < numberOfAnodes; ++i){
      fersLG.at(i+boardID*numberOfAnodes)=lg[i]-dataLG[i+boardID*numberOfAnodes];
      thresholdLG.at(i+boardID*numberOfAnodes)=data1LG[i+boardID*numberOfAnodes];
      fersHG.at(i+boardID*numberOfAnodes)=hg[i]-dataHG[i+boardID*numberOfAnodes];
      thresholdHG.at(i+boardID*numberOfAnodes)=data1HG[i+boardID*numberOfAnodes];
    }
    for(const auto &v : map){
      int fersIndex = v.first;
      double weight = v.second;
      double dEeq=fersLG.at(fersIndex)/weight;
      double dEeq1=fersHG.at(fersIndex)/weight;
      if (dEeq>thresholdLG.at(fersIndex))
      // && dEeq1>thresholdHG.at(fersIndex))
      {
	r->push_back(mRR[fersIndex]);
	zp->push_back(mZZ[fersIndex]);
	phi->push_back(mPP[fersIndex]);
	LG->push_back(dEeq);
	HG->push_back(dEeq1);
      }
    }    
    if(triggerTimeStamp-triggerTS<thresholdTrigger){
      k++;
    }
    else{
      k=0;
    }
    triggerTS=triggerTimeStamp;
    if(k==0){
      Timestamp=triggerTimeStamp;
      if(r->size()>majority)
        tree->Fill();
      delete r;
      delete zp;
      delete phi;
      delete LG;
      delete HG;
    }
  }
  tree->Write();
  file->Close();
  fc->Close(); 
}
