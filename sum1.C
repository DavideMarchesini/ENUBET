void sum1(){
    string p = "sum_data.txt";
    std::ifstream inputFile(p);
    std::string line;
    std::vector<float> R_data;
    std::vector<float> Z_data;
    std::vector<float> PHI_data;
    std::vector<float> E_data;
    int n0=0, n1 = 0;
    float mean, mean0 = 0, mean1=0; 
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        double col1, col2, col3, col4;
        if (!(iss >> col1 >> col2 >> col3 >> col4)) {
            std::cerr << "Error reading line: " << line << std::endl;
            continue; 
        }
        int R = static_cast<int>(col1); 
        int Z = static_cast<int>(col2); 
        int PHI = static_cast<int>(col3)+8; 
        float E = static_cast<float>(col4);
        R_data.push_back(R);
        Z_data.push_back(Z);
        PHI_data.push_back(PHI);
        E_data.push_back(E);
        if(R==0 || R==1){
          mean0 += E;
          n0++;
        }
        else{
          mean1 += E;
          n1++;
        }
    }
    mean0 = mean0 / n0;
    mean1 = mean1 / n1;
    cout<<mean0<<endl<<mean1<<endl;
    inputFile.close();

    p = "u.txt";
    std::ofstream df(p);
    for(int k = 0; k < 5; k++){
        TCanvas* c = new TCanvas("c", "Canvas", 800, 600);
        
        c->SetRightMargin(0.15);  
        
        string t = "Layer " + std::to_string(k);
        TH2D* h = new TH2D("h", (t + "; Z; PHI").c_str(), 15, 0, 15, 26, 0, 26);
        float j;
        for (size_t i = 0; i < R_data.size(); ++i) {
            if(R_data[i] == k){
                if(k==0 || k==1) mean=mean0;
                else mean=mean1;
                j = (E_data[i] - mean) / (E_data[i] + mean);
                //j=E_data[i];
                h->Fill(Z_data[i], PHI_data[i], j);
                df << k << "  " << Z_data[i] << "  " << PHI_data[i] << "  " << j << "\n";   
            }
        }
        gStyle->SetOptStat(0);
        h->SetContour(99);
        h->Draw("COL1Z");

        for (int i = 0; i <= 26; ++i) {
            TLine* line = new TLine(8, i, 15, i);
            line->SetLineColor(kBlack);
            line->Draw();
        }
        for (int i = 0; i <= 10; ++i) {
            TLine* line = new TLine(0, 8 + i, 15, 8 + i);
            line->SetLineColor(kBlack);
            line->Draw();
        }
        for (int i = 0; i <= 8; ++i) {
            TLine* line = new TLine(8 + i, 0, 8 + i, 26);
            line->SetLineColor(kBlack);
            line->Draw();
        }
        for (int i = 0; i <= 11; ++i) {
            TLine* line = new TLine(i, 8, i, 18);
            line->SetLineColor(kBlack);
            line->Draw();
        }
        
        string result = "sum/u_" + std::to_string(k) + ".png";
        c->SaveAs(result.c_str());
    } 
    df.close(); 
}

