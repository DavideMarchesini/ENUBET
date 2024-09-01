void sum0() {
    auto open = [](const char* filename, const char* leaf) {
        std::vector<float> data;
        TFile file(filename);
        if (!file.IsOpen()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return data;
        }
        TTree* tree = (TTree*)file.Get("tree");
        if (!tree) {
            std::cerr << "Error: Could not find tree in file " << filename << std::endl;
            return data;
        }
        std::vector<float>* branch = nullptr;
        tree->SetBranchAddress(leaf, &branch);
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);
            data.insert(data.end(), branch->begin(), branch->end());
        }
        file.Close();
        return data;
    };

    auto get_energy = [](int x, int y, int layer, const std::vector<float>& caloPhi, const std::vector<float>& caloZ, const std::vector<float>& caloR, const std::vector<float>& ucmEDep) {
        float e = 0;
        int k=0;
        for (size_t i = 0; i < ucmEDep.size(); ++i) {
            if (caloZ[i] == x && caloPhi[i] == y && caloR[i] == layer) {
                e += ucmEDep[i];
                k += 1;
            }
        }
        return k;
    };

    const char* filename = "output.root";
    std::vector<float> r = open(filename, "r");
    std::vector<float> zp = open(filename, "zp");
    std::vector<float> phi = open(filename, "phi");
    std::vector<float> energy = open(filename, "LG");

    std::string p = "sum_data.txt";
    std::ofstream df(p);

    for (int k = 0; k < 5; ++k) {
        for (int p = 0; p < 15; ++p) {
            for (int t = 0; t < 26; ++t) {
                float e_sum = get_energy(p, t - 8, k, phi, zp, r, energy);
                cout << k << "  " << p << "  " << t - 8 << "  " << e_sum << endl;
                if(e_sum>0)
                  df << k << "  " << p << "  " << t - 8 << "  " << e_sum << "\n";
            }
        }
    }

    df.close();
}

