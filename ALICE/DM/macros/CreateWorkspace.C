using namespace RooFit;

void CreateWorkspace()
{
  gSystem->Exec("rm -f ws_raw.root ws_pT_int.root ws_C_int.root ws_Full.root");
  string w_name = "ws_raw.root";
  string w_pT_name = "ws_pT_int.root";
  string w_C_name = "ws_C_int.root";
  string w_Full_name = "ws_Full.root";

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  TString root_file_dir = "../ROOT_FILES/";
  TString root_file_name = "DM_LHC15o_pDCAcutAndPcut025GeV_tree_ptFine_5Mix";
  //TString root_file_name = "DM_LHC16b1_e2";
  TString full_file_name = root_file_dir + root_file_name + ".root";
  TFile *f = new TFile(full_file_name, "OLD");

  vector<string> C_bins = {"020", "2050", "5090"};
  //vector<string> pT_bins = {"23", "34", "45", "56", "67", "78", "89", "90"};

  vector<string> pT_bins = {"56", "67", "78", "89", "90"};

  vector<string> pair_bins = {"PM", "PP", "MM"};
  vector<string> dm_bins = {"", "Mix"};

  RooWorkspace w("w", root_file_name);
  for (auto i_dm = 0; i_dm < dm_bins.size(); ++i_dm)
  {
    for (auto i_pair = 0; i_pair < pair_bins.size(); ++i_pair)
    {
      for (auto i_c = 0; i_c < C_bins.size(); ++i_c)
      {
        for (auto i_pt = 0; i_pt < pT_bins.size(); ++i_pt)
        {
          string tree_name = "ftree" + dm_bins[i_dm] + pair_bins[i_pair] + C_bins[i_c] + pT_bins[i_pt];
          string data_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_C" + C_bins[i_c] + "_pT" + pT_bins[i_pt];
          cout << "Processing " << tree_name << "..." << endl;
          TTree *tree = (TTree *)f->Get(tree_name.c_str());
          RooRealVar tMass("tMass", "tMass", 1.5, 8);
          RooDataSet data(data_name.c_str(), data_name.c_str(), RooArgSet(tMass), Import(*tree));
          string nevents_name = "entries" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_C" + C_bins[i_c] + "_pT" + pT_bins[i_pt];
          RooRealVar nevents(nevents_name.c_str(), nevents_name.c_str(), data.numEntries());
          w.import(nevents);
          w.import(data);
        }
      }
    }
  }
  // w.Print();
  // cout << ">>> writing the workspace for the raw data...";
  // w.writeToFile(w_name.c_str());
  // cout << "done.\n";

  cout << "-------------------------------------------\n";

  // pT integrated datasets
  RooWorkspace w_pT("w_pT", root_file_name);
  for (auto i_dm = 0; i_dm < dm_bins.size(); ++i_dm)
  {
    for (auto i_pair = 0; i_pair < pair_bins.size(); ++i_pair)
    {
      for (auto i_c = 0; i_c < C_bins.size(); ++i_c)
      {
        RooDataSet *data_pT_int[pT_bins.size()];
        string data_pT_int_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_C" + C_bins[i_c];
        cout << "Creating " << data_pT_int_name << "..." << endl;
        for (auto i_pt = 0; i_pt < pT_bins.size(); ++i_pt)
        {
          string data_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_C" + C_bins[i_c] + "_pT" + pT_bins[i_pt];
          data_pT_int[i_pt] = (RooDataSet *)w.data(data_name.c_str());
        }
        RooDataSet data = RooDataSet(*data_pT_int[0], data_pT_int_name.c_str());
        for (auto i = 0; i < pT_bins.size(); ++i)
          data.append(*data_pT_int[i]);
        string nevents_name = "entries" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_C" + C_bins[i_c];
        RooRealVar nevents(nevents_name.c_str(), nevents_name.c_str(), data.numEntries());
        w_pT.import(nevents);
        w_pT.import(data);
      }
    }
  }
  // w_pT.Print();
  // cout << ">>> writing the workspace for the pT integrated data...";
  // w_pT.writeToFile(w_pT_name.c_str());
  // cout << "done.\n";

  cout << "-------------------------------------------\n";

  // centrality integrated datasets
  RooWorkspace w_C("w_C", root_file_name);
  for (auto i_dm = 0; i_dm < dm_bins.size(); ++i_dm)
  {
    for (auto i_pair = 0; i_pair < pair_bins.size(); ++i_pair)
    {
      for (auto i_pt = 0; i_pt < pT_bins.size(); ++i_pt)
      {
        RooDataSet *data_C_int[C_bins.size()];
        string data_C_int_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_pT" + pT_bins[i_pt];
        cout << "Creating " << data_C_int_name << "..." << endl;
        for (auto i_c = 0; i_c < C_bins.size(); ++i_c)
        {
          string data_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_C" + C_bins[i_c] + "_pT" + pT_bins[i_pt];
          data_C_int[i_c] = (RooDataSet *)w.data(data_name.c_str());
        }
        RooDataSet data = RooDataSet(*data_C_int[0], data_C_int_name.c_str());
        for (auto i = 0; i < C_bins.size(); ++i)
          data.append(*data_C_int[i]);
        string nevents_name = "entries" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_pT" + pT_bins[i_pt];
        RooRealVar nevents(nevents_name.c_str(), nevents_name.c_str(), data.numEntries());
        w_C.import(nevents);
        w_C.import(data);
      }
    }
  }
  w_C.Print();
  cout << ">>> writing the workspace for the centrality integrated data...";
  w_C.writeToFile(w_C_name.c_str());
  cout << "done.\n";

  cout << "-------------------------------------------\n";

  // Full datasets
  RooWorkspace w_Full("w_Full", root_file_name);
  for (auto i_dm = 0; i_dm < dm_bins.size(); ++i_dm)
  {
    for (auto i_pair = 0; i_pair < pair_bins.size(); ++i_pair)
    {

      string nametmp = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_pT" + pT_bins[0]; 
      RooDataSet *datatmp = (RooDataSet *)w_C.data(nametmp.c_str());
      for (auto i_pt = 1; i_pt < pT_bins.size(); ++i_pt) 
      {
        nametmp = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_pT" + pT_bins[i_pt];
        RooDataSet *datatmp2 = (RooDataSet *)w_C.data(nametmp.c_str());
        datatmp->append(*datatmp2);
      }
      nametmp = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_full";
      datatmp->SetName(nametmp.c_str());
      w_Full.import(*datatmp);
    }
  }
  w_Full.Print();
  cout << ">>> writing the workspace for the full dataset...";
  w_Full.writeToFile(w_Full_name.c_str());
  cout << "done.\n";
  gApplication->Terminate();
}
