using namespace RooFit;
using namespace RooStats;
bool do_draw = true;

double ComputeR(double mixPM, double mixPP, double mixMM)
{
  return mixPM / (2 * sqrt(mixPP * mixMM));
}

double ComputeF(double R, double dataPP, double dataMM, double mixPM)
{
  return (2 * R * sqrt(dataPP * dataMM)) / mixPM;
}

TH1F *SubtractMixing(int nBins, double F, RooDataSet *data, RooDataSet *mixing)
{
  TH1F *h_data = (TH1F *)data->createHistogram("tMass", nBins);
  TH1F *h_mixing = (TH1F *)mixing->createHistogram("tMass", nBins);
  TH1F *subtraction = (TH1F *)h_data->Clone();
  h_mixing->Scale(F);
  subtraction->Add(h_mixing, -1);
  subtraction->SetLineColor(1);
  subtraction->SetMarkerStyle(20);
  return subtraction;
}

void MergeRooDataSet(vector<RooDataSet *> v)
{
  for (auto i = 1; i < v.size(); ++i)
    v[0]->append(*v[i]);
}

void ViewFullDataset(int n_bins = 130)
{
  gStyle->SetOptStat("");
  TGaxis::SetMaxDigits(2);
  TString plot_file_name;
  TString plot_out_folder = "Plots/";
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  cout << ">> Opening the ws files" << endl;

  TFile *file = TFile::Open("../workspaces/ws_Full.root");
  RooWorkspace *w = (RooWorkspace *)file->Get("w_Full");
  RooRealVar *mass = w->var("tMass");
  //mass->setRange(3, 8);

  RooWorkspace out_w("processed_w_Full");

  cout << ">> Workspace w corresponds to file: " << w->GetTitle() << endl;

  vector<string> pair_bins = {"PM", "PP", "MM"};
  vector<string> dm_bins = {"", "Mix"};

  map<string, RooDataSet *> raw_data;
  map<string, TH1F *> h_raw_data;
  map<string, double> data_entries;
  TH1F *h_data_mix;

  vector<string> data_type = {"Data", "Mixing"};
  vector<string> data_pair = {"#mu^{#pm}#mu^{#mp}", "#mu^{+}#mu^{+}", "#mu^{-}#mu^{-}"};

  cout << ">> Reading the workspace and creating the histograms" << endl;
  for (auto i_dm = 0; i_dm < dm_bins.size(); ++i_dm)
  {
    for (auto i_pair = 0; i_pair < pair_bins.size(); ++i_pair)
    {
      string data_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_full";
      string h_data_name = "h_data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_full";
      cout << data_name << endl;
      raw_data[data_name] = (RooDataSet *)w->data(data_name.c_str());
      data_entries[data_name] = raw_data[data_name]->numEntries();
      TH1F *tmp = (TH1F *)raw_data[data_name]->createHistogram("tMass", n_bins);
      tmp->SetName(h_data_name.c_str());
      h_raw_data[data_name] = tmp;
    }
  }

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // plot full dataset
  if (do_draw)
  {
    new TCanvas;
    h_raw_data["data_PM_full"]->SetName("h_raw_data");
    h_raw_data["data_PM_full"]->SetTitle("Data");
    h_raw_data["data_PM_full"]->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    h_raw_data["data_PM_full"]->SetLineColor(1);
    h_raw_data["data_PM_full"]->SetMarkerStyle(20);
    h_raw_data["data_PM_full"]->Draw();
  }

  // plot full dataset (mixing)
  if (do_draw)
  {
    new TCanvas;
    h_raw_data["dataMix_PM_full"]->SetName("h_raw_data_mix");
    h_raw_data["dataMix_PM_full"]->SetTitle("Mixing");
    h_raw_data["dataMix_PM_full"]->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    h_raw_data["dataMix_PM_full"]->SetLineColor(1);
    h_raw_data["dataMix_PM_full"]->SetMarkerStyle(20);
    h_raw_data["dataMix_PM_full"]->Draw();
  }

  // data - mixing
  int mixPM = data_entries["dataMix_PM_full"];
  int mixPP = data_entries["dataMix_PP_full"];
  int mixMM = data_entries["dataMix_MM_full"];
  int dataPP = data_entries["data_PP_full"];
  int dataMM = data_entries["data_MM_full"];
  double R = ComputeR(mixPM, mixPP, mixMM);
  double F = ComputeF(R, dataPP, dataMM, mixPM);
  h_data_mix = SubtractMixing(n_bins, F, raw_data["data_PM_full"], raw_data["dataMix_PM_full"]);
  RooDataHist *h_data_mix_roo = new RooDataHist("h_data_mix", "", *mass, h_data_mix);
  out_w.import(*h_data_mix_roo);
  string h_data_mix_name = "Histogram of Data (Mixing subtracted)";
  h_data_mix->SetTitle(h_data_mix_name.c_str());
  if (do_draw)
  {
    new TCanvas;
    h_data_mix->Draw();
  }

  RooRealVar weight("weight", "weight", -1, 1);
  weight.setVal(1);
  RooDataSet *data_mix = new RooDataSet("data_mix", "data_mix", RooArgSet(*mass, weight), Import(*raw_data["data_PM_full"]), WeightVar(weight));

  weight.setVal(-F);
  RooDataSet *mixtmp = new RooDataSet("mixtmp", "mixtmp", RooArgSet(*mass, weight), Import(*raw_data["dataMix_PM_full"]), WeightVar(weight));

  data_mix->append(*mixtmp);
  data_mix->Print();
  out_w.import(*data_mix);
  out_w.import(*raw_data["data_PM_full"]);

  new TCanvas;
  RooPlot *plot_data_mix = mass->frame(Name("plot_data_mix"), Title("Data (Mixing subtracted)"), Bins(n_bins));
  data_mix->plotOn(plot_data_mix);
  plot_data_mix->Draw();

  out_w.Print();
  out_w.writeToFile("processed_w_Full.root");
  if (!do_draw)
    gApplication->Terminate();
}
