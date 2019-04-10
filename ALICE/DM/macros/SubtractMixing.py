#! /usr/bin/env python
import os, sys, getopt, multiprocessing, copy, math, itertools
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom, Double
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis, TH1D, TGaxis
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TMath
from ROOT import RooFit, RooRealVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooExtendPdf, RooAddPdf, RooProdPdf, RooNumConvPdf, RooFFTConvPdf

from drawUtils import *

########## OPTIONS ##########

import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-b', '--batch',      action='store_true', default=False,           dest='batch'                    )
parser.add_option('-v', '--verbose',    action='store_true', default=False,           dest='verbose'                  )
parser.add_option('-d', '--datapath',   action='store',      default='data/2015/',    dest='datapath',  type='string' )
parser.add_option('-w', '--wspath',     action='store',      default='../workspaces/',   dest='wspath',    type='string' )
parser.add_option('-D', '--do_draw',    action='store_true', default=False,           dest='do_draw'                  )
parser.add_option('-n', '--n_bins',     action='store',      default=130,             dest='n_bins',    type='int'    )
(options, args) = parser.parse_args()
if options.batch: gROOT.SetBatch(True)


def ComputeR(mixPM, mixPP, mixMM):
  return mixPM / (2 * math.sqrt(mixPP * mixMM))

def ComputeF(R, dataPP, dataMM, mixPM):
  return 2 * R * math.sqrt(dataPP * dataMM) / mixPM

def SubtractMixing(F, h_data, h_mixing):
  subtraction = h_data.Clone()
  h_mixing.Scale(F)
  subtraction.Add(h_mixing, -1)
  return subtraction

def MergeRooDataSet(dsvect):
  for ids in dsvect:
    dsvect[0].append(ids)

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)

TGaxis.SetMaxDigits(2)
plot_file_name  = ''
plot_out_folder = '../plots/'
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

C_bins    = ["020", "2050", "5090"]
pT_bins   = ["23", "34", "45", "56", "67", "78", "89", "90"]
pair_bins = ["PM", "PP", "MM"]
dm_bins   = ["", "Mix"]

raw_data      = {}
h_raw_data    = {}
data_entries  = {}
h_data_mix    = {}
h_data_mix_roo= {}
weight        = {}
data_mix      = {}

print ("=== Opening the ws files - Integrated over Centrality")

file = TFile("../workspaces/ws_C_int.root",'read');
w = file.Get("w_C");
print ("    Workspace w corresponds to file:", w.GetTitle())
out_w = RooWorkspace("processed_w_C")
w.Print()

mass = w.var("tMass")
mass_range = (1.5, 8.)
binsmass = RooBinning(options.n_bins, mass_range[0], mass_range[1])
mass.setRange('rangemass',mass_range[0], mass_range[1])

# All bins - Integrated over Centrality
for it in itertools.product(dm_bins,pair_bins,pT_bins):
  data_name     = "data{}_{}_pT{}".format(it[0],it[1],it[2])
  h_data_name   = 'h_%s'%data_name
  raw_data[data_name] = w.data(data_name)
  data_entries[data_name] = raw_data[data_name].numEntries()
  h_raw_data[data_name]   = raw_data[data_name].createHistogram(mass, mass, options.n_bins, options.n_bins, "", data_name).ProjectionX(data_name)
  if options.do_draw:
    c = TCanvas("c_full_%s"%data_name, "c_full_%s"%data_name, 800, 600)
    h_raw_data[data_name].Draw("PE")
    h_raw_data[data_name].GetYaxis().SetTitle("Events / x GeV")
    h_raw_data[data_name].GetXaxis().SetTitle("Mass #mu#mu (GeV)")
    setHistStyle(h_raw_data[data_name], 1.1)
    h_raw_data[data_name].GetYaxis().SetTitleOffset(h_raw_data[data_name].GetYaxis().GetTitleOffset()*0.5)
    h_raw_data[data_name].GetYaxis().SetTitleOffset(1.1)
    h_raw_data[data_name].GetYaxis().SetRangeUser(0., h_raw_data[data_name].GetMaximum()*1.2)
    drawALICE("","Internal")
    drawRegion(data_name)
    c.SaveAs(plot_out_folder+data_name+'.pdf')
    c.SaveAs(plot_out_folder+data_name+'.png')

# data-mixing --- All pT bins
for it in pT_bins:
  mixPM   = data_entries["dataMix_PM_pT%s"%it]
  mixPP   = data_entries["dataMix_PP_pT%s"%it]
  mixMM   = data_entries["dataMix_MM_pT%s"%it]
  dataPP  = data_entries["data_PP_pT%s"%it]
  dataMM  = data_entries["data_MM_pT%s"%it]
  R = ComputeR(mixPM, mixPP, mixMM)
  F = ComputeF(R, dataPP, dataMM, mixPM)
  h_data_mix[it]     = SubtractMixing(F, h_raw_data["data_PM_pT%s"%it], h_raw_data["dataMix_PM_pT%s"%it])
  h_data_mix_roo[it] = RooDataHist("h_data_mix_pT%s"%it, "", RooArgList(mass), RooFit.Import(h_data_mix[it]))
  if options.do_draw:
    c = TCanvas("c_data_mix_pT%s"%it, "c_data_mix_pT%s"%it, 800, 600)
    h_data_mix[it].Draw("PE")
    drawALICE("","Internal")
    drawRegion('pT%s'%it)
    c.SaveAs(plot_out_folder+'data_mix_pT%s.pdf'%it)
    c.SaveAs(plot_out_folder+'data_mix_pT%s.png'%it)
  getattr(out_w, 'import')(h_data_mix_roo[it])
  weight[it] = RooRealVar("weight_pT%s"%it, "weight_pT%s"%it, -1, 1)
  weight[it].setVal(1)
  data_mix[it] = RooDataSet("data_mix_pT%s"%it, "data_mix_pT%s"%it, RooArgSet(mass, weight[it]), RooFit.Import(raw_data["data_PM_pT%s"%it]), RooFit.WeightVar(weight[it]))
  weight[it].setVal(-F)
  mixtmp = RooDataSet("mixtmp", "mixtmp", RooArgSet(mass, weight[it]), RooFit.Import(raw_data["dataMix_PM_pT%s"%it]), RooFit.WeightVar(weight[it]))
  data_mix[it].append(mixtmp)
  getattr(out_w, 'import')(data_mix[it])
  getattr(out_w, 'import')(raw_data["data_PM_pT%s"%it])
  out_w.writeToFile('{}{}.root'.format(options.wspath,out_w.GetName()))

out_w.Print()





print ("=== Opening the ws files - Fully Integrated")

file_Full = TFile("../workspaces/ws_Full.root",'read');
w_Full = file_Full.Get("w_Full");
print ("    Workspace w corresponds to file:", w_Full.GetTitle())
out_w_Full = RooWorkspace("processed_w_Full")
w_Full.Print()

mass = w_Full.var("tMass")
mass_range = (1.5, 8.)
binsmass = RooBinning(options.n_bins, mass_range[0], mass_range[1])
mass.setRange('rangemass',mass_range[0], mass_range[1])

# All bins - Integrated over Centrality
for it in itertools.product(dm_bins,pair_bins):
  data_name     = "data{}_{}".format(it[0],it[1])
  h_data_name   = 'h_%s'%data_name
  raw_data[data_name] = w_Full.data(data_name)
  data_entries[data_name] = raw_data[data_name].numEntries()
  h_raw_data[data_name]   = raw_data[data_name].createHistogram(mass, mass, options.n_bins, options.n_bins, "", data_name).ProjectionX(data_name)
  if options.do_draw:
    c = TCanvas("c_full_%s"%data_name, "c_full_%s"%data_name, 800, 600)
    h_raw_data[data_name].Draw("PE")
    h_raw_data[data_name].GetYaxis().SetTitle("Events / x GeV")
    h_raw_data[data_name].GetXaxis().SetTitle("Mass #mu#mu (GeV)")
    setHistStyle(h_raw_data[data_name], 1.1)
    h_raw_data[data_name].GetYaxis().SetTitleOffset(h_raw_data[data_name].GetYaxis().GetTitleOffset()*0.5)
    h_raw_data[data_name].GetYaxis().SetTitleOffset(1.1)
    h_raw_data[data_name].GetYaxis().SetRangeUser(0., h_raw_data[data_name].GetMaximum()*1.2)
    drawALICE("","Internal")
    drawRegion(data_name)
    c.SaveAs(plot_out_folder+data_name+'.pdf')
    c.SaveAs(plot_out_folder+data_name+'.png')

# data-mixing --- All pT bins
mixPM   = data_entries["dataMix_PM"]
mixPP   = data_entries["dataMix_PP"]
mixMM   = data_entries["dataMix_MM"]
dataPP  = data_entries["data_PP"]
dataMM  = data_entries["data_MM"]
R = ComputeR(mixPM, mixPP, mixMM)
F = ComputeF(R, dataPP, dataMM, mixPM)
h_data_mix     = SubtractMixing(F, h_raw_data["data_PM"], h_raw_data["dataMix_PM"])
h_data_mix_roo = RooDataHist("h_data_mix_pT", "", RooArgList(mass), RooFit.Import(h_data_mix))
if options.do_draw:
  c = TCanvas("c_data_mix", "c_data_mix", 800, 600)
  h_data_mix.Draw("PE")
  drawALICE("","Internal")
  drawRegion('')
  c.SaveAs(plot_out_folder+'data_mix.pdf')
  c.SaveAs(plot_out_folder+'data_mix.png')
getattr(out_w_Full, 'import')(h_data_mix_roo)
weight = RooRealVar("weight", "weight", -1, 1)
weight.setVal(1)
data_mix = RooDataSet("data_mix", "data_mix", RooArgSet(mass, weight), RooFit.Import(raw_data["data_PM"]), RooFit.WeightVar(weight))
weight.setVal(-F)
mixtmp = RooDataSet("mixtmp", "mixtmp", RooArgSet(mass, weight), RooFit.Import(raw_data["dataMix_PM"]), RooFit.WeightVar(weight))
data_mix.append(mixtmp)
getattr(out_w_Full, 'import')(data_mix)
getattr(out_w_Full, 'import')(raw_data["data_PM"])
out_w_Full.writeToFile('{}{}.root'.format(options.wspath,out_w_Full.GetName()))

out_w_Full.Print()





















# void SubtractMixing(int n_bins = 130)
# {
#   gStyle->SetOptStat("");
#   TGaxis::SetMaxDigits(2);
#   TString plot_file_name;
#   TString plot_out_folder = "Plots/";
#   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
#   cout << ">> Opening the ws files" << endl;

#   TFile *file = TFile::Open("../workspaces/ws_C_int.root");
#   RooWorkspace *w = (RooWorkspace *)file->Get("w_C");
#   w->Print();
#   RooRealVar *mass = w->var("tMass");
#   //mass->setRange(3, 8);

#   RooWorkspace out_w("processed_w_Full");

#   cout << ">> Workspace w corresponds to file: " << w->GetTitle() << endl;

#   vector<string> pT_bins = {"56", "67", "78", "89", "90"};
#   vector<string> pair_bins = {"PM", "PP", "MM"};
#   vector<string> dm_bins = {"", "Mix"};

#   map<string, RooDataSet *> raw_data;
#   map<string, TH1F *> h_raw_data;
#   map<string, double> data_entries;
#   TH1F *h_data_mix;

#   vector<string> data_type = {"Data", "Mixing"};
#   vector<string> data_pair = {"#mu^{#pm}#mu^{#mp}", "#mu^{+}#mu^{+}", "#mu^{-}#mu^{-}"};

#   cout << ">> Reading the workspace and creating the histograms" << endl;
#   for (auto i_pt = 0; i_pt < pT_bins.size(); ++i_pt)
#   {
#     for (auto i_dm = 0; i_dm < dm_bins.size(); ++i_dm)
#     {
#       for (auto i_pair = 0; i_pair < pair_bins.size(); ++i_pair)
#       {
#         string data_name = "data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_pT" + pT_bins[i_pt];
#         string h_data_name = "h_data" + dm_bins[i_dm] + "_" + pair_bins[i_pair] + "_pT" + pT_bins[i_pt];
#         cout << data_name << endl;
#         raw_data[data_name] = (RooDataSet *)w->data(data_name.c_str());
#         data_entries[data_name] = raw_data[data_name]->numEntries();
#         TH1F *tmp = (TH1F *)raw_data[data_name]->createHistogram("tMass", n_bins);
#         tmp->SetName(h_data_name.c_str());
#         h_raw_data[data_name] = tmp;
#       }
#     }
#   }

#   //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

#   // plot full dataset
#   if (do_draw)
#   {
#     TCanvas *c_full = new TCanvas("c_full", "c_full", 1900, 1000);
#     c_full->Divide(3, 2);
#     for (auto i_pt = 0; i_pt < pT_bins.size(); ++i_pt)
#     {
#       c_full->cd(i_pt + 1);
#       string data_name = "data_PM_pT" + pT_bins[i_pt];
#       h_raw_data[data_name]->SetName("h_raw_data");
#       h_raw_data[data_name]->SetTitle("Data");
#       h_raw_data[data_name]->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
#       h_raw_data[data_name]->SetLineColor(1);
#       h_raw_data[data_name]->SetMarkerStyle(20);
#       h_raw_data[data_name]->Draw();
#     }
#   }

#   return;

#   // plot full dataset (mixing)
#   if (do_draw)
#   {
#     TCanvas *c_mix = new TCanvas("c_mix", "c_mix", 1900, 1000);
#     c_mix->Divide(3, 2);
#     for (auto i_pt = 0; i_pt < pT_bins.size(); ++i_pt)
#     {
#       c_mix->cd(i_pt + 1);
#       string data_name = "dataMix_PM_pT" + pT_bins[i_pt];
#       h_raw_data[data_name]->SetName("h_raw_data_mix");
#       h_raw_data[data_name]->SetTitle("Mixing");
#       h_raw_data[data_name]->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
#       h_raw_data[data_name]->SetLineColor(1);
#       h_raw_data[data_name]->SetMarkerStyle(20);
#       h_raw_data[data_name]->Draw();
#     }
#   }

#   // data - mixing
#   int mixPM = data_entries["dataMix_PM_full"];
#   int mixPP = data_entries["dataMix_PP_full"];
#   int mixMM = data_entries["dataMix_MM_full"];
#   int dataPP = data_entries["data_PP_full"];
#   int dataMM = data_entries["data_MM_full"];
#   double R = ComputeR(mixPM, mixPP, mixMM);
#   double F = ComputeF(R, dataPP, dataMM, mixPM);
#   h_data_mix = SubtractMixing(n_bins, F, raw_data["data_PM_full"], raw_data["dataMix_PM_full"]);
#   RooDataHist *h_data_mix_roo = new RooDataHist("h_data_mix", "", *mass, h_data_mix);
#   out_w.import(*h_data_mix_roo);
#   string h_data_mix_name = "Histogram of Data (Mixing subtracted)";
#   h_data_mix->SetTitle(h_data_mix_name.c_str());
#   if (do_draw)
#   {
#     new TCanvas;
#     h_data_mix->Draw();
#   }

#   RooRealVar weight("weight", "weight", -1, 1);
#   weight.setVal(1);
#   RooDataSet *data_mix = new RooDataSet("data_mix", "data_mix", RooArgSet(*mass, weight), Import(*raw_data["data_PM_full"]), WeightVar(weight));

#   weight.setVal(-F);
#   RooDataSet *mixtmp = new RooDataSet("mixtmp", "mixtmp", RooArgSet(*mass, weight), Import(*raw_data["dataMix_PM_full"]), WeightVar(weight));

#   data_mix->append(*mixtmp);
#   data_mix->Print();
#   out_w.import(*data_mix);
#   out_w.import(*raw_data["data_PM_full"]);

#   new TCanvas;
#   RooPlot *plot_data_mix = mass->frame(Name("plot_data_mix"), Title("Data (Mixing subtracted)"), Bins(n_bins));
#   data_mix->plotOn(plot_data_mix);
#   plot_data_mix->Draw();

#   out_w.Print();
#   out_w.writeToFile("processed_w_Full.root");
#   if (!do_draw)
#     gApplication->Terminate();
# }
