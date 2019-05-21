#! /usr/bin/env python
import os, sys, getopt, multiprocessing, copy, math, itertools
import numpy as np
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom, Double
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, TH3F, THStack, TGraph, TGaxis, TH1D, TLorentzVector
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TMath
from ROOT import RooFit, RooRealVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooExtendPdf, RooAddPdf, RooProdPdf, RooNumConvPdf, RooFFTConvPdf

mumass = 0.10565837 # GeV

########## OPTIONS ##########

import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-b', '--batch',      action='store_true', default=False,           dest='batch'                    )
parser.add_option('-v', '--verbose',    action='store_true', default=False,           dest='verbose'                  )
parser.add_option('-d', '--datapath',   action='store',      default='data/2015/',    dest='datapath',  type='string' )
parser.add_option('-w', '--wspath',     action='store',      default='workspaces/',   dest='wspath',    type='string' )
(options, args) = parser.parse_args()
if options.batch: gROOT.SetBatch(True)

########## SETTINGS ##########

# Silent RooFit
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

# Remove stats and title from plots
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

# Input filenames
root_file_name  = 'DM_LHC15o_pDCAcutAndPcut025GeV_tree_ptFine_5Mix_new_var'
# root_file_name  = 'DM_LHC15o_pDCAcutAndPcut025GeV_tree_ptFine_5Mix'

# Output filenames
w_name          = 'ws_raw.root'
w_pT_name       = 'ws_pT_int.root'
w_C_name        = 'ws_C_int.root'
w_Full_name     = 'ws_Full.root'

########## CREATE WSs ##########

f       = TFile(options.datapath + root_file_name + '.root', "OLD")
w       = RooWorkspace("w",       root_file_name);
w_pT    = RooWorkspace("w_pT",    root_file_name);
w_C     = RooWorkspace("w_C",     root_file_name);
w_Full  = RooWorkspace("w_Full",  root_file_name);

C_bins    = ["020", "2050", "5090"]
pT_bins   = ["23"]
# pT_bins   = ["23", "34", "45", "56", "67", "78", "89", "90"]
pair_bins = ["PM"]
dm_bins   = [""]
#pT_bins   = ["34", "45", "56", "67", "78", "89", "90"]
# pair_bins = ["PM", "PP", "MM"]
# dm_bins   = ["", "Mix"]

# Chain all the trees together
tree = TChain()
for it in itertools.product(dm_bins,pair_bins,C_bins,pT_bins):
  tree_name     = "ftree{}{}{}{}".format(it[0],it[1],it[2],it[3])
  print('adding tree',tree_name)
  tree.Add(options.datapath + root_file_name + '.root/'+tree_name)
print('total n. entries =',tree.GetEntries())

# Define TLorentzVector lists for m1 and m2
list1 = []
list2 = []

# Instantiate a number of mass_shuffledos
c = TCanvas('c','c',800,600)
mass_standard = TH1F('mass_standard','mass_standard',1000,0,10)
pt_standard = TH1F('pt_standard','pt_standard',1000,0,20)
rapidity_standard = TH1F('rapidity_standard','rapidity_standard',1000,-5,0)
mass_vs_pt_standard = TH2F('mass_vs_pt_standard','mass_vs_pt_standard',100,0,10,100,0,20)
pt_vs_pt_standard = TH2F('pt_vs_pt_standard','pt_vs_pt_standard',100,0,10,100,0,10)

mass_tree = TH1F('mass_tree','mass_tree',1000,0,10)

mass_shuffled = TH1F('mass_shuffled','mass_shuffled',1000,0,10)
pt_shuffled = TH1F('pt_shuffled','pt_shuffled',1000,0,20)
rapidity_shuffled = TH1F('rapidity_shuffled','rapidity_shuffled',1000,-5,0)
mass_vs_pt_shuffled = TH2F('mass_vs_pt_shuffled','mass_vs_pt_shuffled',100,0,10,100,0,20)
pt_shuffled_tmp = TH1F('pt_shuffled_tmp','pt_shuffled_tmp',1000,0,20)
pt_vs_pt_shuffled = TH2F('pt_vs_pt_shuffled','pt_vs_pt_shuffled',100,0,10,100,0,10)

mass_shuffled_ptcut = TH1F('mass_shuffled_ptcut','mass_shuffled_ptcut',1000,0,10)
pt_shuffled_ptcut = TH1F('pt_shuffled_ptcut','pt_shuffled_ptcut',1000,0,20)
rapidity_shuffled_ptcut = TH1F('rapidity_shuffled_ptcut','rapidity_shuffled_ptcut',1000,-5,0)
mass_vs_pt_shuffled_ptcut = TH2F('mass_vs_pt_shuffled_ptcut','mass_vs_pt_shuffled_ptcut',100,0,10,100,0,20)
pt_shuffled_tmp_ptcut = TH1F('pt_shuffled_tmp_ptcut','pt_shuffled_tmp_ptcut',1000,0,20)
pt_vs_pt_shuffled_ptcut = TH2F('pt_vs_pt_shuffled_ptcut','pt_vs_pt_shuffled_ptcut',100,0,10,100,0,10)

# pt_ratio = TH1F('pt_ratio','pt_ratio',1000,0,20)

test = TH3F('test','test',1000,0,20,1000,0,20,1000,0,20)


# Loop over tree entries and fill histos and lists
print('filling standard lists')
for it in tree:
  l1 = TLorentzVector()
  l2 = TLorentzVector()
  l1.SetXYZM(it.tp1x,it.tp1y,it.tp1z,mumass)
  l2.SetXYZM(it.tp2x,it.tp2y,it.tp2z,mumass)
  # Fill standard lists
  list1.append(l1)
  list2.append(l2)  
  # Fill standard histos
  ll = l1+l2
  mass_tree.Fill          (it.tMass)
  mass_standard.Fill      (ll.M())
  pt_standard.Fill        (ll.Pt())
  rapidity_standard.Fill  (ll.Rapidity())
  mass_vs_pt_standard.Fill(ll.M(),ll.Pt())
  pt_vs_pt_standard.Fill(l1.Pt(),l2.Pt())
print('done')

print(len(list1),'elements found')

# First reshuffling - estimate of pT-dependent weight
print('shuffling both lists independently')
np.random.shuffle(list1)
np.random.shuffle(list2)

# Loop over shuffled lists - estimate of pT-dependent weight
for e1, e2 in zip(list1,list2):
  # Fill histo only if pT > 2 
  ll = e1+e2
  pt_shuffled_tmp.Fill(ll.Pt())
  # if ll.Pt() > 2:
  #   pt_shuffled_tmp.Fill(ll.Pt())
print('done')

# Evaluate and plot the pT-dependent weight
pt_ratio = pt_standard.Clone('pt_ratio')
pt_ratio.Divide(pt_shuffled_tmp)
pt_ratio.Draw("HIST")
c.SaveAs('pt_ratio.png')
c.SaveAs('pt_ratio.pdf')
c.SaveAs('pt_ratio.root')
input('\nHit return to continue...')

# Second reshuffling - apply pT-dependent weight and evaluate reshuffled distributions
print('shuffling both lists independently')
np.random.shuffle(list1)
np.random.shuffle(list2)

# Loop over shuffled lists - apply pT-dependent weight
for e1, e2 in zip(list1,list2):
  ll = e1+e2
  # w = pt_ratio.GetBinContent(pt_ratio.FindBin(ll.Pt()))
  w = 1.
  rapidity_shuffled.Fill  (ll.Rapidity(),   w)
  mass_shuffled.Fill      (ll.M(),          w)
  pt_shuffled.Fill        (ll.Pt(),         w)
  mass_vs_pt_shuffled.Fill(ll.M(),ll.Pt(),  w)
  pt_vs_pt_shuffled.Fill(e1.Pt(),e2.Pt())  
  if ll.Pt() > 2:
    rapidity_shuffled_ptcut.Fill  (ll.Rapidity())
    mass_shuffled_ptcut.Fill      (ll.M())
    pt_shuffled_ptcut.Fill        (ll.Pt())
    mass_vs_pt_shuffled_ptcut.Fill(ll.M(),ll.Pt())
    pt_vs_pt_shuffled_ptcut.Fill(e1.Pt(),e2.Pt())  
print('done')

# Draw and save plots
leg = TLegend(0.6, 0.9-0.05*3, 0.95, 0.9)
leg.SetBorderSize(0)
leg.SetFillStyle(0) #1001
leg.SetFillColor(0)
leg.AddEntry(pt_standard,'ordered TLorentzVector','f')
leg.AddEntry(pt_shuffled,'shuffled TLorentzVector','l')

pt_standard.SetLineColor(4)
pt_standard.SetFillColor(4)
pt_standard.GetXaxis().SetTitle('p_{T}(#mu#mu) [GeV]')
pt_standard.Draw()
pt_shuffled.SetLineColor(2)
pt_shuffled.Draw("E1, SAMES")
leg.Draw('SAME')
c.SaveAs('pt.pdf')
c.SaveAs('pt.png')
c.SaveAs('pt.root')

rapidity_standard.SetLineColor(4)
rapidity_standard.SetFillColor(4)
rapidity_standard.GetXaxis().SetTitle('y(#mu#mu)')
rapidity_standard.Draw()
rapidity_shuffled.SetLineColor(2)
rapidity_shuffled.Draw("E1, SAMES")
leg.Draw('SAME')
c.SaveAs('rapidity.pdf')
c.SaveAs('rapidity.png')
c.SaveAs('rapidity.root')

mass_vs_pt_standard.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
mass_vs_pt_standard.GetYaxis().SetTitle('p_{T}(#mu#mu) [GeV]')
mass_vs_pt_standard.Draw("COLZ")
mass_vs_pt_shuffled.Draw("BOX, SAMES")
leg.Draw('SAME')
c.SaveAs('mass_vs_pt.pdf')
c.SaveAs('mass_vs_pt.png')
c.SaveAs('mass_vs_pt.root')

pt_vs_pt_standard.GetXaxis().SetTitle('p_{T}(#mu_{1}) [GeV]')
pt_vs_pt_standard.GetYaxis().SetTitle('p_{T}(#mu_{2}) [GeV]')
pt_vs_pt_standard.Draw("COLZ")
pt_vs_pt_shuffled.Draw("BOX, SAMES")
leg.Draw('SAME')
c.SaveAs('pt_vs_pt.pdf')
c.SaveAs('pt_vs_pt.png')
c.SaveAs('pt_vs_pt.root')







pt_standard.SetLineColor(4)
pt_standard.SetFillColor(4)
pt_standard.GetXaxis().SetTitle('p_{T}(#mu#mu) [GeV]')
pt_standard.Draw()
pt_shuffled_ptcut.SetLineColor(2)
pt_shuffled_ptcut.Draw("E1, SAMES")
leg.Draw('SAME')
c.SaveAs('pt_ptcut.pdf')
c.SaveAs('pt_ptcut.png')
c.SaveAs('pt_ptcut.root')

rapidity_standard.SetLineColor(4)
rapidity_standard.SetFillColor(4)
rapidity_standard.GetXaxis().SetTitle('y(#mu#mu)')
rapidity_standard.Draw()
rapidity_shuffled_ptcut.SetLineColor(2)
rapidity_shuffled_ptcut.Draw("E1, SAMES")
leg.Draw('SAME')
c.SaveAs('rapidity_ptcut.pdf')
c.SaveAs('rapidity_ptcut.png')
c.SaveAs('rapidity_ptcut.root')

mass_vs_pt_standard.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
mass_vs_pt_standard.GetYaxis().SetTitle('p_{T}(#mu#mu) [GeV]')
mass_vs_pt_standard.Draw("COLZ")
mass_vs_pt_shuffled_ptcut.Draw("BOX, SAMES")
leg.Draw('SAME')
c.SaveAs('mass_vs_pt_ptcut.pdf')
c.SaveAs('mass_vs_pt_ptcut.png')
c.SaveAs('mass_vs_pt_ptcut.root')

pt_vs_pt_standard.GetXaxis().SetTitle('p_{T}(#mu_{1}) [GeV]')
pt_vs_pt_standard.GetYaxis().SetTitle('p_{T}(#mu_{2}) [GeV]')
pt_vs_pt_standard.Draw("COLZ")
pt_vs_pt_shuffled_ptcut.Draw("BOX, SAMES")
leg.Draw('SAME')
c.SaveAs('pt_vs_pt_ptcut.pdf')
c.SaveAs('pt_vs_pt_ptcut.png')
c.SaveAs('pt_vs_pt_ptcut.root')





leg.AddEntry(mass_tree,'tMass tree variable','l')

mass_standard.SetLineColor(4)
mass_standard.SetFillColor(4)
mass_standard.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
mass_standard.Draw()
mass_tree.SetLineColor(1)
mass_tree.Draw("E1, SAMES")
mass_shuffled.SetLineColor(2)
mass_shuffled.Draw("E1, SAMES")
leg.Draw('SAME')
c.SaveAs('mass.pdf')
c.SaveAs('mass.png')
c.SaveAs('mass.root')

mass_standard.SetLineColor(4)
mass_standard.SetFillColor(4)
mass_standard.GetXaxis().SetTitle('m(#mu#mu) [GeV]')
mass_standard.Draw()
mass_tree.SetLineColor(1)
mass_tree.Draw("E1, SAMES")
mass_shuffled_ptcut.SetLineColor(2)
mass_shuffled_ptcut.Draw("E1, SAMES")
leg.Draw('SAME')
c.SaveAs('mass_ptcut.pdf')
c.SaveAs('mass_ptcut.png')
c.SaveAs('mass_ptcut.root')

# input('\nHit return to continue...')

exit()




tMass = RooRealVar("tMass", "tMass", 1.5, 8, 'GeV')

data      = {}
data_pT   = {}
data_C    = {}
data_Full = {}
nevents   = {}

# All bins
for it in itertools.product(dm_bins,pair_bins,C_bins,pT_bins):
  tree_name     = "ftree{}{}{}{}".format(it[0],it[1],it[2],it[3])
  bin_name      = "{}_{}_C{}_pT{}".format(it[0],it[1],it[2],it[3])
  tree          = f.Get(tree_name)
  data[bin_name]    = RooDataSet('data%s'%bin_name,     'data%s'%bin_name,    RooArgSet(tMass), RooFit.Import(tree))
  nevents[bin_name] = RooRealVar('entries%s'%bin_name,  'entries%s'%bin_name, data[bin_name].numEntries())
  getattr(w, 'import')(data[bin_name])
  getattr(w, 'import')(nevents[bin_name])
  w.writeToFile('{}{}'.format(options.wspath,w_name))


# Integrate over pT
for it in pT_bins:
  pattern = '_pT{}'.format(it)
  for oldbin in [key for key in data if pattern in key]:
    newbin = oldbin.replace(pattern,'')
    if newbin in data_pT.keys(): 
      data_pT[newbin].append(data[oldbin])
    else:
      data_pT[newbin] = data[oldbin].Clone('data%s'%newbin)
for key in data_pT:
  nevents[key] = RooRealVar('entries%s'%key,  'entries%s'%key, data_pT[key].numEntries())
  getattr(w_pT, 'import')(data_pT[key])
  getattr(w_pT, 'import')(nevents[key])
  w_pT.writeToFile('{}{}'.format(options.wspath,w_pT_name))


# Integrate over Centrality
for it in C_bins:
  pattern = '_C{}'.format(it)
  for oldbin in [key for key in data if pattern in key]:
    newbin = oldbin.replace(pattern,'')
    if newbin in data_C.keys(): 
      data_C[newbin].append(data[oldbin])
    else:
      data_C[newbin] = data[oldbin].Clone('data%s'%newbin)
for key in data_C:
  nevents[key] = RooRealVar('entries%s'%key,  'entries%s'%key, data_C[key].numEntries())
  getattr(w_C, 'import')(data_C[key])
  getattr(w_C, 'import')(nevents[key])
  w_C.writeToFile('{}{}'.format(options.wspath,w_C_name))


# Integrate over pT and Centrality
for it in itertools.product(C_bins,pT_bins):
  pattern = '_C%s_pT%s'%(it[0],it[1])
  for oldbin in [key for key in data if pattern in key]:
    newbin = oldbin.replace(pattern,'')
    if newbin in data_Full.keys(): 
      data_Full[newbin].append(data[oldbin])
    else:
      data_Full[newbin] = data[oldbin].Clone('data%s'%newbin)
for key in data_Full:
  nevents[key] = RooRealVar('entries%s'%key,  'entries%s'%key, data_Full[key].sumEntries())
  getattr(w_Full, 'import')(data_Full[key])
  getattr(w_Full, 'import')(nevents[key])
  w_Full.writeToFile('{}{}'.format(options.wspath,w_Full_name))


w_Full.Print()
