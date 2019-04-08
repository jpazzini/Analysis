#! /usr/bin/env python
import os, sys, getopt, multiprocessing, copy, math, itertools
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom, Double
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis, TH1D
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TMath
from ROOT import RooFit, RooRealVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooExtendPdf, RooAddPdf, RooProdPdf, RooNumConvPdf, RooFFTConvPdf

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

# Input filenames
root_file_name  = 'DM_LHC15o_pDCAcutAndPcut025GeV_tree_ptFine_5Mix'

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
pT_bins   = ["23", "34", "45", "56", "67", "78", "89", "90"]
pair_bins = ["PM", "PP", "MM"]
dm_bins   = ["", "Mix"]

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
