#! /usr/bin/env python
import os, sys, getopt, multiprocessing, copy, math, itertools
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom, Double
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis, TH1D, TF1
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TMath, TGraphAsymmErrors
from ROOT import RooFit, RooRealVar, RooPolyVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooExtendPdf, RooAddPdf, RooProdPdf, RooNumConvPdf, RooFFTConvPdf

from drawUtils import *

print ('LOADING CUSTOM PDFs')
rc = gSystem.Load('../pdfs/HWWLVJRooPdfs_cxx.so') 
if not rc==0:
  print ('SOMETHING WENT WRONG')
  exit()

from ROOT import RooDoubleCrystalBall, RooExpNPdf, RooExpTailPdf, Roo2ExpPdf, RooErfExpPdf, RooErfExpNPdf, RooErfExpTailPdf, RooAlpha, RooAlpha4ErfExpPdf, RooExpN4Pdf, RooAlpha4ExpN4Pdf

RATIO = 3.75

########## OPTIONS ##########

import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-b', '--batch',      action='store_true', default=False,           dest='batch'                    )
parser.add_option('-v', '--verbose',    action='store_true', default=False,           dest='verbose'                  )
parser.add_option('-d', '--datapath',   action='store',      default='data/2015/',    dest='datapath',  type='string' )
parser.add_option('-w', '--wspath',     action='store',      default='../workspaces/',dest='wspath',    type='string' )
parser.add_option('-x', '--exclude_fit',action='store_true', default=False           ,dest='exclude_fit'              )
parser.add_option('-D', '--do_draw',    action='store_true', default=False,           dest='do_draw'                  )
parser.add_option('-n', '--n_bins',     action='store',      default=-1,              dest='n_bins',    type='int'    )
parser.add_option('-m', '--mass_min',   action='store',      default=2.5,             dest='mass_min',  type='float'  )
(options, args) = parser.parse_args()
if options.batch: gROOT.SetBatch(True)

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)

TGaxis.SetMaxDigits(2)
plot_out_folder = '../plots/'
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

def ComputeR(mixPM, mixPP, mixMM):
  return mixPM.numEntries() / (2 * math.sqrt(mixPP.numEntries() * mixMM.numEntries()))

def ComputeF(R, dataPP, dataMM, mixPM):
  return 2 * R * math.sqrt(dataPP.numEntries() * dataMM.numEntries()) / mixPM.numEntries()

def SubtractMixing(F, h_data, h_mixing):
  subtraction = h_data.Clone()
  h_mixing.Scale(F)
  subtraction.Add(h_mixing, -1)
  return subtraction

def MergeRooDataSet(dsvect):
  for ids in dsvect:
    dsvect[0].append(ids)

# ingest workspace (centrality-integrated / pT-binned)
file  = TFile("../workspaces/ws_C_int.root",'read')
w     = file.Get("w_C");
out_w = RooWorkspace("processed_w_C")
print ("    Workspace w corresponds to file:", w.GetTitle())
w.Print()

# define the mass-depentend features in ROOFIT
mass                  = w.var("tMass")
mass_range            = (options.mass_min, 8.)
min_mass_bin          = 0.025 # half of min mass resolution --> m_res ~ 60 MeV @ m = 2 GeV
n_bins                = options.n_bins if options.n_bins > 0 else round((8.-options.mass_min)/min_mass_bin)
n_bins_fullrange      = round((8.-1.5)/min_mass_bin)
binsmass_fullrange    = RooBinning(n_bins_fullrange, 1.5, 8.)
binsmass              = RooBinning(n_bins, mass_range[0], mass_range[1])
mass.setRange('rangemass',mass_range[0], mass_range[1])
mass.setRange('preJPsi',mass_range[0], 2.8)
mass.setRange('postJPsi',3.3,mass_range[1])

pT_bins   = ['34', '45', '56', '67', '78', '89', '90'] # pT_bins   = ['23', '34', '45', '56', '67', '78', '89', '90']

data_full   = w.data('data_PM_pT{}'.format(pT_bins[0]))
mix_full    = w.data('dataMix_PM_pT{}'.format(pT_bins[0]))

dataPP_full = w.data('data_PP_pT{}'.format(pT_bins[0]))
mixPP_full  = w.data('dataMix_PP_pT{}'.format(pT_bins[0]))

dataMM_full = w.data('data_MM_pT{}'.format(pT_bins[0]))
mixMM_full  = w.data('dataMix_MM_pT{}'.format(pT_bins[0]))

print(data_full.numEntries(),end = '')

for it in pT_bins[1:]:
  
  print('->',data_full.numEntries(),end = '')

  data_full.append( w.data('data_PM_pT{}'.format(it)) )
  mix_full.append( w.data('dataMix_PM_pT{}'.format(it)) )

  dataPP_full.append( w.data('data_PP_pT{}'.format(it)) )
  mixPP_full.append( w.data('dataMix_PP_pT{}'.format(it)) )

  dataMM_full.append( w.data('data_MM_pT{}'.format(it)) )
  mixMM_full.append( w.data('dataMix_MM_pT{}'.format(it)) )

print()
print()
print('TOTAL NUMBER OF EVENTS')
print('total number of entries in data   PM =',data_full.numEntries())
print('total number of entries in mixing PM =',mix_full.numEntries())
R = ComputeR(mix_full, mixPP_full, mixMM_full)
F = ComputeF(R, dataPP_full, dataMM_full, mix_full)
print('R',R)
print('F',F)
weight = RooRealVar("weight", "weight", -1, 1)
weight.setVal(1)
data_mix_full = RooDataSet("data_mix_full", "data_mix_full", RooArgSet(mass, weight), RooFit.Import(data_full), RooFit.WeightVar(weight))
weight.setVal(-F)
mixtmp_full   = RooDataSet("mixtmp_full", "mixtmp_full", RooArgSet(mass, weight), RooFit.Import(mix_full), RooFit.WeightVar(weight))
data_mix_full.append(mixtmp_full)
frame = mass.frame()
data_mix_full.plotOn(frame,
                     RooFit.Binning(binsmass_fullrange),
                     RooFit.DataError(RooAbsData.SumW2))
data_full.plotOn(frame,
                 RooFit.Binning(binsmass_fullrange),
                 RooFit.DataError(RooAbsData.Poisson),
                 RooFit.LineColor(2),
                 RooFit.MarkerColor(2))
mix_full.plotOn(frame,
                RooFit.Binning(binsmass_fullrange),
                RooFit.DataError(RooAbsData.Poisson),
                RooFit.LineColor(4),
                RooFit.MarkerColor(4))
# frame.Draw()
print()
input('hit return to continue...')
print()

print('EVALUATING R AND F ONLY ON THE RESTRICTED MASS RANGE USED FOR THE ANALYSIS')
print('number of entries after mass cut in data PM =',data_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])).numEntries())
print('number of entries after mass cut in mix  PM =',mix_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])).numEntries())
R_cut = ComputeR(mix_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])), mixPP_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])), mixMM_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])))
F_cut = ComputeF(R_cut, dataPP_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])), dataMM_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])), mix_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1])))
print('R',R_cut)
print('F',F_cut)
print('difference between R values = {:.3f}%'.format((R_cut-R)/R*100))
print('difference between F values = {:.3f}%'.format((F_cut-F)/F*100))

data_cut = data_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1]))
mix_cut  = mix_full.reduce('tMass > {} && tMass < {}'.format(mass_range[0],mass_range[1]))
weight_cut = RooRealVar("weight_cut", "weight_cut", -1, 1)
weight_cut.setVal(1)
data_mix_cut = RooDataSet("data_mix_cut", "data_mix_cut", RooArgSet(mass, weight_cut), RooFit.Import(data_cut), RooFit.WeightVar(weight_cut))
weight_cut.setVal(-F_cut)
mixtmp_cut   = RooDataSet("mixtmp_cut", "mixtmp_cut", RooArgSet(mass, weight_cut), RooFit.Import(mix_cut), RooFit.WeightVar(weight_cut))
data_mix_cut.append(mixtmp_cut)
data_mix_cut.plotOn(frame,
                    RooFit.Binning(binsmass),
                    RooFit.DataError(RooAbsData.SumW2),
                    RooFit.LineColor(3),
                    RooFit.MarkerColor(3))

data_cut.plotOn(frame,
                RooFit.Binning(binsmass),
                RooFit.DataError(RooAbsData.Poisson),
                RooFit.LineColor(5),
                RooFit.MarkerColor(5))
mix_cut.plotOn(frame,
               RooFit.Binning(binsmass),
               RooFit.DataError(RooAbsData.Poisson),
               RooFit.LineColor(6),
               RooFit.MarkerColor(6))
frame.Draw()

print()
input('hit return to continue...')
print()











"""
### FIT ON MIXING-SUBTRACTED DATA

### JPSI MODEL ###
meanJPsi    = RooRealVar('meanJPsi',  'meanJPsi',   3.09,  3.0,  3.15)
sigmaJPsi   = RooRealVar('sigmaJPsi', 'sigmaJPsi',  0.07, 0.03,  0.14)
alpha1JPsi  = RooRealVar('alpha1JPsi','alpha1JPsi', 0.8,   0.5,  1.5 )
n1JPsi      = RooRealVar('n1JPsi',    'n1JPsi',     4,       2,  7   )
alpha2JPsi  = RooRealVar('alpha2JPsi','alpha2JPsi', 0.8,   0.5,  2   )
n2JPsi      = RooRealVar('n2JPsi',    'n2JPsi',     4,     1,   15   )
# JPsiModel   = RooGaussian('JPsiModel','JPsiModel',mass,meanJPsi,sigmaJPsi)
JPsiModel   = RooCBShape('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi)
# JPsiModel   = RooDoubleCrystalBall('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

## PSIp MODEL ###
DmeanJPsiPsip  = RooRealVar   ('DmeanJPsiPsip', 'DmeanJPsiPsip',   0.590, 0.580, 0.600)
RsigmaJPsiPsip = RooRealVar   ('RsigmaJPsiPsip','RsigmaJPsiPsip',   1.02, 1.00,  1.10 )
meanPsip       = RooFormulaVar('meanPsip',      'meanPsip',   '@0+@1',  RooArgList(meanJPsi, DmeanJPsiPsip) )
sigmaPsip      = RooFormulaVar('sigmaPsip',     'sigmaPsip',  '@0*@1',  RooArgList(sigmaJPsi,RsigmaJPsiPsip))
# PsipModel      = RooGaussian('PsipModel','PsipModel',mass,meanPsip,sigmaPsip)
PsipModel      = RooCBShape('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi)
# PsipModel      = RooDoubleCrystalBall('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)
DmeanJPsiPsip.setConstant(True)
RsigmaJPsiPsip.setConstant(True)

### COMB BKG MODEL ###
## f(x) = exp(p_0 * x)
p0Comb        = RooRealVar    ('p0Comb',       'p0Comb',      -0.9, -2,   0)
BkgCombModel  = RooExponential('BkgCombModel', 'BkgCombModel', mass, p0Comb)

# ## f(x) = x**p_0
# p0Comb        = RooRealVar    ('p0Comb',       'p0Comb',      -3.5, -4.5, -2.5)
# BkgCombModel  = RooGenericPdf ('BkgCombModel', 'BkgCombModel', 'pow(@0,@1)', RooArgList(mass,p0Comb))

# ## f(x) = exp(p_0 * x) * (1 + Erf( (x - p_1) / p_2) )/2
# p0Comb        = RooRealVar    ('p0Comb',       'p0Comb',      -0.9, -2,    0  )
# p1Comb        = RooRealVar    ('p1Comb',       'p1Comb',       0.8,  0.5,  1. )
# p2Comb        = RooRealVar    ('p2Comb',       'p2Comb',       0.3,  0.1,  0.8)
# BkgCombModel  = RooErfExpPdf  ('BkgCombModel', 'BkgCombModel', mass, p0Comb, p1Comb, p2Comb)

# ## f(x) = exp(p_0 * x) + p_2 * exp(p_1 * x)
# p0Comb        = RooRealVar    ('p0Comb',       'p0Comb',      -2,   -3,   -0.5  )
# p1Comb        = RooRealVar    ('p1Comb',       'p1Comb',      -0.5, -1,    0.01 )
# p2Comb        = RooRealVar    ('p2Comb',       'p2Comb',       0.1,  0.05, 0.5  )
# BkgCombModel  = Roo2ExpPdf    ('BkgCombModel', 'BkgCombModel', mass, p0Comb, p1Comb, p2Comb)

# ## f(x) = exp(p_0 * x + p_1 / x)
# p0Comb        = RooRealVar    ('p0Comb',       'p0Comb',      -1,   -2,    0    )
# p1Comb        = RooRealVar    ('p1Comb',       'p1Comb',      10,    0,   20    )
# BkgCombModel  = RooExpNPdf    ('BkgCombModel', 'BkgCombModel', mass, p0Comb, p1Comb)

# ## f(x) = exp(-x / (p_0 + p_1 * x))
# p0Comb        = RooRealVar    ('p0Comb',       'p0Comb',     3e-2,  1e-3,   1e-1 )
# p1Comb        = RooRealVar    ('p1Comb',       'p1Comb',     4e-2,  1e-3,   1e-1 )
# BkgCombModel  = RooExpTailPdf ('BkgCombModel', 'BkgCombModel', mass, p0Comb, p1Comb)

### RELATIVE FRACTIONS ###
fJPsi     = RooRealVar('fJPsi', 'fJPsi', 4e-1, 1e-2, 1e0 )
fPsip     = RooRealVar('fPsip', 'fPsip', 1e-5, 1e-6, 1e-1)

baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel),RooArgList(fJPsi))
# baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel,PsipModel),RooArgList(fJPsi,fPsip))
fitRes    = baseModel.fitTo(data_mix_cut,
                            # RooFit.Extended(True),
                            RooFit.SumW2Error(True),
                            RooFit.Strategy(2),
                            RooFit.Minimizer('Minuit'),
                            RooFit.Range('rangemass'),
                            RooFit.Save(1),
                            RooFit.NumCPU(12),
                            RooFit.PrintLevel(1 if options.verbose else -1) )

print ('=== OVERALL FIT RESULTS ===')
fitRes.Print()
print ('===========================')

c = TCanvas('c','c',800,800)
c.Divide(1,2)
setTopPad(c.GetPad(1),RATIO)
setBotPad(c.GetPad(2),RATIO)
# draw top pad
c.cd(1)
frame = mass.frame()
frame.GetYaxis().SetTitle('Events / x GeV')
frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
setHistStyle(frame,1.1)
data_obs = data_mix_cut.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(1), RooFit.Range('rangemass'),RooFit.Name('tot'))
chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
res  = frame.pullHist()
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(1), RooFit.Components('BkgCombModel,JPsiModel'), RooFit.Range('rangemass'),RooFit.Name('JPsi'))
baseModel.plotOn(frame, RooFit.LineColor(4), RooFit.LineStyle(2), RooFit.Components('BkgCombModel'), RooFit.Range('rangemass'),RooFit.Name('comb'))
frame.GetYaxis().SetTitleOffset(1.1)
frame.Draw()
drawALICE('','Internal')
# draw bottom pad
c.cd(2)
frame_res = mass.frame()
frame_res.GetYaxis().SetTitle('Events / x GeV')
frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
setHistStyle(frame_res,1.1)
setBotStyle(frame_res,RATIO)
frame_res.addPlotable(res,'P')
frame_res.GetYaxis().SetRangeUser(-5,5)
frame_res.GetYaxis().SetTitle('(data-fit)/#sigma_{data}')
frame_res.GetYaxis().CenterTitle()
frame_res.GetYaxis().SetTitleOffset(0.4)
frame_res.Draw()
line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)
c.SaveAs(plot_out_folder+'simpleFit.pdf')
c.SaveAs(plot_out_folder+'simpleFit.png')
c.SaveAs(plot_out_folder+'simpleFit.root')

c.SetLogy()
c.SaveAs(plot_out_folder+'simpleFit_LOG.pdf')
c.SaveAs(plot_out_folder+'simpleFit_LOG.png')

print()
input('hit return to continue...')
print()
"""

















"""
### FIT ON MIXING

# ## COMB BKG MODEL ###
# ## f(x) = exp(p_0 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.9,    -1.5,      -0.5   )
# BkgCombModel  = RooExponential('BkgCombModel','BkgCombModel', mass, p0Comb)

# ## f(x) = x**p_0
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.5,  -5,    0   )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'pow(@0,@1)', RooArgList(mass, p0Comb))

# ## f(x) = exp(p_0 * x) + p_2 * exp(p_1 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.85, -0.92,  -0.5 )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -1.00, -5,      0 )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',    0.1,   0.,   0.95)
# BkgCombModel  = Roo2ExpPdf('BkgCombModel',    'BkgCombModel', mass, p0Comb, p1Comb, p2Comb)

# ## f(x) = x**p_0 * exp(p_0 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -6,  -100,    100  )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -1,    -5,      0  )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'pow(@0,@1)*exp(pow(@0,1)*@2)', RooArgList(mass, p0Comb, p1Comb))

# ## f(x) = x**(p_0 + p_1 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb', -3.7e-01, -2,  0)
# p1Comb        = RooRealVar('p1Comb',          'p1Comb', -3.3e-01, -2,  0)
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'pow(@0, @1 + @2 * @0 )', RooArgList(mass, p0Comb, p1Comb))

# ## f(x) = exp(p_0 + p_1 * x + p_2 / x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   14.2,     7,    28    )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -1.05,   -2.1,  -0.05 )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',   -2.68,   -5.3,  -1.2  )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'exp(@1 + pow(@0,1)*@2 + pow(@0,-1)*@3)', RooArgList(mass, p0Comb, p1Comb, p2Comb))

# ## f(x) = exp(p_0 + p_1 * x + p_2 * x**2)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   12.4,       5,      25      )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -0.65,     -1.4,    -0.05   )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',   -2.68e-2,  -5.3e-2, -1.2e-2 )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3)', RooArgList(mass, p0Comb, p1Comb, p2Comb))

## f(x) = exp(p_0 + p_1 * x + p_2 * x**2 + p_3 / x)
p0Comb        = RooRealVar('p0Comb',          'p0Comb',   12.7,       5,        25    )
p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -7.3e-1,  -15.e-1,  -3.2e-1 )
p2Comb        = RooRealVar('p2Comb',          'p2Comb',   -2.2e-2,  -5e-2,    -1e-2   )
p3Comb        = RooRealVar('p3Comb',          'p3Comb',   -5.1e-1,  -10e-1,   -1.5e-1 )
BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3 + pow(@0,-1)*@4)', RooArgList(mass, p0Comb, p1Comb, p2Comb, p3Comb))

baseModel = BkgCombModel
fitRes    = baseModel.fitTo(mix_cut,
                            # RooFit.Extended(True),
                            # RooFit.SumW2Error(True),
                            RooFit.Strategy(2),
                            RooFit.Minimizer('Minuit'),
                            RooFit.Range('rangemass'),
                            RooFit.Save(1),
                            RooFit.NumCPU(12),
                            RooFit.PrintLevel(1 if options.verbose else -1) )

print ('=== OVERALL FIT RESULTS ===')
fitRes.Print()
print ('===========================')

c = TCanvas('c','c',800,800)
c.Divide(1,2)
setTopPad(c.GetPad(1),RATIO)
setBotPad(c.GetPad(2),RATIO)
# draw top pad
c.cd(1)
frame = mass.frame()
frame.GetYaxis().SetTitle('Events / x GeV')
frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
setHistStyle(frame,1.1)
mix_obs = mix_cut.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mix_obs'))
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(1), RooFit.Range('rangemass'),RooFit.Name('tot'))
chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
res  = frame.pullHist()
frame.GetYaxis().SetTitleOffset(1.1)
frame.Draw()
drawALICE('','Internal')
# draw bottom pad
c.cd(2)
frame_res = mass.frame()
frame_res.GetYaxis().SetTitle('Events / x GeV')
frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
setHistStyle(frame_res,1.1)
setBotStyle(frame_res,RATIO)
frame_res.addPlotable(res,'P')
frame_res.GetYaxis().SetRangeUser(-5,5)
frame_res.GetYaxis().SetTitle('(data-fit)/#sigma_{data}')
frame_res.GetYaxis().CenterTitle()
frame_res.GetYaxis().SetTitleOffset(0.4)
frame_res.Draw()
line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)
c.SaveAs(plot_out_folder+'simpleFit.pdf')
c.SaveAs(plot_out_folder+'simpleFit.png')
c.SaveAs(plot_out_folder+'simpleFit.root')

c.SetLogy()
c.SaveAs(plot_out_folder+'simpleFit_LOG.pdf')
c.SaveAs(plot_out_folder+'simpleFit_LOG.png')

print()
input('hit return to continue...')
print()
"""



























### FIT ON DATA W/O MIXING SUBTRACTION

### JPSI MODEL ###
meanJPsi    = RooRealVar('meanJPsi',  'meanJPsi',   3.09,  3.05,  3.15)
sigmaJPsi   = RooRealVar('sigmaJPsi', 'sigmaJPsi',  0.07, 0.03,  0.14)
alpha1JPsi  = RooRealVar('alpha1JPsi','alpha1JPsi', 0.8,   0.5,  1.5 )
n1JPsi      = RooRealVar('n1JPsi',    'n1JPsi',     4,       2,  7   )
alpha2JPsi  = RooRealVar('alpha2JPsi','alpha2JPsi', 0.8,   0.5,  2   )
n2JPsi      = RooRealVar('n2JPsi',    'n2JPsi',     4,     1,   15   )
# JPsiModel   = RooGaussian('JPsiModel','JPsiModel',mass,meanJPsi,sigmaJPsi)
JPsiModel   = RooCBShape('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi)
# JPsiModel   = RooDoubleCrystalBall('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)
# sigmaJPsi.setConstant(True)
# n1JPsi.setConstant(True)
# alpha1JPsi.setConstant(True)


## PSIp MODEL ###
DmeanJPsiPsip  = RooRealVar   ('DmeanJPsiPsip', 'DmeanJPsiPsip',   0.590, 0.580, 0.600)
RsigmaJPsiPsip = RooRealVar   ('RsigmaJPsiPsip','RsigmaJPsiPsip',   1.02, 1.00,  1.10 )
meanPsip       = RooFormulaVar('meanPsip',      'meanPsip',   '@0+@1',  RooArgList(meanJPsi, DmeanJPsiPsip) )
sigmaPsip      = RooFormulaVar('sigmaPsip',     'sigmaPsip',  '@0*@1',  RooArgList(sigmaJPsi,RsigmaJPsiPsip))
# PsipModel      = RooGaussian('PsipModel','PsipModel',mass,meanPsip,sigmaPsip)
PsipModel      = RooCBShape('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi)
# PsipModel      = RooDoubleCrystalBall('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)
DmeanJPsiPsip.setConstant(True)
RsigmaJPsiPsip.setConstant(True)

### COMB BKG MODEL ###
# ## f(x) = exp(p_0 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.9,    -1.5,      -0.5   )
# BkgCombModel  = RooExponential('BkgCombModel','BkgCombModel', mass, p0Comb)

# ## f(x) = x**p_0
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.5,  -5,    0 )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'pow(@0,@1)', RooArgList(mass, p0Comb))

# ## f(x) = exp(p_0 * x) + p_2 * exp(p_1 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.85, -0.92,  -0.5 )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -1.00, -5,      0 )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',    0.1,   0.,   0.95)
# BkgCombModel  = Roo2ExpPdf('BkgCombModel',    'BkgCombModel', mass, p0Comb, p1Comb, p2Comb)

# ## f(x) = x**p_0 * exp(p_0 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   -0.4,  -2,    0  )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -0.9,  -3,    0  )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'pow(@0,@1)*exp(pow(@0,1)*@2)', RooArgList(mass, p0Comb, p1Comb))

# ## f(x) = x**(p_0 + p_1 * x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb', -3.7e-01, -2,  0)
# p1Comb        = RooRealVar('p1Comb',          'p1Comb', -3.3e-01, -2,  0)
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'pow(@0, @1 + @2 * @0 )', RooArgList(mass, p0Comb, p1Comb))

# ## f(x) = exp(p_0 + p_1 * x + p_2 / x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   14.2,     7,    28    )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -1.05,   -2.1,  -0.05 )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',   -2.68,   -5.3,  -1.2  )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'exp(@1 + pow(@0,1)*@2 + pow(@0,-1)*@3)', RooArgList(mass, p0Comb, p1Comb, p2Comb))

# ## f(x) = exp(p_0 + p_1 * x + p_2 * x**2)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   12.4,       5,      25      )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -0.65,     -1.4,    -0.05   )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',   -2.68e-2,  -5.3e-2, -1.2e-2 )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3)', RooArgList(mass, p0Comb, p1Comb, p2Comb))

# ## f(x) = exp(p_0 + p_1 * x + p_2 * x**2 + p_3 / x)
# p0Comb        = RooRealVar('p0Comb',          'p0Comb',   12.7,       5,        25    )
# p1Comb        = RooRealVar('p1Comb',          'p1Comb',   -7.3e-1,  -15.e-1,  -3.2e-1 )
# p2Comb        = RooRealVar('p2Comb',          'p2Comb',   -1.5e-2,  -1e-1,    -1e-3   )
# p3Comb        = RooRealVar('p3Comb',          'p3Comb',   -5.1e-1,  -10e-1,   -1.5e-1 )
# BkgCombModel  = RooGenericPdf('BkgCombModel', 'BkgCombModel', 'exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3 + pow(@0,-1)*@4)', RooArgList(mass, p0Comb, p1Comb, p2Comb, p3Comb))

### RELATIVE FRACTIONS ###
fJPsi     = RooRealVar('fJPsi', 'fJPsi', 4e-1, 1e-2, 1e0 )
fPsip     = RooRealVar('fPsip', 'fPsip', 1e-5, 1e-6, 1e-1)

baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel),RooArgList(fJPsi))
# baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel,PsipModel),RooArgList(fJPsi,fPsip))
fitRes    = baseModel.fitTo(data_cut,
                            # RooFit.Extended(True),
                            # RooFit.SumW2Error(True),
                            RooFit.Strategy(2),
                            RooFit.Minimizer('Minuit'),
                            RooFit.Range('rangemass'),
                            RooFit.Save(1),
                            RooFit.NumCPU(12),
                            RooFit.PrintLevel(1 if options.verbose else -1) )

print ('=== OVERALL FIT RESULTS ===')
fitRes.Print()
print ('===========================')

c = TCanvas('c','c',800,800)
c.Divide(1,2)
setTopPad(c.GetPad(1),RATIO)
setBotPad(c.GetPad(2),RATIO)
# draw top pad
c.cd(1)
frame = mass.frame()
frame.GetYaxis().SetTitle('Events / x GeV')
frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
setHistStyle(frame,1.1)
data_nomix_obs = data_cut.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data_nomix_obs'))
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(1), RooFit.Range('rangemass'),RooFit.Name('tot'))
chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
res  = frame.pullHist()
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(1), RooFit.Components('BkgCombModel,JPsiModel'), RooFit.Range('rangemass'),RooFit.Name('JPsi'))
baseModel.plotOn(frame, RooFit.LineColor(4), RooFit.LineStyle(2), RooFit.Components('BkgCombModel'), RooFit.Range('rangemass'),RooFit.Name('comb'))
frame.GetYaxis().SetTitleOffset(1.1)
frame.Draw()
drawALICE('','Internal')
# draw bottom pad
c.cd(2)
frame_res = mass.frame()
frame_res.GetYaxis().SetTitle('Events / x GeV')
frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
setHistStyle(frame_res,1.1)
setBotStyle(frame_res,RATIO)
frame_res.addPlotable(res,'P')
frame_res.GetYaxis().SetRangeUser(-5,5)
frame_res.GetYaxis().SetTitle('(data-fit)/#sigma_{data}')
frame_res.GetYaxis().CenterTitle()
frame_res.GetYaxis().SetTitleOffset(0.4)
frame_res.Draw()
line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)
c.SaveAs(plot_out_folder+'simpleFit.pdf')
c.SaveAs(plot_out_folder+'simpleFit.png')
c.SaveAs(plot_out_folder+'simpleFit.root')

c.SetLogy()
c.SaveAs(plot_out_folder+'simpleFit_LOG.pdf')
c.SaveAs(plot_out_folder+'simpleFit_LOG.png')

print()
input('hit return to continue...')
print()










































































































#########################################################





# ### OVERALL FIT 

# # file = TFile('../workspaces/processed_w_Full.root','read')
# # w = file.Get('processed_w_Full')
# # w.Print()

# # mass = w.var('tMass')
# # data = w.data('data_mix')


# def fitter(pTbin=None):

#   file = TFile('../workspaces/ws_C_int.root' if pTbin else '../workspaces/ws_Full.root','read')
#   w = file.Get('w_C' if pTbin else 'w_Full')
#   w.Print()

#   suffix = '_pT%s'%pTbin if pTbin else ''

#   mass = w.var('tMass')
#   data = w.data('data_PM%s'%suffix)
#   mix  = w.data('dataMix_PM%s'%suffix)
#   Nmix = w.var('entriesMix_PM%s'%suffix).getVal()
#   Ndata= w.var('entries_PM%s'%suffix).getVal()


#   dataPP  = w.data('data_PP%s'%suffix)
#   mixPP   = w.data('dataMix_PP%s'%suffix)
#   NdataPP = w.var('entries_PP%s'%suffix).getVal()
#   NmixPP  = w.var('entriesMix_PP%s'%suffix).getVal()

#   dataMM  = w.data('data_MM%s'%suffix)
#   mixMM   = w.data('dataMix_MM%s'%suffix)
#   NdataMM = w.var('entries_MM%s'%suffix).getVal()
#   NmixMM  = w.var('entriesMix_MM%s'%suffix).getVal()

#   mass.setRange('preJPsi',mass_range[0],2.7)
#   mass.setRange('postJPsi',3.3,mass_range[1])

#   mass.setRange('rangemass',mass_range[0], mass_range[1])
#   mass.setRange('fitrangemass',fit_mass_range[0], fit_mass_range[1])

#   R   = 0.5 * Nmix / math.sqrt(NmixPP*NmixMM)
#   dRx = 0.5 / math.sqrt(NmixPP*NmixMM)
#   dRy = - 0.25 * Nmix * NmixMM / pow(NmixPP*NmixMM,3/2)
#   dRz = - 0.25 * Nmix * NmixPP / pow(NmixPP*NmixMM,3/2)
#   dR  = math.sqrt ( dRx**2*Nmix + dRy**2*NmixPP + dRz**2*NmixMM ) 

#   F     = 2 * R * math.sqrt(NdataPP*NdataMM) / Nmix
#   F_bis = math.sqrt( (NdataPP*NdataMM) / (NmixPP*NmixMM)) 
#   F_ter = (NdataPP + NdataMM) / (NmixPP + NmixMM)

#   print("R    = %f +/- %f"%(R,dR))
#   print("F    = %f"%F)
#   print("F_bis= %f"%F_bis)
#   print("F_ter= %f"%F_ter)

#   c = TCanvas('c','c',800,600)
#   c.cd()
#   frame = mass.frame()
#   dataPP.plotOn(frame, RooFit.Normalization(NdataPP, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataPP'),RooFit.LineColor(2),RooFit.MarkerColor(2))
#   dataMM.plotOn(frame, RooFit.Normalization(NdataMM, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
#   frame.Draw()
#   c.SaveAs('c_data.root')

#   cMix = TCanvas('cMix','cMix',800,600)
#   cMix.cd()
#   frameMix = mass.frame()
#   mixPP.plotOn(frameMix, RooFit.Normalization(NmixPP, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixPP'),RooFit.LineColor(2),RooFit.MarkerColor(2))
#   mixMM.plotOn(frameMix, RooFit.Normalization(NmixMM, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
#   frameMix.Draw()
#   cMix.SaveAs('c_mix.root')

#   cPP = TCanvas('cPP','cPP',800,600)
#   cPP.cd()
#   framePP = mass.frame()
#   mixPP.plotOn(framePP, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixPP'),RooFit.LineColor(2),RooFit.MarkerColor(2))
#   dataPP.plotOn(framePP, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataPP'),RooFit.LineColor(4),RooFit.MarkerColor(4))
#   framePP.Draw()
#   cPP.SaveAs('c_PP.root')

#   cMM = TCanvas('cMM','cMM',800,600)
#   cMM.cd()
#   frameMM = mass.frame()
#   mixMM.plotOn(frameMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixMM'),RooFit.LineColor(2),RooFit.MarkerColor(2))
#   dataMM.plotOn(frameMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
#   frameMM.Draw()
#   cMM.SaveAs('c_MM.root')

#   drawComparisonPlot(mass,dataPP,NdataPP,dataMM,NdataMM,'data_PP-MM%s'%suffix)
#   drawComparisonPlot(mass,mixPP, NmixPP, mixMM, NmixMM, 'mix_PP-MM%s'%suffix)
#   drawComparisonPlot(mass,dataPP,NdataPP,mixPP,NmixPP,'data-mix_PP%s'%suffix)
#   drawComparisonPlot(mass,dataMM,NdataMM,mixMM,NmixMM,'data-mix_MM%s'%suffix)

#   drawDoubleComparisonPlot(mass,dataPP,NdataPP,dataMM,NdataMM,mixPP, NmixPP, mixMM, NmixMM,'doubleComparison%s'%suffix)

#   if options.verbose: input("SINGLE-DATASET COMPARISONS DONE\nMERGING PP AND MM TOGETHER")


#   dataPPMM = RooDataSet(dataPP,'dataPPMM')
#   dataPPMM.append(dataMM)
#   mixPPMM  = RooDataSet(mixPP, 'mixPPMM')
#   mixPPMM.append(mixMM)

#   NdataPPMM= dataPPMM.numEntries()
#   NmixPPMM = mixPPMM.numEntries()


#   cPPMM = TCanvas('cPPMM','cPPMM',800,600)
#   cPPMM.cd()
#   framePPMM = mass.frame()
#   mixPPMM.plotOn(framePPMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixPPMM'),RooFit.LineColor(2),RooFit.MarkerColor(2))
#   dataPPMM.plotOn(framePPMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataPPMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
#   framePPMM.Draw()
#   cPPMM.SaveAs('c_PPMM.root')

#   drawComparisonPlot(mass,dataPPMM,NdataPPMM,mixPPMM,NmixPPMM,'data-mix_PPMM%s'%suffix)


#   histMixPPMM  = mixPPMM.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histMixPPMM').ProjectionX('histMixPPMM')
#   histDataPPMM = dataPPMM.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histDataPPMM').ProjectionX('histDataPPMM')
#   histDataPPMM.Sumw2()
#   histDataPPMM.Divide(histMixPPMM)

#   ca = TCanvas('ca','ca',800,600)
#   ca.cd()
#   histDataPPMM.Draw()

#   histMixPM    = mix.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histMixPM').ProjectionX('histMixPM')
#   histMixPM.Sumw2()
#   histMixPM.Multiply(histDataPPMM)

#   histFMixPM    = mix.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histFMixPM').ProjectionX('histFMixPM')
#   histFMixPM.Sumw2()
#   histFMixPM.Scale(F)

#   cb = TCanvas('cb','cb',800,600)
#   cb.cd()
#   histMixPM.SetLineColor(4)
#   histMixPM.Draw()
#   histFMixPM.SetLineColor(2)
#   histFMixPM.Draw('SAME')

#   histDataPM    = data.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histDataPM').ProjectionX('histDataPM')
#   histDataPM.Sumw2()
  
#   cc = TCanvas('cc','cc',800,600)
#   cc.cd()



#   histMixPM.SetLineColor(4)
#   histMixPM.Draw('HIST')

#   histFMixPM.SetLineColor(2)
#   histFMixPM.Draw('HIST, SAME')

#   histDataPM.SetLineColor(1)
#   histDataPM.SetMarkerColor(1)
#   histDataPM.SetMarkerStyle(21)
#   histDataPM.Draw('E1, SAME')

#   cd = TCanvas('cd','cd',800,600)
#   cd.cd()

#   histDataPM_histMixPM  = histDataPM.Clone('histDataPM_histMixPM')
#   histDataPM_histMixPM.Add(histMixPM,-1)
#   histDataPM_histMixPM.Sumw2()
#   histDataPM_histMixPM.SetMarkerColor(2)
#   histDataPM_histMixPM.SetLineColor(2)

#   histDataPM_histFMixPM = histDataPM.Clone('histDataPM_histFMixPM')
#   histDataPM_histFMixPM.Add(histFMixPM,-1)
#   histDataPM_histFMixPM.Sumw2()
#   histDataPM_histFMixPM.SetMarkerColor(4)
#   histDataPM_histFMixPM.SetLineColor(4)  

#   histDataPM_histMixPM.Draw('E1')
#   histDataPM_histFMixPM.Draw('E1, SAME')

#   ce = TCanvas('ce','ce',800,600)
#   ce.cd()

#   subhist = [histDataPM_histMixPM,histDataPM_histFMixPM]
#   minhist = 0 if subhist[0].GetMinimum() < subhist[1].GetMinimum() else 1
#   minbin  = subhist[minhist].GetMinimumBin()
#   minval  = subhist[minhist].GetMinimum() - subhist[minhist].GetBinError(minbin)

#   flatval = abs(minval)*1.2
#   flatfun = TF1("flatfun","[0]",mass_range[0], mass_range[1])
#   flatfun.SetParameter(0, flatval)

#   histDataPM_histMixPM.Add(flatfun)
#   histDataPM_histFMixPM.Add(flatfun)

#   histDataPM_histMixPM.Draw('E1')
#   histDataPM_histFMixPM.Draw('E1, SAME')

#   input('---PAUSA---')

























#   datasub = RooDataSet("datasub", "datasub", RooArgSet(mass), RooFit.Import(data))
#   weight  = RooRealVar("weight", "weight", -1, 1)
#   weight.setVal(-F)
#   mixtmp  = RooDataSet("mixtmp", "mixtmp", RooArgSet(mass, weight), RooFit.Import(mix), RooFit.WeightVar(weight))
#   datasub.append(mixtmp)
#   Ndatasub= datasub.numEntries()
#   print('weight',weight.getVal(),'data', data.numEntries(), 'mix',mixtmp.numEntries(), 'datasub', datasub.numEntries())
#   getattr(w, 'import')(datasub)
#   datasub.Print()


#   input('STOP')

#   ### JPSI MODEL ###
#   meanJPsi    = RooRealVar('meanJPsi',  'meanJPsi',   3.1,  3.0,  3.2 )
#   sigmaJPsi   = RooRealVar('sigmaJPsi', 'sigmaJPsi',  5e-3, 1e-3, 5e-2)
#   alpha1JPsi  = RooRealVar('alpha1JPsi','alpha1JPsi', 0.8,   0.5,  1.5 )
#   n1JPsi      = RooRealVar('n1JPsi',    'n1JPsi',     4,    1.5,  7   )
#   alpha2JPsi  = RooRealVar('alpha2JPsi','alpha2JPsi', 0.8,   0.5,  2   )
#   n2JPsi      = RooRealVar('n2JPsi',    'n2JPsi',     4,    1,    15  )

#   # JPsiModel  = RooGaussian('JPsiModel','JPsiModel',mass,meanJPsi,sigmaJPsi)
#   JPsiModel = RooCBShape('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi)
#   # JPsiModel   = RooDoubleCrystalBall('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

#   ## PSIp MODEL ###
#   DmeanJPsiPsip  = RooRealVar('DmeanJPsiPsip','DmeanJPsiPsip',0.590,0.580,0.600)
#   RsigmaJPsiPsip = RooRealVar('RsigmaJPsiPsip','RsigmaJPsiPsip',1.02,1.00,1.10)
#   meanPsip  = RooFormulaVar('meanPsip','meanPsip','@0+@1',RooArgList(meanJPsi,DmeanJPsiPsip))
#   sigmaPsip = RooFormulaVar('sigmaPsip','sigmaPsip','@0*@1',RooArgList(sigmaJPsi,RsigmaJPsiPsip))

#   # PsipModel  = RooGaussian('PsipModel','PsipModel',mass,meanPsip,sigmaPsip)
#   PsipModel   = RooCBShape('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi)
#   # PsipModel   = RooDoubleCrystalBall('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

#   # ### COMB BKG EXPO MODEL ###
#   # constComb = RooRealVar('constComb','constComb',-1,-2,0)
#   # BkgCombModel   = RooExponential('BkgCombModel','BkgCombModel',mass,constComb)

#   # # ### COMB BKG EXPO MODEL ###
#   # # constComb  = RooRealVar('constComb','constComb',-2,-2,-0.5)
#   # # offsetComb = RooRealVar('offsetComb','offsetComb',1.1,0.8,1.4)
#   # # widthComb  = RooRealVar('widthComb','widthComb',0.5,0.3,0.8)
#   # # BkgCombModel   = RooErfExpPdf('BkgCombModel','BkgCombModel',mass,constComb,offsetComb,widthComb)

#   # # ### COMB BKG 2 EXPO MODEL ###
#   # # const1Comb = RooRealVar('const1Comb','const1Comb',-0.94,-1.3,-0.5)
#   # # const2Comb = RooRealVar('const2Comb','const2Comb',-1.79,-2.5,-1.3)
#   # # fComb = RooRealVar('fComb','fComb',0.23,0.15,0.45)
#   # # BkgCombModel   = Roo2ExpPdf('BkgCombModel','BkgCombModel',mass,const1Comb,const2Comb,fComb)

#   # # ### COMB BKG EXPN MODEL ###
#   # # numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
#   # # denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e-3,   1e-6,   1e2)
#   # # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)

#   # # ### COMB BKG EXPTAIL MODEL ###
#   # # numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
#   # # denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e3,   1e0,   1e5)
#   # # # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)
#   # # BkgCombModel= RooExpTailPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)



#   # ============ GOOD FUNCTIONS ON dataMMsub --> TO BE TESTED ON PM DATA SUBTRACTED

#   ### QUITE BAD

#   # ### f(x) = exp(p_0 * x)
#   # constComb = RooRealVar('constComb','constComb',-1,-2,0)
#   # BkgCombModel   = RooExponential('BkgCombModel','BkgCombModel',mass,constComb)

#   # ### f(x) = x**p_0
#   # powComb = RooRealVar('powComb','powComb',-6,-100,100)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','pow(@0,@1)',RooArgList(mass,powComb))

#   ### NOT SO BAD

#   # ### f(x) = exp(p_0 * x) + f * exp(p_1 * x)
#   # const1Comb = RooRealVar('const1Comb','const1Comb',-0.85,-0.92,-0.5)
#   # const2Comb = RooRealVar('const2Comb','const2Comb',-1.00,-5, -0.9)
#   # fComb = RooRealVar('fComb','fComb',0.1,0.05,0.95)
#   # BkgCombModel   = Roo2ExpPdf('BkgCombModel','BkgCombModel',mass,const1Comb,const2Comb,fComb)

#   # ### f(x) = x**p_0 * exp(p_0 * x)
#   # constComb = RooRealVar('constComb','constComb',-1,-5,0)
#   # powComb = RooRealVar('powComb','powComb',-6,-100,100)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','pow(@0,@1)*exp(pow(@0,1)*@2)',RooArgList(mass,powComb,constComb))

#   # ### f(x) = x**(p_0 + p_1 * x)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',-3.6719e-01,-2,0)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-3.2799e-01,-2,0)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','pow(@0, @1 + @2 * @0 )',RooArgList(mass,p_0_Comb,p_1_Comb))

#   # ### f(x) = exp(p_0 + p_1 * x + p_2 / x)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',14.2,7,28)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-1.05,-2.1,-0.05)
#   # p_m1_Comb = RooRealVar('p_m1_Comb','p_m1_Comb',-2.68,-5.3,-1.2)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + pow(@0,1)*@2 + pow(@0,-1)*@3)',RooArgList(mass,p_0_Comb,p_1_Comb,p_m1_Comb))

#   # ### f(x) = exp(p_0 + p_1 * x + p_2 * x**2)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',12.4,5,25)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-0.65,-1.4,-0.05)
#   # p_2_Comb = RooRealVar('p_2_Comb','p_2_Comb',-2.68e-2,-5.3e-2,-1.2e-2)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3)',RooArgList(mass,p_0_Comb,p_1_Comb,p_2_Comb))

#   ### f(x) = exp(p_0 + p_1 * x + p_2 * x**2 + p_3 / x)
#   p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',12.7,5,25)
#   p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-7.3e-1,-15.e-1,-3.2e-1)
#   p_2_Comb = RooRealVar('p_2_Comb','p_2_Comb',-2.2e-2,-5e-2,-1e-2)
#   p_m1_Comb = RooRealVar('p_m1_Comb','p_m1_Comb',-5.1e-1,-10e-1,-1.5e-1)
#   BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3 + pow(@0,-1)*@4)',RooArgList(mass,p_0_Comb,p_1_Comb,p_2_Comb,p_m1_Comb))



#   fJPsi     = RooRealVar('fJPsi','fJPsi',1e-2,1e-3,1e0)
#   fPsip     = RooRealVar('fPsip','fPsip',1e-5,1e-6,1e-2)

#   DmeanJPsiPsip.setConstant(True)
#   RsigmaJPsiPsip.setConstant(True)

#   fitRes    = BkgCombModel.fitTo(datasub,
#                             RooFit.SumW2Error(True),
#                             RooFit.Strategy(2),
#                             RooFit.Minimizer('Minuit'),
#                             # RooFit.Range('rangemass'),
#                             RooFit.Range('preJPsi,postJPsi'),
#                             RooFit.Save(1),
#                             RooFit.NumCPU(12),
#                             RooFit.PrintLevel(1))

#   baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel),RooArgList(fJPsi))
#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt,PsipModelExt),RooArgList(NBkgComb,NJPsi,NPsip))
#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel,PsipModel),RooArgList(fJPsi,fPsip))
#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt),RooArgList(NBkgComb,NJPsi))
#   fitRes    = baseModel.fitTo(datasub,
#                             RooFit.SumW2Error(True),
#                             RooFit.Strategy(2),
#                             RooFit.Minimizer('Minuit'),
#                             RooFit.Range('rangemass'),
#                             # RooFit.Range('preJPsi,postJPsi'),
#                             RooFit.Save(1),
#                             RooFit.NumCPU(12),
#                             RooFit.PrintLevel(1))




#   print ('===   BASE FIT    ===')
#   fitRes.Print()
#   print ('===========================')

#   drawPlot(mass,datasub,Ndatasub,BkgCombModel,fitRes,'testFit%s'%suffix)



#   input('---PAUSA---')























#   dataMMsub = RooDataSet("dataMMsub", "dataMMsub", RooArgSet(mass), RooFit.Import(dataMM))
#   weight  = RooRealVar("weight", "weight", -1, 1)
#   weight.setVal(-F)
#   mixMMtmp  = RooDataSet("mixMMtmp", "mixMMtmp", RooArgSet(mass, weight), RooFit.Import(mixMM), RooFit.WeightVar(weight))
#   dataMMsub.append(mixMMtmp)
#   NdataMMsub= dataMMsub.numEntries()


#   can = TCanvas('can','can',800,600)
#   can.cd()
#   framecan = mass.frame()
#   dataMMsub.plotOn(framecan, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'))
#   framecan.Draw()
#   can.SaveAs('can.root')


#   ### JPSI MODEL ###
#   meanJPsi    = RooRealVar('meanJPsi',  'meanJPsi',   3.1,  3.0,  3.2 )
#   sigmaJPsi   = RooRealVar('sigmaJPsi', 'sigmaJPsi',  5e-3, 1e-3, 5e-2)
#   alpha1JPsi  = RooRealVar('alpha1JPsi','alpha1JPsi', 0.8,   0.5,  1.5 )
#   n1JPsi      = RooRealVar('n1JPsi',    'n1JPsi',     4,    1.5,  7   )
#   alpha2JPsi  = RooRealVar('alpha2JPsi','alpha2JPsi', 0.8,   0.5,  2   )
#   n2JPsi      = RooRealVar('n2JPsi',    'n2JPsi',     4,    1,    15  )

#   # JPsiModel  = RooGaussian('JPsiModel','JPsiModel',mass,meanJPsi,sigmaJPsi)
#   JPsiModel = RooCBShape('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi)
#   # JPsiModel   = RooDoubleCrystalBall('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

#   ## PSIp MODEL ###
#   DmeanJPsiPsip  = RooRealVar('DmeanJPsiPsip','DmeanJPsiPsip',0.590,0.580,0.600)
#   RsigmaJPsiPsip = RooRealVar('RsigmaJPsiPsip','RsigmaJPsiPsip',1.02,1.00,1.10)
#   meanPsip  = RooFormulaVar('meanPsip','meanPsip','@0+@1',RooArgList(meanJPsi,DmeanJPsiPsip))
#   sigmaPsip = RooFormulaVar('sigmaPsip','sigmaPsip','@0*@1',RooArgList(sigmaJPsi,RsigmaJPsiPsip))

#   # PsipModel  = RooGaussian('PsipModel','PsipModel',mass,meanPsip,sigmaPsip)
#   PsipModel   = RooCBShape('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi)
#   # PsipModel   = RooDoubleCrystalBall('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

#   # ### COMB BKG EXPO MODEL ###
#   # constComb = RooRealVar('constComb','constComb',-1,-2,0)
#   # BkgCombModel   = RooExponential('BkgCombModel','BkgCombModel',mass,constComb)

#   # ### COMB BKG EXPO MODEL ###
#   # constComb  = RooRealVar('constComb','constComb',-2,-2,-0.5)
#   # offsetComb = RooRealVar('offsetComb','offsetComb',1.1,0.8,1.4)
#   # widthComb  = RooRealVar('widthComb','widthComb',0.5,0.3,0.8)
#   # BkgCombModel   = RooErfExpPdf('BkgCombModel','BkgCombModel',mass,constComb,offsetComb,widthComb)

#   # ### COMB BKG 2 EXPO MODEL ###
#   # const1Comb = RooRealVar('const1Comb','const1Comb',-0.85,-0.92,-0.5)
#   # const2Comb = RooRealVar('const2Comb','const2Comb',-1.00,-5, -0.9)
#   # fComb = RooRealVar('fComb','fComb',0.1,0.05,0.95)
#   # BkgCombModel   = Roo2ExpPdf('BkgCombModel','BkgCombModel',mass,const1Comb,const2Comb,fComb)

#   # ### COMB BKG EXPN MODEL ###
#   # numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
#   # denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e-1,   1e-3,   1e2)
#   # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)






#   # ============ GOOD FUNCTIONS ON dataMMsub --> TO BE TESTED ON PM DATA SUBTRACTED

#   ### QUITE BAD

#   # ### f(x) = exp(p_0 * x)
#   # constComb = RooRealVar('constComb','constComb',-1,-2,0)
#   # BkgCombModel   = RooExponential('BkgCombModel','BkgCombModel',mass,constComb)

#   # ### f(x) = x**p_0
#   # powComb = RooRealVar('powComb','powComb',-6,-100,100)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','pow(@0,@1)',RooArgList(mass,powComb))

#   ### NOT SO BAD

#   # ### f(x) = exp(p_0 * x) + f * exp(p_1 * x)
#   # const1Comb = RooRealVar('const1Comb','const1Comb',-0.85,-0.92,-0.5)
#   # const2Comb = RooRealVar('const2Comb','const2Comb',-1.00,-5, -0.9)
#   # fComb = RooRealVar('fComb','fComb',0.1,0.05,0.95)
#   # BkgCombModel   = Roo2ExpPdf('BkgCombModel','BkgCombModel',mass,const1Comb,const2Comb,fComb)

#   # ### f(x) = x**p_0 * exp(p_0 * x)
#   # constComb = RooRealVar('constComb','constComb',-1,-5,0)
#   # powComb = RooRealVar('powComb','powComb',-6,-100,100)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','pow(@0,@1)*exp(pow(@0,1)*@2)',RooArgList(mass,powComb,constComb))

#   # ### f(x) = x**(p_0 + p_1 * x)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',-3.6719e-01,-2,0)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-3.2799e-01,-2,0)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','pow(@0, @1 + @2 * @0 )',RooArgList(mass,p_0_Comb,p_1_Comb))

#   # ### f(x) = exp(p_0 + p_1 * x + p_2 / x)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',14.2,7,28)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-1.05,-2.1,-0.05)
#   # p_m1_Comb = RooRealVar('p_m1_Comb','p_m1_Comb',-2.68,-5.3,-1.2)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + pow(@0,1)*@2 + pow(@0,-1)*@3)',RooArgList(mass,p_0_Comb,p_1_Comb,p_m1_Comb))

#   # ### f(x) = exp(p_0 + p_1 * x + p_2 * x**2)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',12.4,5,25)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-0.65,-1.4,-0.05)
#   # p_2_Comb = RooRealVar('p_2_Comb','p_2_Comb',-2.68e-2,-5.3e-2,-1.2e-2)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3)',RooArgList(mass,p_0_Comb,p_1_Comb,p_2_Comb))

#   # ### f(x) = exp(p_0 + p_1 * x + p_2 * x**2 + p_3 / x)
#   # p_0_Comb = RooRealVar('p_0_Comb','p_0_Comb',12.7,5,25)
#   # p_1_Comb = RooRealVar('p_1_Comb','p_1_Comb',-7.3e-1,-15.e-1,-3.2e-1)
#   # p_2_Comb = RooRealVar('p_2_Comb','p_2_Comb',-2.2e-2,-5e-2,-1e-2)
#   # p_m1_Comb = RooRealVar('p_m1_Comb','p_m1_Comb',-5.1e-1,-10e-1,-1.5e-1)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + pow(@0,1)*@2 + pow(@0,2)*@3 + pow(@0,-1)*@4)',RooArgList(mass,p_0_Comb,p_1_Comb,p_2_Comb,p_m1_Comb))

















#   # ### COMB BKG EXPTAIL MODEL ###
#   # numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
#   # denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e3,   1e0,   1e5)
#   # # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)
#   # BkgCombModel= RooExpTailPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)


#   # NJPsi     = RooRealVar('NJPsi','NJPsi',1e3,1e1,1e6)
#   # JPsiModelExt = RooExtendPdf('JPsiModelExt','JPsiModelExt',JPsiModel,NJPsi)
#   # NPsip     = RooRealVar('NPsip','NPsip',1e1,1e-1,1e4)
#   # PsipModelExt = RooExtendPdf('PsipModelExt','PsipModelExt',PsipModel,NPsip)
#   # NBkgComb = RooRealVar('NBkgComb','NBkgComb',  1e6,   1e1,   1e8)
#   # BkgCombModelExt = RooExtendPdf('BkgCombModelExt','BkgCombModelExt',BkgCombModel,NBkgComb)

#   fJPsi     = RooRealVar('fJPsi','fJPsi',1e-2,1e-3,1e0)
#   fPsip     = RooRealVar('fPsip','fPsip',1e-5,1e-6,1e-2)

#   DmeanJPsiPsip.setConstant(True)
#   RsigmaJPsiPsip.setConstant(True)

#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel),RooArgList(fJPsi))
#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt,PsipModelExt),RooArgList(NBkgComb,NJPsi,NPsip))
#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel,PsipModel),RooArgList(fJPsi,fPsip))
#   # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt),RooArgList(NBkgComb,NJPsi))

#   # cMixPM  = RooRealVar('cMixPM','cMixPM',   1.2e+1, 0.6e+1, 2.4e+1)
#   # n1MixPM = RooRealVar('n1MixPM','n1MixPM',-4.0e-1,-8.0e-1,-2.0e-1)
#   # n2MixPM = RooRealVar('n2MixPM','n2MixPM',-1.8e-1,-3.6e-1,-0.9e-1)
#   # n3MixPM = RooRealVar('n3MixPM','n3MixPM', 1.0e-2, 0.5e-2, 2.0e-2)
#   # # d1MixPM = RooRealVar('d1MixPM','d1MixPM', 2.4e0 , 1.2e0 , 4.8e0 )
#   # # d2MixPM = RooRealVar('d2MixPM','d2MixPM',-1.6e-2,-3.2e-2,-1.6e-2)
#   # BkgCombModel   = RooGenericPdf('BkgCombModel','BkgCombModel','exp(@1 + @0*@2 + pow(@0,2)*@3 + pow(@0,3)*@4)',RooArgList(mass,cMixPM,n1MixPM,n2MixPM,n3MixPM))

#   fitRes    = BkgCombModel.fitTo(dataMMsub,
#                             RooFit.SumW2Error(True),
#                             RooFit.Strategy(2),
#                             RooFit.Minimizer('Minuit'),
#                             RooFit.Range('rangemass'),
#                             RooFit.Save(1),
#                             RooFit.NumCPU(12),
#                             RooFit.PrintLevel(1) )




#   print ('===   BASE FIT    ===')
#   fitRes.Print()
#   print ('===========================')

#   drawPlot(mass,dataMMsub,NdataMMsub,BkgCombModel,fitRes,'testFit%s'%suffix)



#   input('---PAUSA---')












#   if options.verbose: input("COMBINED-DATASET COMPARISONS DONE\nMOVING TO FITTING DATA")


#   # ### BKG EXPO MODEL ###
#   # constDataPPMM = RooRealVar('constDataPPMM','constDataPPMM',-1,-2,0)
#   # DataPPMMModel   = RooExponential('DataPPMMModel','DataPPMMModel',mass,constDataPPMM)
#   ### BKG ERF+EXPO MODEL ###
#   constDataPPMM  = RooRealVar('constDataPPMM','constDataPPMM',-1,-3,0)
#   offsetDataPPMM = RooRealVar('offsetDataPPMM','offsetDataPPMM',0.3,0.,4.)
#   widthDataPPMM  = RooRealVar('widthDataPPMM','widthDataPPMM',0.5,0.3,5.)
#   DataPPMMModel  = RooErfExpPdf('DataPPMMModel','DataPPMMModel',mass,constDataPPMM,offsetDataPPMM,widthDataPPMM)
#   # ### BKG 2 EXPO MODEL ###
#   # const1DataPPMM = RooRealVar('const1DataPPMM','const1DataPPMM',-2,-2,2)
#   # const2DataPPMM = RooRealVar('const2DataPPMM','const2DataPPMM',-0.5,-1,-1e-2)
#   # fDataPPMM = RooRealVar('fDataPPMM','fDataPPMM',0.1,0,0.5)
#   # DataPPMMModel   = Roo2ExpPdf('DataPPMMModel','DataPPMMModel',mass,const1DataPPMM,const2DataPPMM,fDataPPMM)
#   # ### BKG EXPN MODEL ###
#   # numBkgDataPPMM  = RooRealVar('numBkgDataPPMM','numBkgDataPPMM',  -1,  -2, 0)
#   # denBkgDataPPMM  = RooRealVar('denBkgDataPPMM','denBkgDataPPMM',  1e-3,   1e-6,   1e2)
#   # DataPPMMModel   = RooExpNPdf('DataPPMMModel','DataPPMMModel',mass,numBkgDataPPMM,denBkgDataPPMM)
  

#   # meangaus  = RooRealVar('meangaus','meangaus',-7,-10,0)
#   # sigmagaus = RooRealVar('sigmagaus','sigmagaus',1,0,2)
#   # gaus      = RooGaussian('gaus', 'gaus', mass, meangaus, sigmagaus)
#   # # exp1const = RooRealVar('exp1const','exp1const',-0.9,-1.,-0.5)
#   # # exp1      = RooExponential('exp1','exp1',mass,exp1const)

#   cheb0     = RooRealVar('cheb0','cheb0',-0.2, -5, 5)
#   cheb1     = RooRealVar('cheb1','cheb1',0.2,  -5, 5)
#   cheb2     = RooRealVar('cheb2','cheb2',-0.2, -1, 1)    
#   cheb      = RooChebychev('cheb','cheb',mass,RooArgList(cheb0,cheb1,cheb2))
#   exp2const = RooRealVar('exp2const','exp2const',-1.2,-1.5,-0.8)
#   exp2      = RooExponential('exp2','exp2',mass,exp2const)
#   exp3const = RooRealVar('exp3const','exp3const',-1.5,-2,-1.1)
#   exp3      = RooExponential('exp3','exp3',mass,exp3const)
#   f1        = RooRealVar('f1','f1',0.8,0.1,0.8)
#   f2        = RooRealVar('f2','f2',0.05,0.001,0.3)
#   # DataPPMMModel = RooAddPdf('DataPPMMModel','Test fit model',RooArgList(cheb,exp2,exp3),RooArgList(f1,f2))


#   # n1DataPPMM = RooRealVar('n1DataPPMM','n1DataPPMM',-9.2e-1,-10.e-1,-8.0e-1)
#   # d1DataPPMM = RooRealVar('d1DataPPMM','d1DataPPMM', 2.4e0 , 1.2e0 , 4.8e0 )
#   # n2DataPPMM = RooRealVar('n2DataPPMM','n2DataPPMM',-2.6e0 ,-5.2e0 , 1.3e0 )
#   # d2DataPPMM = RooRealVar('d2DataPPMM','d2DataPPMM',-1.6e-2,-3.2e-2,-1.6e-2)
#   # # DataPPMMModel   = RooExpN4Pdf('DataPPMMModel','DataPPMMModel',mass,n1DataPPMM,d1DataPPMM,n2DataPPMM,d2DataPPMM)


#   cDataPPMM  = RooRealVar('cDataPPMM','cDataPPMM',   1.2e+1, 0.6e+1, 2.4e+1)
#   n1DataPPMM = RooRealVar('n1DataPPMM','n1DataPPMM',-4.0e-1,-8.0e-1,-2.0e-1)
#   n2DataPPMM = RooRealVar('n2DataPPMM','n2DataPPMM',-1.8e-1,-3.6e-1,-0.9e-1)
#   n3DataPPMM = RooRealVar('n3DataPPMM','n3DataPPMM', 1.0e-2, 0.5e-2, 2.0e-2)
#   # d1DataPPMM = RooRealVar('d1DataPPMM','d1DataPPMM', 2.4e0 , 1.2e0 , 4.8e0 )
#   # d2DataPPMM = RooRealVar('d2DataPPMM','d2DataPPMM',-1.6e-2,-3.2e-2,-1.6e-2)
#   DataPPMMModel   = RooGenericPdf('DataPPMMModel','DataPPMMModel','exp(@1 + @0*@2 + pow(@0,2)*@3 + pow(@0,3)*@4)',RooArgList(mass,cDataPPMM,n1DataPPMM,n2DataPPMM,n3DataPPMM))

#   DataPPMMFitRes    = DataPPMMModel.fitTo(dataPPMM,
#                               RooFit.SumW2Error(True),
#                               RooFit.Strategy(2),
#                               RooFit.Minimizer('Minuit'),
#                               RooFit.Range('fitrangemass'),
#                               RooFit.Save(1),
#                               RooFit.NumCPU(12),
#                               RooFit.PrintLevel(1 if options.verbose else -1) )

#   print ('===   DataPPMM 1D FIT    ===')
#   DataPPMMFitRes.Print()
#   print ('===========================')

#   drawPlot(mass,dataPPMM,NdataPPMM,DataPPMMModel,DataPPMMFitRes,'dataPPMM%s'%suffix)




#   if options.verbose: input("FITTING DATA DONE\nMOVING TO FITTING MIXING")




#   # ### BKG EXPO MODEL ###
#   # constMixPPMM = RooRealVar('constMixPPMM','constMixPPMM',-1,-2,0)
#   # MixPPMMModel   = RooExponential('MixPPMMModel','MixPPMMModel',mass,constMixPPMM)
#   ### BKG ERF+EXPO MODEL ###
#   constMixPPMM  = RooRealVar('constMixPPMM','constMixPPMM',-1,-3,0)
#   offsetMixPPMM = RooRealVar('offsetMixPPMM','offsetMixPPMM',0.3,0.,4.)
#   widthMixPPMM  = RooRealVar('widthMixPPMM','widthMixPPMM',0.5,0.3,5.)
#   # MixPPMMModel  = RooErfExpPdf('MixPPMMModel','MixPPMMModel',mass,constMixPPMM,offsetMixPPMM,widthMixPPMM)
#   # ### BKG 2 EXPO MODEL ###
#   # const1MixPPMM = RooRealVar('const1MixPPMM','const1MixPPMM',-2,-2,2)
#   # const2MixPPMM = RooRealVar('const2MixPPMM','const2MixPPMM',-0.5,-1,-1e-2)
#   # fMixPPMM = RooRealVar('fMixPPMM','fMixPPMM',0.1,0,0.5)
#   # MixPPMMModel   = Roo2ExpPdf('MixPPMMModel','MixPPMMModel',mass,const1MixPPMM,const2MixPPMM,fMixPPMM)
#   # ### BKG EXPN MODEL ###
#   # numBkgMixPPMM  = RooRealVar('numBkgMixPPMM','numBkgMixPPMM',  -1,  -2, 0)
#   # denBkgMixPPMM  = RooRealVar('denBkgMixPPMM','denBkgMixPPMM',  1e-3,   1e-6,   1e2)
#   # MixPPMMModel   = RooExpNPdf('MixPPMMModel','MixPPMMModel',mass,numBkgMixPPMM,denBkgMixPPMM)
  
#   cMixPPMM  = RooRealVar('cMixPPMM','cMixPPMM',   1.2e+1, 0.6e+1, 2.4e+1)
#   n1MixPPMM = RooRealVar('n1MixPPMM','n1MixPPMM',-4.0e-1,-8.0e-1,-2.0e-1)
#   n2MixPPMM = RooRealVar('n2MixPPMM','n2MixPPMM',-1.8e-1,-3.6e-1,-0.9e-1)
#   n3MixPPMM = RooRealVar('n3MixPPMM','n3MixPPMM', 1.0e-2, 0.5e-2, 2.0e-2)
#   # d1MixPPMM = RooRealVar('d1MixPPMM','d1MixPPMM', 2.4e0 , 1.2e0 , 4.8e0 )
#   # d2MixPPMM = RooRealVar('d2MixPPMM','d2MixPPMM',-1.6e-2,-3.2e-2,-1.6e-2)
#   MixPPMMModel   = RooGenericPdf('MixPPMMModel','MixPPMMModel','exp(@1 + @0*@2 + pow(@0,2)*@3 + pow(@0,3)*@4)',RooArgList(mass,cMixPPMM,n1MixPPMM,n2MixPPMM,n3MixPPMM))


#   MixPPMMFitRes    = MixPPMMModel.fitTo(mixPPMM,
#                               RooFit.SumW2Error(True),
#                               RooFit.Strategy(2),
#                               RooFit.Minimizer('Minuit'),
#                               RooFit.Range('fitrangemass'),
#                               RooFit.Save(1),
#                               RooFit.NumCPU(12),
#                               RooFit.PrintLevel(1 if options.verbose else -1) )

#   print ('===   MixPPMM 1D FIT    ===')
#   MixPPMMFitRes.Print()
#   print ('===========================')

#   drawPlot(mass,mixPPMM,NmixPPMM,MixPPMMModel,MixPPMMFitRes,'mixPPMM%s'%suffix)




#   if options.verbose: input("FITTING MIXING DONE\nMOVING TO FITTING MIXING IN +- EVENTS")




#   # ### BKG EXPO MODEL ###
#   # constMixPM = RooRealVar('constMixPM','constMixPM',-1,-2,0)
#   # MixPMModel   = RooExponential('MixPMModel','MixPMModel',mass,constMixPM)
#   ### BKG ERF+EXPO MODEL ###
#   constMixPM  = RooRealVar('constMixPM','constMixPM',-1,-3,0)
#   offsetMixPM = RooRealVar('offsetMixPM','offsetMixPM',0.3,0.,4.)
#   widthMixPM  = RooRealVar('widthMixPM','widthMixPM',0.5,0.3,5.)
#   # MixPMModel  = RooErfExpPdf('MixPMModel','MixPMModel',mass,constMixPM,offsetMixPM,widthMixPM)
#   # ### BKG 2 EXPO MODEL ###
#   # const1MixPM = RooRealVar('const1MixPM','const1MixPM',-2,-2,2)
#   # const2MixPM = RooRealVar('const2MixPM','const2MixPM',-0.5,-1,-1e-2)
#   # fMixPM = RooRealVar('fMixPM','fMixPM',0.1,0,0.5)
#   # MixPMModel   = Roo2ExpPdf('MixPMModel','MixPMModel',mass,const1MixPM,const2MixPM,fMixPM)
#   # ### BKG EXPN MODEL ###
#   # numBkgMixPM  = RooRealVar('numBkgMixPM','numBkgMixPM',  -1,  -2, 0)
#   # denBkgMixPM  = RooRealVar('denBkgMixPM','denBkgMixPM',  1e-3,   1e-6,   1e2)
#   # MixPMModel   = RooExpNPdf('MixPMModel','MixPMModel',mass,numBkgMixPM,denBkgMixPM)
  
#   cMixPM  = RooRealVar('cMixPM','cMixPM',   1.2e+1, 0.6e+1, 2.4e+1)
#   n1MixPM = RooRealVar('n1MixPM','n1MixPM',-4.0e-1,-8.0e-1,-2.0e-1)
#   n2MixPM = RooRealVar('n2MixPM','n2MixPM',-1.8e-1,-3.6e-1,-0.9e-1)
#   n3MixPM = RooRealVar('n3MixPM','n3MixPM', 1.0e-2, 0.5e-2, 2.0e-2)
#   # d1MixPM = RooRealVar('d1MixPM','d1MixPM', 2.4e0 , 1.2e0 , 4.8e0 )
#   # d2MixPM = RooRealVar('d2MixPM','d2MixPM',-1.6e-2,-3.2e-2,-1.6e-2)
#   MixPMModel   = RooGenericPdf('MixPMModel','MixPMModel','exp(@1 + @0*@2 + pow(@0,2)*@3 + pow(@0,3)*@4)',RooArgList(mass,cMixPM,n1MixPM,n2MixPM,n3MixPM))

#   MixPMFitRes    = MixPMModel.fitTo(mix,
#                               RooFit.SumW2Error(True),
#                               RooFit.Strategy(2),
#                               RooFit.Minimizer('Minuit'),
#                               RooFit.Range('fitrangemass'),
#                               RooFit.Save(1),
#                               RooFit.NumCPU(12),
#                               RooFit.PrintLevel(1 if options.verbose else -1) )

#   print ('===   mixPM 1D FIT    ===')
#   MixPMFitRes.Print()
#   print ('===========================')

#   drawPlot(mass,mix,Nmix,MixPMModel,MixPMFitRes,'mixPM%s'%suffix)





#   if options.verbose: input("FITTING MIXING IN +- EVENTS DONE\nATTEMPTING THE SIMULTANEOUS FIT")





#   variables = RooArgSet(mass,'variables')
#   # define the three categories (not ranges) to be used in the fit
#   reg = RooCategory('reg', 'reg')
#   reg.defineType('dataPPMM')
#   reg.defineType('mixPPMM')
#   reg.defineType('mixPM')
#   # combine all the datasets with the corresponding categories together in one big sample
#   setMulti = RooDataSet('setMulti', 'setMulti', variables, RooFit.Index(reg), RooFit.Import('dataPPMM', dataPPMM), RooFit.Import('mixPPMM', mixPPMM), RooFit.Import('mixPM', mix))

#   # define the simultaneous object by linking to the big-dataset the three PDF that will be used in each of the regions
#   simObj = RooSimultaneous('simObj', 'simultaneous pdf', RooArgList(DataPPMMModel, MixPPMMModel, MixPMModel), reg)

#   # define the alpha-ratio object and combinatory bg in SR
#   # alpha = RooAlpha4ErfExpPdf('alpha', '#alpha function (ErfExp)', mass, constDataPPMM, offsetDataPPMM, widthDataPPMM, constMixPPMM, offsetMixPPMM, widthMixPPMM)
#   alpha   = RooGenericPdf('alpha','#alpha function','exp((@1-@5) + @0*(@2-@6) + pow(@0,2)*(@3-@7) + pow(@0,3)*(@4-@8))',RooArgList(mass,cDataPPMM,n1DataPPMM,n2DataPPMM,n3DataPPMM,cMixPPMM,n1MixPPMM,n2MixPPMM,n3MixPPMM))


#   DataCombPM = RooProdPdf('DataCombPM',  'Comb background estimation with alpha', alpha, MixPMModel)
     
#   frSim = simObj.fitTo(setMulti, RooFit.Range('fitrangemass'), RooFit.Strategy(2), RooFit.Minimizer('Minuit'), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if options.verbose else -1), RooFit.NumCPU(12))



#   print ('===   SIMULTANEOUS FIT  ===')
#   frSim.Print()
#   print ('===========================')


#   drawPlotSimultaneous(mass,setMulti,NmixPPMM,simObj,frSim,reg,'mixPPMM', 'Sim-mixPPMM%s'%suffix)
#   drawPlotSimultaneous(mass,setMulti,NmixPPMM,simObj,frSim,reg,'dataPPMM','Sim-dataPPMM%s'%suffix)
#   drawPlotSimultaneous(mass,setMulti,NmixPPMM,simObj,frSim,reg,'mixPM',   'Sim-mixPM%s'%suffix)
#   drawAlphaPlot(mass, alpha, MixPMModel, DataCombPM, frSim)
#   drawPlot(mass,data,Ndata,DataCombPM,frSim,'dataPMComb%s'%suffix,True)


#   file.Close()


# def drawPlot(mass,data,norm,model,fitRes,title,useResiduals=False):

#   c = TCanvas('c','c',800,800)
#   c.Divide(1,2)
#   setTopPad(c.GetPad(1),RATIO)
#   setBotPad(c.GetPad(2),RATIO)

#   # draw top pad
#   c.cd(1)
#   frame = mass.frame()
#   frame.GetYaxis().SetTitle('Events / x GeV')
#   frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame,1.1)
#   # data_obs = data.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
#   # model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('fitrangemass'),RooFit.Name('tot'))

# # data_obs = data.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))

  
#   data_obs = data.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
#   model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('preJPsi,postJPsi'),RooFit.Name('tot'))
#   chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
#   res  = frame.residHist() if useResiduals else frame.pullHist()
#   # model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.LineColor(2), RooFit.Range('fitrangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.LineColor(0),RooFit.FillColor(1),RooFit.FillStyle(3005))
#   # model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('fitrangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.LineColor(0),RooFit.FillColor(1),RooFit.FillStyle(3005))
#   # model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('fitrangemass'),RooFit.Name('tot'))
#   model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('preJPsi,postJPsi'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.LineColor(0),RooFit.FillColor(1),RooFit.FillStyle(3005))
#   model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('preJPsi,postJPsi'),RooFit.Name('tot'))
#   frame.GetYaxis().SetTitleOffset(1.1)
#   frame.SetMaximum(frame.GetMaximum()*1.1)
#   frame.Draw()
#   drawALICE('','Internal')
#   leg = TLegend(0.5, 0.65, 0.95, 0.9)
#   leg.SetBorderSize(0)
#   leg.SetFillStyle(0) #1001
#   leg.SetFillColor(0)
#   leg.AddEntry('data_obs', data.GetName().replace('_',' '), 'LP')
#   leg.AddEntry('tot',      'fit',                            'L' )
#   leg.AddEntry('Unc',      'fit uncertainty',                'F' )
#   leg.Draw()

#   # draw bottom pad
#   c.cd(2)
#   frame_res = mass.frame()
#   frame_res.GetYaxis().SetTitle('Events / x GeV')
#   frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame_res,1.1)
#   setBotStyle(frame_res,RATIO)
#   frame_res.addPlotable(res,'P')
#   frame_res.GetYaxis().SetRangeUser(-100 if useResiduals else -5, 100 if useResiduals else 5) 
#   frame_res.GetYaxis().SetTitle('data-fit' if useResiduals else '(data-fit)/#sigma_{data}')
#   frame_res.GetYaxis().CenterTitle()
#   frame_res.GetYaxis().SetTitleOffset(0.4)
#   frame_res.Draw()
#   line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
#   drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)

#   # save canvas in LIN scale
#   c.SaveAs(plot_out_folder+title+'.pdf')
#   c.SaveAs(plot_out_folder+title+'.png')
#   c.SaveAs(plot_out_folder+title+'.root')

#   # update and save canvas in LOG scale
#   c.GetPad(1).SetLogy()
#   frame.SetMinimum(max(frame.GetMinimum(),1.))
#   frame.SetMaximum(frame.GetMaximum()*20)
#   c.Update()
#   c.SaveAs(plot_out_folder+title+'_LOG.pdf')
#   c.SaveAs(plot_out_folder+title+'_LOG.png')



# def drawPlotSimultaneous(mass,multidata,norm,multimodel,fitRes,region,category,title):

#   c = TCanvas('c','c',800,800)
#   c.Divide(1,2)
#   setTopPad(c.GetPad(1),RATIO)
#   setBotPad(c.GetPad(2),RATIO)

#   # draw top pad
#   c.cd(1)
#   frame = mass.frame()
#   frame.GetYaxis().SetTitle('Events / x GeV')
#   frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame,1.1)
#   data_obs = multidata.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Cut('reg==reg::'+category), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
#   multimodel.plotOn(frame, RooFit.Slice(region, category), RooFit.ProjWData(RooArgSet(region), multidata), RooFit.LineColor(2), RooFit.Range('fitrangemass'),RooFit.Name('tot'))
#   chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
#   res  = frame.pullHist()
#   multimodel.plotOn(frame, RooFit.Slice(region, category), RooFit.ProjWData(RooArgSet(region), multidata), RooFit.Range('fitrangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.FillColor(1),RooFit.FillStyle(3005),RooFit.LineColor(0))
#   frame.GetYaxis().SetTitleOffset(1.1)
#   frame.SetMaximum(frame.GetMaximum()*1.1)
#   frame.Draw()
#   drawALICE('','Internal')
#   leg = TLegend(0.5, 0.65, 0.95, 0.9)
#   leg.SetBorderSize(0)
#   leg.SetFillStyle(0) #1001
#   leg.SetFillColor(0)
#   leg.AddEntry('data_obs', category,                         'LP')
#   leg.AddEntry('tot',      'fit',                            'L' )
#   leg.AddEntry('Unc',      'fit uncertainty',                'F' )
#   leg.Draw()

#   # draw bottom pad
#   c.cd(2)
#   frame_res = mass.frame()
#   frame_res.GetYaxis().SetTitle('Events / x GeV')
#   frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame_res,1.1)
#   setBotStyle(frame_res,RATIO)
#   frame_res.addPlotable(res,'P')
#   frame_res.GetYaxis().SetRangeUser(-5,5)
#   frame_res.GetYaxis().SetTitle('(data-fit)/#sigma_{data}')
#   frame_res.GetYaxis().CenterTitle()
#   frame_res.GetYaxis().SetTitleOffset(0.4)
#   frame_res.Draw()
#   line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
#   drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)

#   # save canvas in LIN scale
#   c.SaveAs(plot_out_folder+title+'.pdf')
#   c.SaveAs(plot_out_folder+title+'.png')
#   c.SaveAs(plot_out_folder+title+'.root')

#   # update and save canvas in LOG scale
#   c.GetPad(1).SetLogy()
#   frame.SetMinimum(max(frame.GetMinimum(),1.))
#   frame.SetMaximum(frame.GetMaximum()*20)
#   c.Update()
#   c.SaveAs(plot_out_folder+title+'_LOG.pdf')
#   c.SaveAs(plot_out_folder+title+'_LOG.png')


# def drawAlphaPlot(mass, alpha, crmodel, srmodel, fitRes, alpha2=None, crmodel2=None, srmodel2=None, fitRes2=None):

#   c = TCanvas('c','c',800,800)

#   # draw alpha canvas
#   frame = mass.frame()
#   frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame,1.1)
#   alpha.plotOn(frame, RooFit.VisualizeError(fitRes, 2, False), RooFit.LineColor(800), RooFit.FillColor(800), RooFit.Name('2sigma'))
#   alpha.plotOn(frame, RooFit.VisualizeError(fitRes, 1, False), RooFit.LineColor(413), RooFit.FillColor(413), RooFit.Name('1sigma'))
#   alpha.plotOn(frame, RooFit.LineColor(1), RooFit.Name('alpha'))
#   if alpha2:
#       alpha2.plotOn(frame, RooFit.LineColor(921), RooFit.LineStyle(7), RooFit.Name('alpha2'))
#   frame.GetYaxis().SetTitleOffset(1.1)
#   frame.SetMaximum(frame.GetMaximum()*1.1)
#   frame.SetMinimum(frame.GetMaximum()*0.6)
#   frame.Draw()
#   drawALICE('','Internal')
#   leg = TLegend(0.5, 0.65, 0.95, 0.9)
#   leg.SetBorderSize(0)
#   leg.SetFillStyle(0) #1001
#   leg.SetFillColor(0)
#   leg.AddEntry('alpha', 'F function', 'L')
#   leg.AddEntry('1sigma', 'F function #pm 1 std. deviation', 'F')
#   leg.AddEntry('2sigma', 'F function #pm 2 std. deviations', 'F')
#   if alpha2: 
#       leg.AddEntry('alpha2', alpha2.GetTitle(), 'L')
#   leg.Draw()
#   frame.GetYaxis().SetTitle('a.u.')
#   # save canvas in LIN scale
#   c.SaveAs(plot_out_folder+'alpha.pdf')
#   c.SaveAs(plot_out_folder+'alpha.png')
#   c.SaveAs(plot_out_folder+'alpha.root')

#   # draw alpha + pre-post function canvas
#   c_1 = TCanvas('c_1','c_1',800,800)  
#   frame_1 = mass.frame()
#   frame_1.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame_1, 1.1)
#   alpha.plotOn(frame_1, RooFit.LineColor(1), RooFit.Name('alpha'))
#   crmodel.plotOn(frame_1, RooFit.LineColor(602), RooFit.Name('crmodel')) 
#   srmodel.plotOn(frame_1, RooFit.LineColor(2), RooFit.Name('srmodel')) 
#   if alpha2:
#     alpha2.plotOn(frame_1, RooFit.LineColor(921), RooFit.LineStyle(7), RooFit.Name('alpha2'))
#     crmodel2.plotOn(frame_1, RooFit.LineColor(602), RooFit.LineStyle(7), RooFit.Name('crmodel2')) 
#     srmodel2.plotOn(frame_1, RooFit.LineColor(2), RooFit.LineStyle(7), RooFit.Name('srmodel2')) 
#   frame_1.GetYaxis().SetTitleOffset(1.1)
#   frame_1.SetMinimum(1.e-5)
#   frame_1.Draw()
#   drawALICE('','Internal')
#   leg1 = TLegend(0.5, 0.55, 0.95, 0.9)
#   leg1.SetBorderSize(0)
#   leg1.SetFillStyle(0) #1001
#   leg1.SetFillColor(0)
#   leg1.AddEntry('alpha',    'F function', 'L')
#   leg1.AddEntry('crmodel', 'fit in MixPM', 'L')
#   leg1.AddEntry('srmodel', 'comb. predict. in DataPM', 'L')
#   if alpha2:
#     leg1.AddEntry('alpha2',    'alternative F function', 'L')
#     leg1.AddEntry('crmodel2', 'alternative fit in MixPM', 'L')
#     leg1.AddEntry('srmodel2', 'alternative comb. predict. in DataPM', 'L')
#   leg1.Draw()
#   frame_1.GetYaxis().SetTitle('a.u.')
#   c_1.Update()
#   # save canvas in LIN scale
#   c_1.SaveAs(plot_out_folder+'alpha_control.pdf')
#   c_1.SaveAs(plot_out_folder+'alpha_control.png')
#   c_1.SaveAs(plot_out_folder+'alpha_control.root')
#   # update and save canvas in LOG scale
#   # frame_1.SetMaximum(1.)
#   c_1.GetPad(0).SetLogy()
#   c_1.Update()
#   c_1.SaveAs(plot_out_folder+'alpha_control_LOG.pdf')
#   c_1.SaveAs(plot_out_folder+'alpha_control_LOG.png')


#   # draw predition in SR canvas
#   c_2 = TCanvas('c_2', 'c_2', 800, 800)
#   c_2.cd()
#   frame_2 = mass.frame()
#   frame_2.GetYaxis().SetTitle('a.u.')
#   frame_2.GetXaxis().SetTitle('Mass #mu#mu (GeV)')  
#   setHistStyle(frame_2, 1.1)
#   srmodel.plotOn(frame_2, RooFit.VisualizeError(fitRes, 1, False), RooFit.FillColor(2), RooFit.FillStyle(3004), RooFit.LineColor(2), RooFit.Name('srmodel')) 
#   srmodel.plotOn(frame_2, RooFit.LineColor(2), RooFit.Name('srmodel')) 
#   if alpha2:
#       srmodel2.plotOn(frame_2, RooFit.VisualizeError(fitRes, 1, False), RooFit.FillColor(636), RooFit.FillStyle(3005), RooFit.LineColor(636), RooFit.LineStyle(7), RooFit.Name('srmodel2')) 
#       srmodel2.plotOn(frame_2, RooFit.LineColor(636), RooFit.LineStyle(7), RooFit.Name('srmodel2')) 
#   frame_2.GetYaxis().SetTitleOffset(1.1)
#   frame_2.SetMinimum(1.e-5)
#   frame_2.Draw()
#   drawALICE('','Internal')
#   leg2 = TLegend(0.5, 0.55, 0.95, 0.9)
#   leg2.SetBorderSize(0)
#   leg2.SetFillStyle(0) #1001
#   leg2.SetFillColor(0)
#   leg2.AddEntry('srmodel', 'comb. predict. in DataPM', 'L')
#   if alpha2:
#       leg2.AddEntry('srmodel2', 'alternative comb. predict. in DataPM', 'L')
#   leg2.Draw()
#   frame_2.GetYaxis().SetTitle('a.u.')
#   c_2.Update()  
#   # save canvas in LIN scale  
#   c_2.SaveAs(plot_out_folder+'alpha_SR.pdf')
#   c_2.SaveAs(plot_out_folder+'alpha_SR.png')
#   c_2.SaveAs(plot_out_folder+'alpha_SR.root')
#   # update and save canvas in LOG scale
#   # frame_2.SetMaximum(1.)
#   c_2.GetPad(0).SetLogy()
#   c_2.Update()
#   c_2.SaveAs(plot_out_folder+'alpha_SR_LOG.pdf')
#   c_2.SaveAs(plot_out_folder+'alpha_SR_LOG.png')




    
# def drawComparisonPlot(mass,data1,norm1,data2,norm2,title):

#   c = TCanvas('c','c',800,800)
#   c.Divide(1,2)
#   setTopPad(c.GetPad(1),RATIO)
#   setBotPad(c.GetPad(2),RATIO)

#   # draw top pad
#   c.cd(1)
#   frame = mass.frame()
#   frame.GetYaxis().SetTitle('Events / x GeV')
#   frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame,1.1)
#   data1_obs = data1.plotOn(frame, RooFit.Normalization(norm1, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data1_obs'),RooFit.LineColor(61),RooFit.MarkerColor(61))
#   data2_obs = data2.plotOn(frame, RooFit.Normalization(norm2, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data2_obs'),RooFit.LineColor(98),RooFit.MarkerColor(98))
#   frame.GetYaxis().SetTitleOffset(1.1)
#   frame.Draw("E2")
#   frame.SetMaximum(frame.GetMaximum()*1.1)  
#   leg = TLegend(0.5, 0.65, 0.95, 0.9)
#   leg.SetBorderSize(0)
#   leg.SetFillStyle(0) #1001
#   leg.SetFillColor(0)
#   leg.AddEntry('data1_obs', data1.GetName().replace('_',' '), 'LP')
#   leg.AddEntry('data2_obs', data2.GetName().replace('_',' '), 'LP')
#   leg.Draw()
#   drawALICE('','Internal')

#   # draw bottom pad
#   c.cd(2)
#   hist1 = data1.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'hist1').ProjectionX('hist1')
#   hist2 = data2.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'hist2').ProjectionX('hist2')
#   hist1.Sumw2(True)
#   hist2.Sumw2(True)
#   ratio = hist1.Clone('ratio')
#   ratio.Divide(hist2)
#   ratio.GetYaxis().SetTitle('Events / x GeV')
#   ratio.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(ratio,1.1)
#   setBotStyle(ratio,RATIO,False)
#   ratio.SetFillStyle(3005)
#   ratio.SetFillColor(1)
#   ratio.GetYaxis().SetTitle('ratio')
#   ratio.GetYaxis().CenterTitle()
#   ratio.GetYaxis().SetTitleOffset(0.4)
#   ratio.Draw("E2")
#   r,e = drawRatio(hist1, hist2, isBkg=False)
#   line_res = drawLine(ratio.GetXaxis().GetXmin(),r,ratio.GetXaxis().GetXmax(),r)
#   ratio.GetYaxis().SetRangeUser(0.9*r,1.1*r)
#   if r+0.5 < 1 : ratio.GetYaxis().SetMaxDigits(3)

#   # save canvas in LIN scale
#   c.SaveAs(plot_out_folder+title+'.pdf')
#   c.SaveAs(plot_out_folder+title+'.png')
#   c.SaveAs(plot_out_folder+title+'.root')

#   # update and save canvas in LOG scale
#   c.GetPad(1).SetLogy()
#   frame.SetMinimum(max(frame.GetMinimum(),1.))
#   frame.SetMaximum(frame.GetMaximum()*20)
#   c.Update()
#   c.SaveAs(plot_out_folder+title+'_LOG.pdf')
#   c.SaveAs(plot_out_folder+title+'_LOG.png')


# def drawDoubleComparisonPlot(mass,dataA1,normA1,dataA2,normA2,dataB1,normB1,dataB2,normB2,title):


#   # draw all hists
#   c = TCanvas('c','c',800,800)
#   frame = mass.frame()
#   frame.GetYaxis().SetTitle('Events / x GeV')
#   frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   setHistStyle(frame,1.1)
#   dataA1_obs = dataA1.plotOn(frame, RooFit.Normalization(normA1, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataA1_obs'),RooFit.LineColor(61),RooFit.MarkerColor(61))
#   dataA2_obs = dataA2.plotOn(frame, RooFit.Normalization(normA2, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataA2_obs'),RooFit.LineColor(98),RooFit.MarkerColor(98))
#   dataB1_obs = dataB1.plotOn(frame, RooFit.Normalization(normB1, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataB1_obs'),RooFit.LineColor(66),RooFit.MarkerColor(66))
#   dataB2_obs = dataB2.plotOn(frame, RooFit.Normalization(normB2, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataB2_obs'),RooFit.LineColor(94),RooFit.MarkerColor(94))
#   frame.GetYaxis().SetTitleOffset(1.1)
#   frame.Draw("E2")
#   frame.SetMaximum(frame.GetMaximum()*1.1)  
#   leg = TLegend(0.5, 0.65, 0.95, 0.9)
#   leg.SetBorderSize(0)
#   leg.SetFillStyle(0) #1001
#   leg.SetFillColor(0)
#   nameA1 = dataA1.GetName().replace('_',' ')
#   nameA2 = dataA2.GetName().replace('_',' ')
#   nameB1 = dataB1.GetName().replace('_',' ')
#   nameB2 = dataB2.GetName().replace('_',' ')
#   leg.AddEntry('dataA1_obs', nameA1, 'LP')
#   leg.AddEntry('dataA2_obs', nameA2, 'LP')
#   leg.AddEntry('dataB1_obs', nameB1, 'LP')
#   leg.AddEntry('dataB2_obs', nameB2, 'LP')
#   leg.Draw()
#   drawALICE('','Internal')
#   # save canvas in LIN scale
#   c.SaveAs(plot_out_folder+title+'.pdf')
#   c.SaveAs(plot_out_folder+title+'.png')
#   c.SaveAs(plot_out_folder+title+'.root')
#   # update and save canvas in LOG scale
#   c.GetPad(0).SetLogy()
#   frame.SetMinimum(max(frame.GetMinimum(),1.))
#   frame.SetMaximum(frame.GetMaximum()*20)
#   c.Update()
#   c.SaveAs(plot_out_folder+title+'_LOG.pdf')
#   c.SaveAs(plot_out_folder+title+'_LOG.png')

#   # prepare other histos
#   histA1 = dataA1.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histA1').ProjectionX('histA1')
#   histA2 = dataA2.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histA2').ProjectionX('histA2')
#   histB1 = dataB1.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histB1').ProjectionX('histB1')
#   histB2 = dataB2.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'histB2').ProjectionX('histB2')
#   histA1.Sumw2(True)
#   histA2.Sumw2(True)
#   histB1.Sumw2(True)
#   histB2.Sumw2(True)
  
#   ratioA1A2 = histA1.Clone('ratioA1A2')
#   ratioA1A2.Divide(histA2)

#   ratioB1B2 = histB1.Clone('ratioB1B2')
#   ratioB1B2.Divide(histB2)
  
#   ratioA1B1 = histA1.Clone('ratioA1B1')
#   ratioA1B1.Divide(histB1)
  
#   ratioA2B2 = histA2.Clone('ratioA2B2')
#   ratioA2B2.Divide(histB2)

#   doubleratioA1A2B1B2 = ratioA1A2.Clone('doubleratioA1A2B1B2')
#   doubleratioA1A2B1B2.Divide(ratioB1B2)
  
#   doubleratioA1B1A2B2 = ratioA1B1.Clone('doubleratioA1B1A2B2')
#   doubleratioA1B1A2B2.Divide(ratioA2B2)

#   diffA1A2B1B2 = ratioA1A2.Clone('diffA1A2B1B2')
#   sumA1A2B1B2  = ratioA1A2.Clone('sumA1A2B1B2')
#   diffA1A2B1B2.Add(ratioB1B2,-1.)
#   sumA1A2B1B2.Add(ratioB1B2,+1.)
#   asymmetryA1A2B1B2 = diffA1A2B1B2.Clone('asymmetryA1A2B1B2')
#   asymmetryA1A2B1B2.Divide(sumA1A2B1B2)

#   diffA1B1A2B2 = ratioA1B1.Clone('diffA1B1A2B2')
#   sumA1B1A2B2  = ratioA1B1.Clone('sumA1B1A2B2')
#   diffA1B1A2B2.Add(ratioA2B2,-1.)
#   sumA1B1A2B2.Add(ratioA2B2,+1.)
#   asymmetryA1B1A2B2 = diffA1B1A2B2.Clone('asymmetryA1B1A2B2')
#   asymmetryA1B1A2B2.Divide(sumA1B1A2B2)

#   ptbin  = 'pT{}'.format(nameA1.split('pT')[1]) if 'pT' in title else ''
#   nameA1 = nameA1.replace('data','').split('pT')[0] if 'pT' in title else nameA1
#   nameA2 = nameA2.replace('data','').split('pT')[0] if 'pT' in title else nameA2
#   nameB1 = nameB1.replace('data','').split('pT')[0] if 'pT' in title else nameB1
#   nameB2 = nameB2.replace('data','').split('pT')[0] if 'pT' in title else nameB2

#   # draw ratioA1A2/ratioB1B2
#   c_1 = TCanvas('c_1','c_1',800,800)
#   c_1.cd()
#   setHistStyle(doubleratioA1A2B1B2,1.1)  
#   doubleratioA1A2B1B2.SetFillColor(1)
#   doubleratioA1A2B1B2.SetFillStyle(3005)
#   doubleratioA1A2B1B2.Draw('E2')
#   doubleratioA1A2B1B2.GetYaxis().SetTitle('ratio')
#   doubleratioA1A2B1B2.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   line_1 = drawLine(doubleratioA1A2B1B2.GetXaxis().GetXmin(),1.,doubleratioA1A2B1B2.GetXaxis().GetXmax(),1.)  
#   drawALICE('','Internal')
#   drawText(0.7, 0.85, '#frac{%s/%s}{%s/%s} '%(nameA1,nameA2,nameB1,nameB2))
#   drawText(0.4, 0.85,  '%s'%ptbin)
#   # save canvas in LIN scale
#   c_1.SaveAs(plot_out_folder+title+'_1.pdf')
#   c_1.SaveAs(plot_out_folder+title+'_1.png')
#   c_1.SaveAs(plot_out_folder+title+'_1.root') 
#   # save canvas in ZOOM scale
#   doubleratioA1A2B1B2.GetYaxis().SetRangeUser(0.9,1.1)
#   c_1.SaveAs(plot_out_folder+title+'_1_ZOOM.pdf')
#   c_1.SaveAs(plot_out_folder+title+'_1_ZOOM.png')


#   # draw ratioA1B1/ratioA2B2
#   c_2 = TCanvas('c_2','c_2',800,800)
#   c_2.cd()
#   setHistStyle(doubleratioA1B1A2B2,1.1)  
#   doubleratioA1B1A2B2.SetFillColor(1)
#   doubleratioA1B1A2B2.SetFillStyle(3005)
#   doubleratioA1B1A2B2.Draw('E2')
#   doubleratioA1B1A2B2.GetYaxis().SetTitle('ratio')
#   doubleratioA1B1A2B2.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   line_1 = drawLine(doubleratioA1A2B1B2.GetXaxis().GetXmin(),1.,doubleratioA1A2B1B2.GetXaxis().GetXmax(),1.)  
#   drawALICE('','Internal')
#   drawText(0.7, 0.85, '#frac{%s/%s}{%s/%s} '%(nameA1,nameB1,nameA2,nameB2))
#   drawText(0.4, 0.85,  '%s'%ptbin)
#   # save canvas in LIN scale
#   c_2.SaveAs(plot_out_folder+title+'_2.pdf')
#   c_2.SaveAs(plot_out_folder+title+'_2.png')
#   c_2.SaveAs(plot_out_folder+title+'_2.root')
#   # save canvas in ZOOM scale
#   doubleratioA1B1A2B2.GetYaxis().SetRangeUser(0.9,1.1)
#   c_2.SaveAs(plot_out_folder+title+'_2_ZOOM.pdf')
#   c_2.SaveAs(plot_out_folder+title+'_2_ZOOM.png')


#   # draw asymmetryA1A2B1B2
#   c_3 = TCanvas('c_3','c_3',800,800)
#   c_3.cd()
#   setHistStyle(asymmetryA1A2B1B2,1.1)  
#   asymmetryA1A2B1B2.SetFillColor(1)
#   asymmetryA1A2B1B2.SetFillStyle(3005)
#   asymmetryA1A2B1B2.Draw('E2')
#   asymmetryA1A2B1B2.GetYaxis().SetTitle('asymmetry')
#   asymmetryA1A2B1B2.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   line_0 = drawLine(doubleratioA1A2B1B2.GetXaxis().GetXmin(),0.,doubleratioA1A2B1B2.GetXaxis().GetXmax(),0.)  
#   drawALICE('','Internal')
#   drawText(0.7, 0.85, '#frac{%s-%s}{%s+%s} '%(nameA1,nameA2,nameB1,nameB2))  
#   drawText(0.4, 0.85,  '%s'%ptbin)
#   # save canvas in LIN scale
#   c_3.SaveAs(plot_out_folder+title+'_3.pdf')
#   c_3.SaveAs(plot_out_folder+title+'_3.png')
#   c_3.SaveAs(plot_out_folder+title+'_3.root')
#   # save canvas in ZOOM scale
#   asymmetryA1A2B1B2.GetYaxis().SetRangeUser(-5e-2,5e-2)
#   c_3.SaveAs(plot_out_folder+title+'_3_ZOOM.pdf')
#   c_3.SaveAs(plot_out_folder+title+'_3_ZOOM.png')


#   # draw asymmetryA1B1A2B2
#   c_4 = TCanvas('c_4','c_4',800,800)
#   c_4.cd()
#   setHistStyle(asymmetryA1B1A2B2,1.1)  
#   asymmetryA1B1A2B2.SetFillColor(1)
#   asymmetryA1B1A2B2.SetFillStyle(3005)
#   asymmetryA1B1A2B2.Draw('E2')
#   asymmetryA1B1A2B2.GetYaxis().SetTitle('asymmetry')
#   asymmetryA1B1A2B2.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
#   line_0 = drawLine(doubleratioA1A2B1B2.GetXaxis().GetXmin(),0.,doubleratioA1A2B1B2.GetXaxis().GetXmax(),0.)  
#   drawALICE('','Internal')
#   drawText(0.7, 0.85, '#frac{%s-%s}{%s+%s} '%(nameA1,nameB1,nameA2,nameB2))  
#   drawText(0.4, 0.85,  '%s'%ptbin)
#   # save canvas in LIN scale
#   c_4.SaveAs(plot_out_folder+title+'_4.pdf')
#   c_4.SaveAs(plot_out_folder+title+'_4.png')
#   c_4.SaveAs(plot_out_folder+title+'_4.root')
#   # save canvas in ZOOM scale
#   asymmetryA1B1A2B2.GetYaxis().SetRangeUser(-5e-2,5e-2)
#   c_4.SaveAs(plot_out_folder+title+'_4_ZOOM.pdf')
#   c_4.SaveAs(plot_out_folder+title+'_4_ZOOM.png')









# if __name__ == '__main__':
#   # fitter(pT_bins[0])
#   fitter(None)



# exit()


# # # data.Print()

# # ### JPSI MODEL ###
# # meanJPsi    = RooRealVar('meanJPsi',  'meanJPsi',   3.1,  3.0,  3.2 )
# # sigmaJPsi   = RooRealVar('sigmaJPsi', 'sigmaJPsi',  5e-3, 1e-3, 5e-2)
# # alpha1JPsi  = RooRealVar('alpha1JPsi','alpha1JPsi', 0.8,   0.5,  1.5 )
# # n1JPsi      = RooRealVar('n1JPsi',    'n1JPsi',     4,    1.5,  7   )
# # alpha2JPsi  = RooRealVar('alpha2JPsi','alpha2JPsi', 0.8,   0.5,  2   )
# # n2JPsi      = RooRealVar('n2JPsi',    'n2JPsi',     4,    1,    15  )

# # # JPsiModel  = RooGaussian('JPsiModel','JPsiModel',mass,meanJPsi,sigmaJPsi)
# # JPsiModel = RooCBShape('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi)
# # # JPsiModel   = RooDoubleCrystalBall('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

# # ## PSIp MODEL ###
# # DmeanJPsiPsip  = RooRealVar('DmeanJPsiPsip','DmeanJPsiPsip',0.590,0.580,0.600)
# # RsigmaJPsiPsip = RooRealVar('RsigmaJPsiPsip','RsigmaJPsiPsip',1.02,1.00,1.10)
# # meanPsip  = RooFormulaVar('meanPsip','meanPsip','@0+@1',RooArgList(meanJPsi,DmeanJPsiPsip))
# # sigmaPsip = RooFormulaVar('sigmaPsip','sigmaPsip','@0*@1',RooArgList(sigmaJPsi,RsigmaJPsiPsip))

# # # PsipModel  = RooGaussian('PsipModel','PsipModel',mass,meanPsip,sigmaPsip)
# # PsipModel   = RooCBShape('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi)
# # # PsipModel   = RooDoubleCrystalBall('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

# # ### COMB BKG EXPO MODEL ###
# # constComb = RooRealVar('constComb','constComb',-1,-2,0)
# # BkgCombModel   = RooExponential('BkgCombModel','BkgCombModel',mass,constComb)

# # # ### COMB BKG EXPO MODEL ###
# # # constComb  = RooRealVar('constComb','constComb',-2,-2,-0.5)
# # # offsetComb = RooRealVar('offsetComb','offsetComb',1.1,0.8,1.4)
# # # widthComb  = RooRealVar('widthComb','widthComb',0.5,0.3,0.8)
# # # BkgCombModel   = RooErfExpPdf('BkgCombModel','BkgCombModel',mass,constComb,offsetComb,widthComb)

# # # ### COMB BKG 2 EXPO MODEL ###
# # # const1Comb = RooRealVar('const1Comb','const1Comb',-2,-2,0)
# # # const2Comb = RooRealVar('const2Comb','const2Comb',-0.5,-1,-1e-2)
# # # fComb = RooRealVar('fComb','fComb',0.1,0,0.5)
# # # BkgCombModel   = Roo2ExpPdf('BkgCombModel','BkgCombModel',mass,const1Comb,const2Comb,fComb)

# # # ### COMB BKG EXPN MODEL ###
# # # numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
# # # denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e-3,   1e-6,   1e2)
# # # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)

# # # ### COMB BKG EXPTAIL MODEL ###
# # # numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
# # # denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e3,   1e0,   1e5)
# # # # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)
# # # BkgCombModel= RooExpTailPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)


# # # NJPsi     = RooRealVar('NJPsi','NJPsi',1e3,1e1,1e6)
# # # JPsiModelExt = RooExtendPdf('JPsiModelExt','JPsiModelExt',JPsiModel,NJPsi)
# # # NPsip     = RooRealVar('NPsip','NPsip',1e1,1e-1,1e4)
# # # PsipModelExt = RooExtendPdf('PsipModelExt','PsipModelExt',PsipModel,NPsip)
# # # NBkgComb = RooRealVar('NBkgComb','NBkgComb',  1e6,   1e1,   1e8)
# # # BkgCombModelExt = RooExtendPdf('BkgCombModelExt','BkgCombModelExt',BkgCombModel,NBkgComb)


# # fJPsi     = RooRealVar('fJPsi','fJPsi',1e-2,1e-3,1e0)
# # fPsip     = RooRealVar('fPsip','fPsip',1e-5,1e-6,1e-2)


# # DmeanJPsiPsip.setConstant(True)
# # RsigmaJPsiPsip.setConstant(True)

# # # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel),RooArgList(frac))
# # # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt,PsipModelExt),RooArgList(NBkgComb,NJPsi,NPsip))
# # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(JPsiModel,PsipModel,BkgCombModel),RooArgList(fJPsi,fPsip))
# # # baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt),RooArgList(NBkgComb,NJPsi))
# # fitRes    = baseModel.fitTo(data,
# #                             RooFit.Extended(True),
# #                             RooFit.SumW2Error(True),
# #                             RooFit.Strategy(2),
# #                             RooFit.Minimizer('Minuit'),
# #                             RooFit.Range('rangemass'),
# #                             RooFit.Save(1),
# #                             RooFit.NumCPU(12),
# #                             RooFit.PrintLevel(1 if options.verbose else -1) )

# # print ('=== OVERALL FIT RESULTS ===')
# # fitRes.Print()
# # print ('===========================')

# # RATIO = 3.75
# # c = TCanvas('c','c',800,800)
# # c.Divide(1,2)
# # setTopPad(c.GetPad(1),RATIO)
# # setBotPad(c.GetPad(2),RATIO)
# # # draw top pad
# # c.cd(1)
# # frame = mass.frame()
# # frame.GetYaxis().SetTitle('Events / x GeV')
# # frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
# # setHistStyle(frame,1.1)
# # data_obs = data.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
# # baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('tot'))
# # chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
# # res  = frame.pullHist()
# # # baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(2), RooFit.Components('BkgCombModelExt'), RooFit.Range('rangemass'),RooFit.Name('comb'))
# # # baseModel.plotOn(frame, RooFit.LineColor(4), RooFit.LineStyle(3), RooFit.Components('JPsiModelExt'), RooFit.Range('rangemass'),RooFit.Name('JPsi'))
# # # baseModel.plotOn(frame, RooFit.LineColor(3), RooFit.LineStyle(4), RooFit.Components('PsipModelExt'), RooFit.Range('rangemass'),RooFit.Name('Psip'))
# # baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(2), RooFit.Components('BkgCombModel'), RooFit.Range('rangemass'),RooFit.Name('comb'))
# # baseModel.plotOn(frame, RooFit.LineColor(4), RooFit.LineStyle(3), RooFit.Components('JPsiModel'), RooFit.Range('rangemass'),RooFit.Name('JPsi'))
# # baseModel.plotOn(frame, RooFit.LineColor(3), RooFit.LineStyle(4), RooFit.Components('PsipModel'), RooFit.Range('rangemass'),RooFit.Name('Psip'))
# # # baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.FillColor(1),RooFit.FillStyle(3005))
# # frame.GetYaxis().SetTitleOffset(1.1)
# # frame.Draw()
# # drawALICE('','Internal')
# # # draw bottom pad
# # c.cd(2)
# # frame_res = mass.frame()
# # frame_res.GetYaxis().SetTitle('Events / x GeV')
# # frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
# # setHistStyle(frame_res,1.1)
# # setBotStyle(frame_res,RATIO)
# # frame_res.addPlotable(res,'P')
# # frame_res.GetYaxis().SetRangeUser(-5,5)
# # frame_res.GetYaxis().SetTitle('(data-fit)/#sigma_{data}')
# # frame_res.GetYaxis().CenterTitle()
# # frame_res.GetYaxis().SetTitleOffset(0.4)
# # frame_res.Draw()
# # line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
# # drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)
# # c.SaveAs(plot_out_folder+'simpleFit.pdf')
# # c.SaveAs(plot_out_folder+'simpleFit.png')
# # c.SaveAs(plot_out_folder+'simpleFit.root')

# # c.SetLogy()
# # c.SaveAs(plot_out_folder+'simpleFit_LOG.pdf')
# # c.SaveAs(plot_out_folder+'simpleFit_LOG.png')

# # file.Close()

# # exit()



# # ### pT BINNED FIT 

# # file = TFile('../workspaces/processed_w_C.root','read')
# # w = file.Get('processed_w_C')
# # w.Print()

# # for it in pT_bins:
# #   mass = w.var('tMass')
# #   data = w.data('data_mix_pT%s'%it)
# #   mass_range = (1.5, 8.)
# #   binsmass = RooBinning(options.n_bins, mass_range[0], mass_range[1])
# #   mass.setRange('rangemass',mass_range[0], mass_range[1])

# #   meanJPsi  = RooRealVar('meanJPsi','meanJPsi',3.1,2.5,3.6)
# #   sigmaJPsi = RooRealVar('sigmaJPsi','sigmaJPsi',5e-2,1e-3,1)
# #   gausJPsi  = RooGaussian('gausJPsi','gausJPsi',mass,meanJPsi,sigmaJPsi)

# #   constComb = RooRealVar('constComb','constComb',-1,-2,0)
# #   expComb   = RooExponential('expComb','expComb',mass,constComb)

# #   frac      = RooRealVar('frac','frac',0.5,0,1)
# #   baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(expComb,gausJPsi),RooArgList(frac))
# #   fitRes    = baseModel.fitTo(data,RooFit.SumW2Error(True),RooFit.Range('rangemass'),RooFit.Save(1),RooFit.NumCPU(10),RooFit.PrintLevel(-1))

# #   print ('===  pT%s FIT RESULTS ==='%it)
# #   fitRes.Print()
# #   print ('===========================')
# #   print ('')
# #   print ('')


# #   RATIO = 3.75
# #   c = TCanvas('c','c',800,800)
# #   c.Divide(1,2)
# #   setTopPad(c.GetPad(1),RATIO)
# #   setBotPad(c.GetPad(2),RATIO)
# #   # draw top pad
# #   c.cd(1)
# #   frame = mass.frame()
# #   frame.GetYaxis().SetTitle('Events / x GeV')
# #   frame.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
# #   setHistStyle(frame,1.1)
# #   data_obs = data.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
# #   baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('tot'))
# #   chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
# #   res  = frame.pullHist()
# #   baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(2), RooFit.Components('expComb'), RooFit.Range('rangemass'),RooFit.Name('comb'))
# #   baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.FillColor(1),RooFit.FillStyle(3005))
# #   frame.GetYaxis().SetTitleOffset(1.1)
# #   frame.Draw()
# #   drawALICE('','Internal')
# #   drawRegion('pT%s'%it)
# #   # draw bottom pad
# #   c.cd(2)
# #   frame_res = mass.frame()
# #   frame_res.GetYaxis().SetTitle('Events / x GeV')
# #   frame_res.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
# #   setHistStyle(frame_res,1.1)
# #   setBotStyle(frame_res,RATIO)
# #   frame_res.addPlotable(res,'P')
# #   frame_res.GetYaxis().SetRangeUser(-5,5)
# #   frame_res.GetYaxis().SetTitle('(data-fit)/#sigma_{data}')
# #   frame_res.GetYaxis().CenterTitle()
# #   frame_res.GetYaxis().SetTitleOffset(0.4)
# #   frame_res.Draw()
# #   line_res = drawLine(frame_res.GetXaxis().GetXmin(),0.,frame_res.GetXaxis().GetXmax(),0.)
# #   drawText(0.75, 0.95, '#chi^{2}/ndof = %.2f'%chi2, 0.10)
# #   c.SaveAs(plot_out_folder+'simpleFit_pT%s.pdf'%it)
# #   c.SaveAs(plot_out_folder+'simpleFit_pT%s.png'%it)
























# #     # # define the three categories (not ranges) to be used in the fit
# #     # reg = RooCategory('reg', 'reg')
# #     # reg.defineType('mcSR')
# #     # reg.defineType('mcSB')
# #     # reg.defineType('dataSB')
# #     # # combine all the datasets with the corresponding categories together in one big sample
# #     # setMulti = RooDataSet('setMulti', 'setMulti', variables, RooFit.Index(reg), RooFit.WeightVar(weight), RooFit.Import('mcSR', setVjetSR), RooFit.Import('mcSB', setVjetSB), RooFit.Import('dataSB', setDataSB))
    
# #     # # define the simultaneous object by linking to the big-dataset the three PDF that will be used in each of the regions
# #     # simObj = RooSimultaneous('simObj', 'simultaneous pdf', RooArgList(model['VjetSR'], model['VjetSB'], model['BkgSB']), reg)
       
# #     # frSim = simObj.fitTo(setMulti, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(30))
    
# #     # # Print
# #     # if VERBOSE: print '********** Fit result [SIMULTANEOUS] ***'+'*'*40, '\n', frSim.Print(), '\n', '*'*80
# #     # if WITHPAUSE: raw_input('Press Enter to continue...')

# #     # simObj2 = RooSimultaneous('simObj2', 'alternative simultaneous pdf', RooArgList(model['VjetSR2'], model['VjetSB2'], model['BkgSB2']), reg) 
# #     # if ALTERNATIVE:
# #     #     # Perform alternative simultaneous fit
# #     #     #frSim2 = simObj2.fitTo(setMulti, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(30))
# #     #     frSim2 = simObj2.fitTo(setMulti, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(30))
# #     #     if VERBOSE: print '********** Fit result [SIMULTANEOUS-ALTERNATIVE] ***'+'*'*40, '\n', frSim2.Print(), '\n', '*'*80
# #     #     if WITHPAUSE: raw_input('Press Enter to continue...')