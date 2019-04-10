#! /usr/bin/env python
import os, sys, getopt, multiprocessing, copy, math, itertools
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom, Double
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis, TH1D
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TMath, TGraphAsymmErrors
from ROOT import RooFit, RooRealVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooExtendPdf, RooAddPdf, RooProdPdf, RooNumConvPdf, RooFFTConvPdf

from drawUtils import *

print ('LOADING CUSTOM PDFs')
rc = gSystem.Load('../pdfs/HWWLVJRooPdfs_cxx.so') 
if not rc==0:
  print ('SOMETHING WENT WRONG')
  exit()

from ROOT import RooDoubleCrystalBall, RooExpNPdf, RooExpTailPdf, Roo2ExpPdf, RooErfExpPdf, RooErfExpNPdf

RATIO = 3.75

########## OPTIONS ##########

import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-b', '--batch',      action='store_true', default=False,           dest='batch'                    )
parser.add_option('-v', '--verbose',    action='store_true', default=False,           dest='verbose'                  )
parser.add_option('-d', '--datapath',   action='store',      default='data/2015/',    dest='datapath',  type='string' )
parser.add_option('-w', '--wspath',     action='store',      default='../workspaces/',dest='wspath',    type='string' )
parser.add_option('-D', '--do_draw',    action='store_true', default=False,           dest='do_draw'                  )
parser.add_option('-n', '--n_bins',     action='store',      default=130,             dest='n_bins',    type='int'    )
(options, args) = parser.parse_args()
if options.batch: gROOT.SetBatch(True)

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)

TGaxis.SetMaxDigits(2)
plot_out_folder = '../plots/'
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

pT_bins   = ['23', '34', '45', '56', '67', '78', '89', '90']
mass_range = (1.5, 8.)
binsmass = RooBinning(options.n_bins, mass_range[0], mass_range[1])

### OVERALL FIT 

# file = TFile('../workspaces/processed_w_Full.root','read')
# w = file.Get('processed_w_Full')
# w.Print()

# mass = w.var('tMass')
# data = w.data('data_mix')


def fitter(pTbin):

  file = TFile('../workspaces/ws_C_int.root','read')
  w = file.Get('w_C')
  w.Print()

  mass = w.var('tMass')
  data = w.data('data_PM_pT%s'%pTbin)
  mix  = w.data('dataMix_PM_pT%s'%pTbin)
  Nmix = w.var('entriesMix_PM_pT%s'%pTbin).getVal()


  dataPP  = w.data('data_PP_pT%s'%pTbin)
  mixPP   = w.data('dataMix_PP_pT%s'%pTbin)
  NdataPP = w.var('entries_PP_pT%s'%pTbin).getVal()
  NmixPP  = w.var('entriesMix_PP_pT%s'%pTbin).getVal()

  dataMM  = w.data('data_MM_pT%s'%pTbin)
  mixMM   = w.data('dataMix_MM_pT%s'%pTbin)
  NdataMM = w.var('entries_MM_pT%s'%pTbin).getVal()
  NmixMM  = w.var('entriesMix_MM_pT%s'%pTbin).getVal()

  mass.setRange('rangemass',mass_range[0], mass_range[1])

  R   = 0.5 * Nmix / math.sqrt(NmixPP*NmixMM)
  dRx = 0.5 / math.sqrt(NmixPP*NmixMM)
  dRy = - 0.25 * Nmix * NmixMM / pow(NmixPP*NmixMM,3/2)
  dRz = - 0.25 * Nmix * NmixPP / pow(NmixPP*NmixMM,3/2)
  dR  = math.sqrt ( dRx**2*Nmix + dRy**2*NmixPP + dRz**2*NmixMM ) 

  F     = 2 * R * math.sqrt(NdataPP*NdataMM) / Nmix
  F_bis = math.sqrt( (NdataPP*NdataMM) / (NmixPP*NmixMM)) 
  F_ter = (NdataPP + NdataMM) / (NmixPP + NmixMM)

  print("R    = %f +/- %f"%(R,dR))
  print("F    = %f"%F)
  print("F_bis= %f"%F_bis)
  print("F_ter= %f"%F_ter)

  c = TCanvas('c','c',800,600)
  c.cd()
  frame = mass.frame()
  dataPP.plotOn(frame, RooFit.Normalization(NdataPP, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataPP'),RooFit.LineColor(2),RooFit.MarkerColor(2))
  dataMM.plotOn(frame, RooFit.Normalization(NdataMM, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
  frame.Draw()
  c.SaveAs('c_data.root')

  cMix = TCanvas('cMix','cMix',800,600)
  cMix.cd()
  frameMix = mass.frame()
  mixPP.plotOn(frameMix, RooFit.Normalization(NmixPP, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixPP'),RooFit.LineColor(2),RooFit.MarkerColor(2))
  mixMM.plotOn(frameMix, RooFit.Normalization(NmixMM, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
  frameMix.Draw()
  cMix.SaveAs('c_mix.root')

  cPP = TCanvas('cPP','cPP',800,600)
  cPP.cd()
  framePP = mass.frame()
  mixPP.plotOn(framePP, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixPP'),RooFit.LineColor(2),RooFit.MarkerColor(2))
  dataPP.plotOn(framePP, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataPP'),RooFit.LineColor(4),RooFit.MarkerColor(4))
  framePP.Draw()
  cPP.SaveAs('c_PP.root')

  cMM = TCanvas('cMM','cMM',800,600)
  cMM.cd()
  frameMM = mass.frame()
  mixMM.plotOn(frameMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixMM'),RooFit.LineColor(2),RooFit.MarkerColor(2))
  dataMM.plotOn(frameMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
  frameMM.Draw()
  cMM.SaveAs('c_MM.root')

  drawComparisonPlot(mass,dataPP,NdataPP,dataMM,NdataMM,'data_PP-MM_pT%s'%pTbin)
  drawComparisonPlot(mass,mixPP, NmixPP, mixMM, NmixMM, 'mix_PP-MM_pT%s'%pTbin)
  drawComparisonPlot(mass,dataPP,NdataPP,mixPP,NmixPP,'data-mix_PP_pT%s'%pTbin)
  drawComparisonPlot(mass,dataMM,NdataMM,mixMM,NmixMM,'data-mix_MM_pT%s'%pTbin)

  input("wait here")


  dataPPMM = RooDataSet(dataPP,'dataPPMM')
  dataPPMM.append(dataMM)
  mixPPMM  = RooDataSet(mixPP, 'mixPPMM')
  mixPPMM.append(mixMM)

  NdataPPMM= dataPPMM.numEntries()
  NmixPPMM = mixPPMM.numEntries()


  cPPMM = TCanvas('cPPMM','cPPMM',800,600)
  cPPMM.cd()
  framePPMM = mass.frame()
  mixPPMM.plotOn(framePPMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('mixPPMM'),RooFit.LineColor(2),RooFit.MarkerColor(2))
  dataPPMM.plotOn(framePPMM, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('dataPPMM'),RooFit.LineColor(4),RooFit.MarkerColor(4))
  framePPMM.Draw()
  cPPMM.SaveAs('c_PPMM.root')

  drawComparisonPlot(mass,dataPPMM,NdataPPMM,mixPPMM,NmixPPMM,'data-mix_PPMM_pT%s'%pTbin)

  exit()



  # ### BKG EXPO MODEL ###
  # constDataMix = RooRealVar('constDataMix','constDataMix',-1,-2,0)
  # DataMixModel   = RooExponential('DataMixModel','DataMixModel',mass,constDataMix)
  # ### BKG ERF+EXPO MODEL ###
  # constDataMix  = RooRealVar('constDataMix','constDataMix',-1,-3,0)
  # offsetDataMix = RooRealVar('offsetDataMix','offsetDataMix',0.3,0.,4.)
  # widthDataMix  = RooRealVar('widthDataMix','widthDataMix',0.5,0.3,2.)
  # DataMixModel  = RooErfExpPdf('DataMixModel','DataMixModel',mass,constDataMix,offsetDataMix,widthDataMix)
  # ### BKG 2 EXPO MODEL ###
  # const1DataMix = RooRealVar('const1DataMix','const1DataMix',-2,-2,2)
  # const2DataMix = RooRealVar('const2DataMix','const2DataMix',-0.5,-1,-1e-2)
  # fDataMix = RooRealVar('fDataMix','fDataMix',0.1,0,0.5)
  # DataMixModel   = Roo2ExpPdf('DataMixModel','DataMixModel',mass,const1DataMix,const2DataMix,fDataMix)

  ### BKG EXPN MODEL ###
  numBkgDataMix  = RooRealVar('numBkgDataMix','numBkgDataMix',  -1,  -2, 0)
  denBkgDataMix  = RooRealVar('denBkgDataMix','denBkgDataMix',  1e-3,   1e-6,   1e2)
  DataMixModel   = RooExpNPdf('DataMixModel','DataMixModel',mass,numBkgDataMix,denBkgDataMix)


  # ### BKG ERF+EXPO+TAIL MODEL ###
  # constDataMix  = RooRealVar('constDataMix','constDataMix',-1,-3,0)
  # tailDataMix   = RooRealVar('tailDataMix','tailDataMix',1,-1000.,1000.)
  # offsetDataMix = RooRealVar('offsetDataMix','offsetDataMix',0.3,0.,4.)
  # widthDataMix  = RooRealVar('widthDataMix','widthDataMix',0.5,0.3,2.)
  # DataMixModel  = RooErfExpPdf('DataMixModel','DataMixModel',mass,constDataMix,offsetDataMix,widthDataMix)



  DataMixFitRes    = DataMixModel.fitTo(mix,
                              RooFit.SumW2Error(True),
                              RooFit.Strategy(2),
                              RooFit.Minimizer('Minuit'),
                              RooFit.Range('rangemass'),
                              RooFit.Save(1),
                              RooFit.NumCPU(12),
                              RooFit.PrintLevel(1 if options.verbose else -1) )

  print ('===   DataMix 1D FIT    ===')
  DataMixFitRes.Print()
  print ('===========================')

  drawPlot(mass,mix,Nmix,DataMixModel,DataMixFitRes,'dataMix_PM_pT%s'%pTbin)

  file.Close()


def drawPlot(mass,data,norm,model,fitRes,title):

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
  data_obs = data.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
  model.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('tot'))
  chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
  res  = frame.pullHist()
  model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.FillColor(1),RooFit.FillStyle(3005))
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
  c.SaveAs(plot_out_folder+title+'.pdf')
  c.SaveAs(plot_out_folder+title+'.png')
  c.SaveAs(plot_out_folder+title+'.root')

  c.SetLogy()
  c.SaveAs(plot_out_folder+title+'_LOG.pdf')
  c.SaveAs(plot_out_folder+title+'_LOG.png')



def drawComparisonPlot(mass,data1,norm1,data2,norm2,title):

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
  data1_obs = data1.plotOn(frame, RooFit.Normalization(norm1, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data1_obs'),RooFit.LineColor(61),RooFit.MarkerColor(61))
  data2_obs = data2.plotOn(frame, RooFit.Normalization(norm2, RooAbsReal.NumEvent), RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.Poisson),RooFit.Range('rangemass'),RooFit.Name('data2_obs'),RooFit.LineColor(98),RooFit.MarkerColor(98))
  # res  = frame.pullHist()
  frame.GetYaxis().SetTitleOffset(1.1)
  frame.Draw()
  drawALICE('','Internal')
  # draw bottom pad
  c.cd(2)
  hist1 = data1.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'hist1').ProjectionX('hist1')
  hist2 = data2.createHistogram(mass, mass, options.n_bins, options.n_bins, "", 'hist2').ProjectionX('hist2')
  ratio = hist1.Clone('ratio')
  ratio.Divide(hist2)
  ratio.GetYaxis().SetTitle('Events / x GeV')
  ratio.GetXaxis().SetTitle('Mass #mu#mu (GeV)')
  setHistStyle(ratio,1.1)
  setBotStyle(ratio,RATIO)
  # frame_res.addPlotable(res,'P')
  ratio.GetYaxis().SetTitle('ratio')
  ratio.GetYaxis().CenterTitle()
  ratio.GetYaxis().SetTitleOffset(0.4)
  ratio.Draw()
  r,e = drawRatio(hist1, hist2, isBkg=False)
  # drawKolmogorov(hist1, hist2)
  line_res = drawLine(ratio.GetXaxis().GetXmin(),r,ratio.GetXaxis().GetXmax(),r)
  ratio.GetYaxis().SetRangeUser(max(r-0.5,0),r+0.5)
  c.SaveAs(plot_out_folder+title+'.pdf')
  c.SaveAs(plot_out_folder+title+'.png')
  c.SaveAs(plot_out_folder+title+'.root')

  c.SetLogy()
  c.SaveAs(plot_out_folder+title+'_LOG.pdf')
  c.SaveAs(plot_out_folder+title+'_LOG.png')


if __name__ == '__main__':
  fitter(pT_bins[0])

exit()


# data.Print()

### JPSI MODEL ###
meanJPsi    = RooRealVar('meanJPsi',  'meanJPsi',   3.1,  3.0,  3.2 )
sigmaJPsi   = RooRealVar('sigmaJPsi', 'sigmaJPsi',  5e-3, 1e-3, 5e-2)
alpha1JPsi  = RooRealVar('alpha1JPsi','alpha1JPsi', 0.8,   0.5,  1.5 )
n1JPsi      = RooRealVar('n1JPsi',    'n1JPsi',     4,    1.5,  7   )
alpha2JPsi  = RooRealVar('alpha2JPsi','alpha2JPsi', 0.8,   0.5,  2   )
n2JPsi      = RooRealVar('n2JPsi',    'n2JPsi',     4,    1,    15  )

# JPsiModel  = RooGaussian('JPsiModel','JPsiModel',mass,meanJPsi,sigmaJPsi)
JPsiModel = RooCBShape('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi)
# JPsiModel   = RooDoubleCrystalBall('JPsiModel', 'JPsiModel', mass, meanJPsi, sigmaJPsi, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

## PSIp MODEL ###
DmeanJPsiPsip  = RooRealVar('DmeanJPsiPsip','DmeanJPsiPsip',0.590,0.580,0.600)
RsigmaJPsiPsip = RooRealVar('RsigmaJPsiPsip','RsigmaJPsiPsip',1.02,1.00,1.10)
meanPsip  = RooFormulaVar('meanPsip','meanPsip','@0+@1',RooArgList(meanJPsi,DmeanJPsiPsip))
sigmaPsip = RooFormulaVar('sigmaPsip','sigmaPsip','@0*@1',RooArgList(sigmaJPsi,RsigmaJPsiPsip))

# PsipModel  = RooGaussian('PsipModel','PsipModel',mass,meanPsip,sigmaPsip)
PsipModel   = RooCBShape('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi)
# PsipModel   = RooDoubleCrystalBall('PsipModel', 'PsipModel', mass, meanPsip, sigmaPsip, alpha1JPsi, n1JPsi, alpha2JPsi, n2JPsi)

### COMB BKG EXPO MODEL ###
constComb = RooRealVar('constComb','constComb',-1,-2,0)
BkgCombModel   = RooExponential('BkgCombModel','BkgCombModel',mass,constComb)

# ### COMB BKG EXPO MODEL ###
# constComb  = RooRealVar('constComb','constComb',-2,-2,-0.5)
# offsetComb = RooRealVar('offsetComb','offsetComb',1.1,0.8,1.4)
# widthComb  = RooRealVar('widthComb','widthComb',0.5,0.3,0.8)
# BkgCombModel   = RooErfExpPdf('BkgCombModel','BkgCombModel',mass,constComb,offsetComb,widthComb)

# ### COMB BKG 2 EXPO MODEL ###
# const1Comb = RooRealVar('const1Comb','const1Comb',-2,-2,0)
# const2Comb = RooRealVar('const2Comb','const2Comb',-0.5,-1,-1e-2)
# fComb = RooRealVar('fComb','fComb',0.1,0,0.5)
# BkgCombModel   = Roo2ExpPdf('BkgCombModel','BkgCombModel',mass,const1Comb,const2Comb,fComb)

# ### COMB BKG EXPN MODEL ###
# numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
# denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e-3,   1e-6,   1e2)
# BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)

# ### COMB BKG EXPTAIL MODEL ###
# numBkgComb  = RooRealVar('numBkgComb','numBkgComb',  -1,  -2, 0)
# denBkgComb  = RooRealVar('denBkgComb','denBkgComb',  1e3,   1e0,   1e5)
# # BkgCombModel= RooExpNPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)
# BkgCombModel= RooExpTailPdf('BkgCombModel','BkgCombModel',mass,numBkgComb,denBkgComb)


# NJPsi     = RooRealVar('NJPsi','NJPsi',1e3,1e1,1e6)
# JPsiModelExt = RooExtendPdf('JPsiModelExt','JPsiModelExt',JPsiModel,NJPsi)
# NPsip     = RooRealVar('NPsip','NPsip',1e1,1e-1,1e4)
# PsipModelExt = RooExtendPdf('PsipModelExt','PsipModelExt',PsipModel,NPsip)
# NBkgComb = RooRealVar('NBkgComb','NBkgComb',  1e6,   1e1,   1e8)
# BkgCombModelExt = RooExtendPdf('BkgCombModelExt','BkgCombModelExt',BkgCombModel,NBkgComb)


fJPsi     = RooRealVar('fJPsi','fJPsi',1e-2,1e-3,1e0)
fPsip     = RooRealVar('fPsip','fPsip',1e-5,1e-6,1e-2)


DmeanJPsiPsip.setConstant(True)
RsigmaJPsiPsip.setConstant(True)

# baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModel,JPsiModel),RooArgList(frac))
# baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt,PsipModelExt),RooArgList(NBkgComb,NJPsi,NPsip))
baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(JPsiModel,PsipModel,BkgCombModel),RooArgList(fJPsi,fPsip))
# baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(BkgCombModelExt,JPsiModelExt),RooArgList(NBkgComb,NJPsi))
fitRes    = baseModel.fitTo(data,
                            RooFit.Extended(True),
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

RATIO = 3.75
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
data_obs = data.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('tot'))
chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
res  = frame.pullHist()
# baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(2), RooFit.Components('BkgCombModelExt'), RooFit.Range('rangemass'),RooFit.Name('comb'))
# baseModel.plotOn(frame, RooFit.LineColor(4), RooFit.LineStyle(3), RooFit.Components('JPsiModelExt'), RooFit.Range('rangemass'),RooFit.Name('JPsi'))
# baseModel.plotOn(frame, RooFit.LineColor(3), RooFit.LineStyle(4), RooFit.Components('PsipModelExt'), RooFit.Range('rangemass'),RooFit.Name('Psip'))
baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(2), RooFit.Components('BkgCombModel'), RooFit.Range('rangemass'),RooFit.Name('comb'))
baseModel.plotOn(frame, RooFit.LineColor(4), RooFit.LineStyle(3), RooFit.Components('JPsiModel'), RooFit.Range('rangemass'),RooFit.Name('JPsi'))
baseModel.plotOn(frame, RooFit.LineColor(3), RooFit.LineStyle(4), RooFit.Components('PsipModel'), RooFit.Range('rangemass'),RooFit.Name('Psip'))
# baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.FillColor(1),RooFit.FillStyle(3005))
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

file.Close()

exit()



### pT BINNED FIT 

file = TFile('../workspaces/processed_w_C.root','read')
w = file.Get('processed_w_C')
w.Print()

for it in pT_bins:
  mass = w.var('tMass')
  data = w.data('data_mix_pT%s'%it)
  mass_range = (1.5, 8.)
  binsmass = RooBinning(options.n_bins, mass_range[0], mass_range[1])
  mass.setRange('rangemass',mass_range[0], mass_range[1])

  meanJPsi  = RooRealVar('meanJPsi','meanJPsi',3.1,2.5,3.6)
  sigmaJPsi = RooRealVar('sigmaJPsi','sigmaJPsi',5e-2,1e-3,1)
  gausJPsi  = RooGaussian('gausJPsi','gausJPsi',mass,meanJPsi,sigmaJPsi)

  constComb = RooRealVar('constComb','constComb',-1,-2,0)
  expComb   = RooExponential('expComb','expComb',mass,constComb)

  frac      = RooRealVar('frac','frac',0.5,0,1)
  baseModel = RooAddPdf('baseModel','Test fit model',RooArgList(expComb,gausJPsi),RooArgList(frac))
  fitRes    = baseModel.fitTo(data,RooFit.SumW2Error(True),RooFit.Range('rangemass'),RooFit.Save(1),RooFit.NumCPU(10),RooFit.PrintLevel(-1))

  print ('===  pT%s FIT RESULTS ==='%it)
  fitRes.Print()
  print ('===========================')
  print ('')
  print ('')


  RATIO = 3.75
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
  data_obs = data.plotOn(frame, RooFit.Binning(binsmass), RooFit.DataError(RooAbsData.SumW2),RooFit.Range('rangemass'),RooFit.Name('data_obs'))
  baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('tot'))
  chi2 = frame.chiSquare(fitRes.floatParsFinal().getSize())
  res  = frame.pullHist()
  baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.LineStyle(2), RooFit.Components('expComb'), RooFit.Range('rangemass'),RooFit.Name('comb'))
  baseModel.plotOn(frame, RooFit.LineColor(2), RooFit.Range('rangemass'),RooFit.Name('Unc'),RooFit.VisualizeError(fitRes,1,False),RooFit.SumW2Error(True),RooFit.FillColor(1),RooFit.FillStyle(3005))
  frame.GetYaxis().SetTitleOffset(1.1)
  frame.Draw()
  drawALICE('','Internal')
  drawRegion('pT%s'%it)
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
  c.SaveAs(plot_out_folder+'simpleFit_pT%s.pdf'%it)
  c.SaveAs(plot_out_folder+'simpleFit_pT%s.png'%it)
























    # # define the three categories (not ranges) to be used in the fit
    # reg = RooCategory('reg', 'reg')
    # reg.defineType('mcSR')
    # reg.defineType('mcSB')
    # reg.defineType('dataSB')
    # # combine all the datasets with the corresponding categories together in one big sample
    # setMulti = RooDataSet('setMulti', 'setMulti', variables, RooFit.Index(reg), RooFit.WeightVar(weight), RooFit.Import('mcSR', setVjetSR), RooFit.Import('mcSB', setVjetSB), RooFit.Import('dataSB', setDataSB))
    
    # # define the simultaneous object by linking to the big-dataset the three PDF that will be used in each of the regions
    # simObj = RooSimultaneous('simObj', 'simultaneous pdf', RooArgList(model['VjetSR'], model['VjetSB'], model['BkgSB']), reg)
       
    # frSim = simObj.fitTo(setMulti, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(30))
    
    # # Print
    # if VERBOSE: print '********** Fit result [SIMULTANEOUS] ***'+'*'*40, '\n', frSim.Print(), '\n', '*'*80
    # if WITHPAUSE: raw_input('Press Enter to continue...')

    # simObj2 = RooSimultaneous('simObj2', 'alternative simultaneous pdf', RooArgList(model['VjetSR2'], model['VjetSB2'], model['BkgSB2']), reg) 
    # if ALTERNATIVE:
    #     # Perform alternative simultaneous fit
    #     #frSim2 = simObj2.fitTo(setMulti, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(30))
    #     frSim2 = simObj2.fitTo(setMulti, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.SumW2Error(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(30))
    #     if VERBOSE: print '********** Fit result [SIMULTANEOUS-ALTERNATIVE] ***'+'*'*40, '\n', frSim2.Print(), '\n', '*'*80
    #     if WITHPAUSE: raw_input('Press Enter to continue...')