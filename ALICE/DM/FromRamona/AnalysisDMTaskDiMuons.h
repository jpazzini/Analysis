#ifndef AnalysisDMTaskDiMuons_H
#define AnalysisDMTaskDiMuons_H

#include "AliLog.h"
#include "AliLMREvent.h"
#include "AliLMRMuon.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TAxis.h>
#include <iostream>
#include <THnSparse.h>
#include "AliAnalysisMuMuEventCollection.h"

enum CutShapeType {kLin, kRect, kElps};
enum TriggerOptions {kMSH=1<<0, kMSL=1<<1, kMUL=1<<2, kMLL=1<<3, kMB=1<<4,
  kSingleMuon    = kMSH | kMSL,
  kDimuon        = kMUL | kMLL
};

class AnalysisDMTaskDiMuons : public TObject {
  //class AnalysisDMTaskDiMuons : public AliAnalysisTaskSE {
public:
  enum AnalysisOption {kDirect, kMixing};
  enum Pool{kActivity, kZvertex, kNTPool};

  AnalysisDMTaskDiMuons();
  ~AnalysisDMTaskDiMuons();

  void BuildPoolForMixing();

  // Set the trigger for both the analysis and the mixing
  void SetTrigger(TriggerOptions, TriggerOptions);
  void SetNOTTrigger(TriggerOptions, TriggerOptions);
  void SetPhysicsSelection(Short_t PS){fPS = PS;}
  void SetZvxEvent(Double_t Zmin, Double_t Zmax) {fZvxEventMin = Zmin; fZvxEventMax = Zmax;}
  void SetMatchTriggerLevel(Short_t triggerlevel){fTrigLevel = triggerlevel;}
  void SetEtaSingleMuCut(Double_t etaMin, Double_t etaMax) {fEtaSingleMuMin = etaMin; fEtaSingleMuMax = etaMax;}
  void SetPtSingleMuCut(Double_t ptMin, Double_t ptMax) {fPtSingleMuMin = ptMin; fPtSingleMuMax = ptMax;}
  void SetActivitySelection(TString method){fActivitySelection = method;}
  void SetActivityMax(Double_t activityMax){fActivityMax = activityMax;}
  void SetCentralityBins(Int_t nBins, Double_t bins[]);
  void SetWriteTree(Bool_t writeTree){fWriteTree = writeTree;}
  void SetWriteHisto(Bool_t writeHisto){fWriteHisto = writeHisto;}
  
  Bool_t Init(Int_t);
  void Analize();
  void AnalizeMixing();

  void DoPairs   (const Double_t lcentrality);

private:
  void    Terminate();
  Short_t GetTriggerMask(TString);
  Bool_t  HasTriggerFired(AnalysisOption);
  Bool_t  HasPassedEventCuts(AliLMREvent*);
  Bool_t  HasPassedSingleMuCuts(AliLMRMuon*);
  Bool_t  SelectDimu(TLorentzVector *dimu);

  CutShapeType fCutShape;
  Int_t        fNEvents;
  TTree       *fEvtTree;
  AliLMREvent *fEvt;
  AliLMRMuon  *fMu1;
  AliLMRMuon  *fMu2;
  TH2D        *fUSMassC; //cross-check histo
  TH2D        *fUSMass;
  TH2D        *fUSMassLS;
  TH2D        *fUSMassPP;
  TH2D        *fUSMassMM;
  TH2D        *fUSMassMix;
  TH2D        *fUSMassMixLS;

  TH1D * hRabs1;
  TH1D * hPt1;
  TH1D * hpDCA1;
  TH1D * hChi21;
  TH1D * hChi2Match1;
  TH1D * hRabs2;
  TH1D * hPt2;
  TH1D * hpDCA2;
  TH1D * hChi22;
  TH1D * hChi2Match2;

  TH2D *hpDCA1vsMass0080;
  TH2D *hpDCA2vsMass0080;
  TH2D *hpDCA1vsMass0010;
  TH2D *hpDCA2vsMass0010;
  TH2D *hpDCA1vsMass1040;
  TH2D *hpDCA2vsMass1040;
  TH2D *hpDCA1vsMass4080;
  TH2D *hpDCA2vsMass4080;
  
  TH2D *hpT1vsMass0080;
  TH2D *hpT2vsMass0080;
  TH2D *hpT1vsMass0010;
  TH2D *hpT2vsMass0010;
  TH2D *hpT1vsMass1040;
  TH2D *hpT2vsMass1040;
  TH2D *hpT1vsMass4080;
  TH2D *hpT2vsMass4080;
  
  TH2D *hpT1vspT20080;
  TH2D *hpT1vspT20010;
  TH2D *hpT1vspT21040;
  TH2D *hpT1vspT24080;

  TH2D *hpDCAvsPt10080;
  TH2D *hpDCAvsPt20080;
  TH2D *hpDCAvsPt10010;
  TH2D *hpDCAvsPt20010;
  TH2D *hpDCAvsPt11040;
  TH2D *hpDCAvsPt21040;
  TH2D *hpDCAvsPt14080;
  TH2D *hpDCAvsPt24080;
 
  TH2D * hMassvsPt0080  ;
  TH2D * hMassvsPt0010  ;
  TH2D * hMassvsPt1040  ;
  TH2D * hMassvsPt4080  ;

  TH2D * hLSMassvsPt0080;
  TH2D * hLSMassvsPt0010;
  TH2D * hLSMassvsPt1040;
  TH2D * hLSMassvsPt4080;

  //new

  TH2D *hPhivsMass0080;
  TH2D *hThetavsMass0080;

  TH2D *hPhi1vsMass0080;
  TH2D *hEta1vsMass0080;
  TH2D *hPhi2vsMass0080;
  TH2D *hEta2vsMass0080;
 
  TH2D *hMassvsRapidity0080;

  //----------------

  Bool_t       fWriteTree;
  Bool_t       fWriteHisto;

  /* TTree       *ftreeSC; */
  /* TTree       *ftreeLS; */
  /* TTree       *ftreeMixSC; */
  /* TTree       *ftreeMixLS; */
  /* Double_t     tCentrality; */
  
  //---------------------------

  TTree * ftreePM02023;
  TTree * ftreePM02034;
  TTree * ftreePM02045;
  TTree * ftreePM02056;
  TTree * ftreePM02067;
  TTree * ftreePM02078;
  TTree * ftreePM02089;
  TTree * ftreePM02090;
  
  TTree * ftreePM205023;
  TTree * ftreePM205034;
  TTree * ftreePM205045;
  TTree * ftreePM205056;
  TTree * ftreePM205067;
  TTree * ftreePM205078;
  TTree * ftreePM205089;
  TTree * ftreePM205090;
  
  TTree * ftreePM509023;
  TTree * ftreePM509034;
  TTree * ftreePM509045;
  TTree * ftreePM509056;
  TTree * ftreePM509067;
  TTree * ftreePM509078;
  TTree * ftreePM509089;
  TTree * ftreePM509090;
  
  //----------------------

  TTree * ftreePP02023;
  TTree * ftreePP02034;
  TTree * ftreePP02045;
  TTree * ftreePP02056;
  TTree * ftreePP02067;
  TTree * ftreePP02078;
  TTree * ftreePP02089;
  TTree * ftreePP02090;
  
  TTree * ftreePP205023;
  TTree * ftreePP205034;
  TTree * ftreePP205045;
  TTree * ftreePP205056;
  TTree * ftreePP205067;
  TTree * ftreePP205078;
  TTree * ftreePP205089;
  TTree * ftreePP205090;
  
  TTree * ftreePP509023;
  TTree * ftreePP509034;
  TTree * ftreePP509045;
  TTree * ftreePP509056;
  TTree * ftreePP509067;
  TTree * ftreePP509078;
  TTree * ftreePP509089;
  TTree * ftreePP509090;
  
  //---------------------------

  TTree * ftreeMM02023;
  TTree * ftreeMM02034;
  TTree * ftreeMM02045;
  TTree * ftreeMM02056;
  TTree * ftreeMM02067;
  TTree * ftreeMM02078;
  TTree * ftreeMM02089;
  TTree * ftreeMM02090;
  
  TTree * ftreeMM205023;
  TTree * ftreeMM205034;
  TTree * ftreeMM205045;
  TTree * ftreeMM205056;
  TTree * ftreeMM205067;
  TTree * ftreeMM205078;
  TTree * ftreeMM205089;
  TTree * ftreeMM205090;
  
  TTree * ftreeMM509023;
  TTree * ftreeMM509034;
  TTree * ftreeMM509045;
  TTree * ftreeMM509056;
  TTree * ftreeMM509067;
  TTree * ftreeMM509078;
  TTree * ftreeMM509089;
  TTree * ftreeMM509090;

  //------------------------
  //MIX

  TTree * ftreeMixPM02023;
  TTree * ftreeMixPM02034;
  TTree * ftreeMixPM02045;
  TTree * ftreeMixPM02056;
  TTree * ftreeMixPM02067;
  TTree * ftreeMixPM02078;
  TTree * ftreeMixPM02089;
  TTree * ftreeMixPM02090;
  
  TTree * ftreeMixPM205023;
  TTree * ftreeMixPM205034;
  TTree * ftreeMixPM205045;
  TTree * ftreeMixPM205056;
  TTree * ftreeMixPM205067;
  TTree * ftreeMixPM205078;
  TTree * ftreeMixPM205089;
  TTree * ftreeMixPM205090;
  
  TTree * ftreeMixPM509023;
  TTree * ftreeMixPM509034;
  TTree * ftreeMixPM509045;
  TTree * ftreeMixPM509056;
  TTree * ftreeMixPM509067;
  TTree * ftreeMixPM509078;
  TTree * ftreeMixPM509089;
  TTree * ftreeMixPM509090;
 
  //----------------------

  TTree * ftreeMixPP02023;
  TTree * ftreeMixPP02034;
  TTree * ftreeMixPP02045;
  TTree * ftreeMixPP02056;
  TTree * ftreeMixPP02067;
  TTree * ftreeMixPP02078;
  TTree * ftreeMixPP02089;
  TTree * ftreeMixPP02090;
  
  TTree * ftreeMixPP205023;
  TTree * ftreeMixPP205034;
  TTree * ftreeMixPP205045;
  TTree * ftreeMixPP205056;
  TTree * ftreeMixPP205067;
  TTree * ftreeMixPP205078;
  TTree * ftreeMixPP205089;
  TTree * ftreeMixPP205090;
  
  TTree * ftreeMixPP509023;
  TTree * ftreeMixPP509034;
  TTree * ftreeMixPP509045;
  TTree * ftreeMixPP509056;
  TTree * ftreeMixPP509067;
  TTree * ftreeMixPP509078;
  TTree * ftreeMixPP509089;
  TTree * ftreeMixPP509090;

  //---------------------------

  TTree * ftreeMixMM02023;
  TTree * ftreeMixMM02034;
  TTree * ftreeMixMM02045;
  TTree * ftreeMixMM02056;
  TTree * ftreeMixMM02067;
  TTree * ftreeMixMM02078;
  TTree * ftreeMixMM02089;
  TTree * ftreeMixMM02090;
  
  TTree * ftreeMixMM205023;
  TTree * ftreeMixMM205034;
  TTree * ftreeMixMM205045;
  TTree * ftreeMixMM205056;
  TTree * ftreeMixMM205067;
  TTree * ftreeMixMM205078;
  TTree * ftreeMixMM205089;
  TTree * ftreeMixMM205090;
  
  TTree * ftreeMixMM509023;
  TTree * ftreeMixMM509034;
  TTree * ftreeMixMM509045;
  TTree * ftreeMixMM509056;
  TTree * ftreeMixMM509067;
  TTree * ftreeMixMM509078;
  TTree * ftreeMixMM509089;
  TTree * ftreeMixMM509090;

  //

  Double_t     tMass;
  Double_t tp1x;
  Double_t tp1y;
  Double_t tp1z;
  Double_t tp2x;
  Double_t tp2y;
  Double_t tp2z;

  TFile       *fOutput;
  TAxis       *fCentralityBins; // centrality bins of the analisys
  //  THnSparseI  *fMixingPool; // the key is the centrality bins
  TString      fActivitySelection;
  Int_t        fNCentralityBins;
  Double_t     fActivityMax;
  Short_t      fTriggerIncludedData;
  Short_t      fTriggerIncludedMixing;
  Short_t      fTriggerExcludedData;
  Short_t      fTriggerExcludedMixing;
  Bool_t       fIsTriggerSet;
  Short_t      fPS;
  Double_t     fZvxEventMin;
  Double_t     fZvxEventMax;
  Short_t      fTrigLevel;
  Double_t     fEtaSingleMuMin;
  Double_t     fEtaSingleMuMax;
  Double_t     fPtSingleMuMin;
  Double_t     fPtSingleMuMax;
  TH1D        *hZvertex;
  TH1D        *hCentrality;
  TH1D        *fHistMultiplicityOfMixedEvent;

  AliAnalysisMuMuEventCollection ***fEventColl; 
  AliAnalysisMuMuEvent *fEvent; 

  int fzVertexBins;
  int fnCentBins;
  int fnEventsToMix;
  int fMaxFirstMult;
  int fMaxSecondMult;

  ClassDef(AnalysisDMTaskDiMuons, 1);
};

#endif
