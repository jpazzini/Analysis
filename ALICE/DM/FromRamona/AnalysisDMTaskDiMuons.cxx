#include "AnalysisDMTaskDiMuons.h"
#include "AliAnalysisMuMuEventCollection.h"
#include "AliLog.h"
#include <iomanip>
#include <algorithm>

using namespace std;

ClassImp(AnalysisDMTaskDiMuons)


//__________________________________________________
AnalysisDMTaskDiMuons::AnalysisDMTaskDiMuons() :
  fCutShape(kRect), 
  fNEvents(0), 
  fEvtTree(0), 
  fEvt(0), 
  fEvent(0), 
  fMu1(0), 
  fMu2(0),
  fUSMass(0), 
  fUSMassC(0), 
  fUSMassLS(0), 
  fUSMassPP(0), 
  fUSMassMM(0), 
  fUSMassMix(0), 
  fUSMassMixLS(0), 
  hZvertex(0),
  hCentrality(0),
  fHistMultiplicityOfMixedEvent(0),

  hRabs1(0),
  hPt1(0),
  hpDCA1(0),
  hChi21(0),
  hChi2Match1(0),
  hRabs2(0),
  hPt2(0),
  hpDCA2(0),
  hChi22(0),
  hChi2Match2(0),

  hpDCA1vsMass0080(0),
  hpDCA2vsMass0080(0),
  hpDCA1vsMass0010(0),
  hpDCA2vsMass0010(0),
  hpDCA1vsMass1040(0),
  hpDCA2vsMass1040(0),
  hpDCA1vsMass4080(0),
  hpDCA2vsMass4080(0),
  
  hpT1vsMass0080(0),
  hpT2vsMass0080(0),
  hpT1vsMass0010(0),
  hpT2vsMass0010(0),
  hpT1vsMass1040(0),
  hpT2vsMass1040(0),
  hpT1vsMass4080(0),
  hpT2vsMass4080(0),
 
  hpT1vspT20080(0),
  hpT1vspT20010(0),
  hpT1vspT21040(0),
  hpT1vspT24080(0),

  hpDCAvsPt10010(0),
  hpDCAvsPt20010(0),
  hpDCAvsPt11040(0),
  hpDCAvsPt21040(0),
  hpDCAvsPt14080(0),
  hpDCAvsPt24080(0),
  
  hMassvsPt0080(0),
  hMassvsPt0010(0),
  hMassvsPt1040(0),
  hMassvsPt4080(0),
 
  hLSMassvsPt0010(0),
  hLSMassvsPt1040(0),
  hLSMassvsPt4080(0),

  hPhivsMass0080(0),
  hThetavsMass0080(0),
  hPhi1vsMass0080(0),
  hEta1vsMass0080(0),
  hPhi2vsMass0080(0),
  hEta2vsMass0080(0),
  hMassvsRapidity0080(0),

  fWriteTree(kFALSE),
  fWriteHisto(kFALSE),
  
//---------------------------
  
  ftreePM02023(0),
  ftreePM02034(0),
  ftreePM02045(0),
  ftreePM02056(0),
  ftreePM02067(0),
  ftreePM02078(0),
  ftreePM02089(0),
  ftreePM02090(0),
  
  ftreePM205023(0),
  ftreePM205034(0),
  ftreePM205045(0),
  ftreePM205056(0),
  ftreePM205067(0),
  ftreePM205078(0),
  ftreePM205089(0),
  ftreePM205090(0),
  
  ftreePM509023(0),
  ftreePM509034(0),
  ftreePM509045(0),
  ftreePM509056(0),
  ftreePM509067(0),
  ftreePM509078(0),
  ftreePM509089(0),
  ftreePM509090(0),
  
  //----------------------

  ftreePP02023(0),
  ftreePP02034(0),
  ftreePP02045(0),
  ftreePP02056(0),
  ftreePP02067(0),
  ftreePP02078(0),
  ftreePP02089(0),
  ftreePP02090(0),
  
  ftreePP205023(0),
  ftreePP205034(0),
  ftreePP205045(0),
  ftreePP205056(0),
  ftreePP205067(0),
  ftreePP205078(0),
  ftreePP205089(0),
  ftreePP205090(0),
  
  ftreePP509023(0),
  ftreePP509034(0),
  ftreePP509045(0),
  ftreePP509056(0),
  ftreePP509067(0),
  ftreePP509078(0),
  ftreePP509089(0),
  ftreePP509090(0),
  
  //---------------------------

  ftreeMM02023(0),
  ftreeMM02034(0),
  ftreeMM02045(0),
  ftreeMM02056(0),
  ftreeMM02067(0),
  ftreeMM02078(0),
  ftreeMM02089(0),
  ftreeMM02090(0),
  
  ftreeMM205023(0),
  ftreeMM205034(0),
  ftreeMM205045(0),
  ftreeMM205056(0),
  ftreeMM205067(0),
  ftreeMM205078(0),
  ftreeMM205089(0),
  ftreeMM205090(0),
  
  ftreeMM509023(0),
  ftreeMM509034(0),
  ftreeMM509045(0),
  ftreeMM509056(0),
  ftreeMM509067(0),
  ftreeMM509078(0),
  ftreeMM509089(0),
  ftreeMM509090(0),

//Mix
  
  ftreeMixPM02023(0),
  ftreeMixPM02034(0),
  ftreeMixPM02045(0),
  ftreeMixPM02056(0),
  ftreeMixPM02067(0),
  ftreeMixPM02078(0),
  ftreeMixPM02089(0),
  ftreeMixPM02090(0),
  
  ftreeMixPM205023(0),
  ftreeMixPM205034(0),
  ftreeMixPM205045(0),
  ftreeMixPM205056(0),
  ftreeMixPM205067(0),
  ftreeMixPM205078(0),
  ftreeMixPM205089(0),
  ftreeMixPM205090(0),
  
  ftreeMixPM509023(0),
  ftreeMixPM509034(0),
  ftreeMixPM509045(0),
  ftreeMixPM509056(0),
  ftreeMixPM509067(0),
  ftreeMixPM509078(0),
  ftreeMixPM509089(0),
  ftreeMixPM509090(0),
 
  //----------------------

  ftreeMixPP02023(0),
  ftreeMixPP02034(0),
  ftreeMixPP02045(0),
  ftreeMixPP02056(0),
  ftreeMixPP02067(0),
  ftreeMixPP02078(0),
  ftreeMixPP02089(0),
  ftreeMixPP02090(0),
  
  ftreeMixPP205023(0),
  ftreeMixPP205034(0),
  ftreeMixPP205045(0),
  ftreeMixPP205056(0),
  ftreeMixPP205067(0),
  ftreeMixPP205078(0),
  ftreeMixPP205089(0),
  ftreeMixPP205090(0),
  
  ftreeMixPP509023(0),
  ftreeMixPP509034(0),
  ftreeMixPP509045(0),
  ftreeMixPP509056(0),
  ftreeMixPP509067(0),
  ftreeMixPP509078(0),
  ftreeMixPP509089(0),
  ftreeMixPP509090(0),

  //---------------------------

  ftreeMixMM02023(0),
  ftreeMixMM02034(0),
  ftreeMixMM02045(0),
  ftreeMixMM02056(0),
  ftreeMixMM02067(0),
  ftreeMixMM02078(0),
  ftreeMixMM02089(0),
  ftreeMixMM02090(0),
  
  ftreeMixMM205023(0),
  ftreeMixMM205034(0),
  ftreeMixMM205045(0),
  ftreeMixMM205056(0),
  ftreeMixMM205067(0),
  ftreeMixMM205078(0),
  ftreeMixMM205089(0),
  ftreeMixMM205090(0),
  
  ftreeMixMM509023(0),
  ftreeMixMM509034(0),
  ftreeMixMM509045(0),
  ftreeMixMM509056(0),
  ftreeMixMM509067(0),
  ftreeMixMM509078(0),
  ftreeMixMM509089(0),
  ftreeMixMM509090(0),
  
//---------------
  
// ftreeLS(0),
// ftreeMixSC(0),
// ftreeMixLS(0),
//  tCentrality(0),
  tMass(0),
  tp1x(0),
  tp1y(0),
  tp1z(0),
  tp2x(0),
  tp2y(0),
  tp2z(0),
 
  fOutput(0), 
  fCentralityBins(0), 
  fTriggerIncludedData(0),
  fTriggerIncludedMixing(0), 
  fNCentralityBins(1), 
//  fMixingPool(0),
  fIsTriggerSet(false), 
  fPS(0),
  fZvxEventMin(0), 
  fZvxEventMax(0),
  fActivityMax(0), 
  fTrigLevel(0), 
  fActivitySelection(""),
  fEtaSingleMuMin(0), 
  fEtaSingleMuMax(0), 
  fPtSingleMuMin(0),
  fPtSingleMuMax(0),

  fEventColl(0x0),
  fzVertexBins(4), //(10)
  fnCentBins(20),
  fnEventsToMix(3),
//  fnEventsToMix(10),
  fMaxFirstMult(1000), // 1000 for protons 
  fMaxSecondMult(1000)  // was 100

{
  hZvertex  = new TH1D("hZvertex", "Z vertex", 600, -30, 30);
  //  hCentrality = new TH1D("hCentrality", "centrality", 102, -0.5,101.5);
  hCentrality = new TH1D("hCentrality", "centrality", 100, 0,100);
  fCentralityBins   = new TAxis(1, 0, 100);
  fHistMultiplicityOfMixedEvent = new TH1D("fHistMultiplicityOfMixedEvent", "Mix mult", 101, -0.5, 100.5);

  fUSMass      = new TH2D("fUSMass"     , "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);
  fUSMassC     = new TH2D("fUSMassC"    , "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);
  fUSMassLS    = new TH2D("fUSMassLS"   , "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);
  fUSMassPP    = new TH2D("fUSMassPP"   , "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);
  fUSMassMM    = new TH2D("fUSMassMM"   , "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);
  fUSMassMix   = new TH2D("fUSMassMix"  , "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);
  fUSMassMixLS = new TH2D("fUSMassMixLS", "OS #mu #mu mass vs cent ", 20000, 0, 20, 100, 0, 100);

  hRabs1      = new TH1D("hRabs1"     ,"hRabs1"     ,100,0,100);
  hPt1        = new TH1D("hPt1"       ,"hPt1"       ,1000,-100,100);
  hpDCA1      = new TH1D("hpDCA1"     ,"hpDCA1"     ,100,0,1000);
  hChi21      = new TH1D("hChi21"     ,"hChi21"     ,100,0,10);
  hChi2Match1 = new TH1D("hChi2Match1","hChi2Match1",100,0,10);

  hRabs2      = new TH1D("hRabs2"     ,"hRabs2"     ,100,0,100);
  hPt2        = new TH1D("hPt2"       ,"hPt2"       ,1000,0,100);
  hpDCA2      = new TH1D("hpDCA2"     ,"hpDCA2"     ,100,0,1000);
  hChi22      = new TH1D("hChi22"     ,"hChi22"     ,100,0,10);
  hChi2Match2 = new TH1D("hChi2Match2","hChi2Match2",100,0,10);

  hpDCA1vsMass0080 = new TH2D("hpDCA1vsMass0080"   ,"hpDCA1vsMass0080" ,1000,0,1000,400,0,20);
  hpDCA2vsMass0080 = new TH2D("hpDCA2vsMass0080"   ,"hpDCA2vsMass0080" ,1000,0,1000,400,0,20);
  hpDCA1vsMass0010 = new TH2D("hpDCA1vsMass0010"   ,"hpDCA1vsMass0010" ,1000,0,1000,400,0,20);
  hpDCA2vsMass0010 = new TH2D("hpDCA2vsMass0010"   ,"hpDCA2vsMass0010" ,1000,0,1000,400,0,20);
  hpDCA1vsMass1040 = new TH2D("hpDCA1vsMass1040"   ,"hpDCA1vsMass1040" ,1000,0,1000,400,0,20);
  hpDCA2vsMass1040 = new TH2D("hpDCA2vsMass1040"   ,"hpDCA2vsMass1040" ,1000,0,1000,400,0,20);
  hpDCA1vsMass4080 = new TH2D("hpDCA1vsMass4080"   ,"hpDCA1vsMass4080" ,1000,0,1000,400,0,20);
  hpDCA2vsMass4080 = new TH2D("hpDCA2vsMass4080"   ,"hpDCA2vsMass4080" ,1000,0,1000,400,0,20);
  
  hpT1vsMass0080   = new TH2D("hpT1vsMass0080"   ,"hpT1vsMass0080" ,1000,0,100,400,0,20);
  hpT2vsMass0080   = new TH2D("hpT2vsMass0080"   ,"hpT2vsMass0080" ,1000,0,100,400,0,20);
  hpT1vsMass0010   = new TH2D("hpT1vsMass0010"   ,"hpT1vsMass0010" ,1000,0,100,400,0,20);
  hpT2vsMass0010   = new TH2D("hpT2vsMass0010"   ,"hpT2vsMass0010" ,1000,0,100,400,0,20);
  hpT1vsMass1040   = new TH2D("hpT1vsMass1040"   ,"hpT1vsMass1040" ,1000,0,100,400,0,20);
  hpT2vsMass1040   = new TH2D("hpT2vsMass1040"   ,"hpT2vsMass1040" ,1000,0,100,400,0,20);
  hpT1vsMass4080   = new TH2D("hpT1vsMass4080"   ,"hpT1vsMass4080" ,1000,0,100,400,0,20);
  hpT2vsMass4080   = new TH2D("hpT2vsMass4080"   ,"hpT2vsMass4080" ,1000,0,100,400,0,20);
  
  hpT1vspT20080   = new TH2D("hpT1vspT20080"   ,"hpT1vspT20080" ,1000,0,100,1000,0,100);
  hpT1vspT20010   = new TH2D("hpT1vspT20010"   ,"hpT1vspT20010" ,1000,0,100,1000,0,100);
  hpT1vspT21040   = new TH2D("hpT1vspT21040"   ,"hpT1vspT21040" ,1000,0,100,1000,0,100);
  hpT1vspT24080   = new TH2D("hpT1vspT24080"   ,"hpT1vspT24080" ,1000,0,100,1000,0,100);

  hpDCAvsPt10080 = new TH2D("hpDCAvsPt10080"   ,"hpDCAvsPt10080" ,1000,0,1000,1000,0,100);
  hpDCAvsPt20080 = new TH2D("hpDCAvsPt20080"   ,"hpDCAvsPt20080" ,1000,0,1000,1000,0,100);
  hpDCAvsPt10010 = new TH2D("hpDCAvsPt10010"   ,"hpDCAvsPt10010" ,1000,0,1000,1000,0,100);
  hpDCAvsPt20010 = new TH2D("hpDCAvsPt20010"   ,"hpDCAvsPt20010" ,1000,0,1000,1000,0,100);
  hpDCAvsPt11040 = new TH2D("hpDCAvsPt11040"   ,"hpDCAvsPt11040" ,1000,0,1000,1000,0,100);
  hpDCAvsPt21040 = new TH2D("hpDCAvsPt21040"   ,"hpDCAvsPt21040" ,1000,0,1000,1000,0,100);
  hpDCAvsPt14080 = new TH2D("hpDCAvsPt14080"   ,"hpDCAvsPt14080" ,1000,0,1000,1000,0,100);
  hpDCAvsPt24080 = new TH2D("hpDCAvsPt24080"   ,"hpDCAvsPt24080" ,1000,0,1000,1000,0,100);

  hMassvsPt0080 = new TH2D("hMassvsPt0080"   ,"hMassvsPt0080" ,400,0,20,1000,0,100);
  hMassvsPt0010 = new TH2D("hMassvsPt0010"   ,"hMassvsPt0010" ,400,0,20,1000,0,100);
  hMassvsPt1040 = new TH2D("hMassvsPt1040"   ,"hMassvsPt1040" ,400,0,20,1000,0,100);
  hMassvsPt4080 = new TH2D("hMassvsPt4080"   ,"hMassvsPt4080" ,400,0,20,1000,0,100);

  hLSMassvsPt0080 = new TH2D("hLSMassvsPt0080"   ,"hLSMassvsPt0080" ,400,0,20,1000,0,100);
  hLSMassvsPt0010 = new TH2D("hLSMassvsPt0010"   ,"hLSMassvsPt0010" ,400,0,20,1000,0,100);
  hLSMassvsPt1040 = new TH2D("hLSMassvsPt1040"   ,"hLSMassvsPt1040" ,400,0,20,1000,0,100);
  hLSMassvsPt4080 = new TH2D("hLSMassvsPt4080"   ,"hLSMassvsPt4080" ,400,0,20,1000,0,100);


  hPhivsMass0080   = new TH2D("hPhivsMass0080"  ,"hPhivsMass0080"  ,100,-4, 4 ,400,0,20 );
  hThetavsMass0080 = new TH2D("hThetavsMass0080","hThetavsMass0080",100, 2, 4 ,400,0,20 );
  hPhi1vsMass0080  = new TH2D("hPhi1vsMass0080" ,"hPhi1vsMass0080" ,100,-4, 4 ,400,0,20 );
  hEta1vsMass0080  = new TH2D("hEta1vsMass0080" ,"hEta1vsMass0080" ,100,-5,-1 ,400,0,20 );
  hPhi2vsMass0080  = new TH2D("hPhi2vsMass0080" ,"hPhi2vsMass0080" ,100,-4, 4 ,400,0,20 );
  hEta2vsMass0080  = new TH2D("hEta2vsMass0080" ,"hEta2vsMass0080" ,100,-5,-1 ,400,0,20 );
  hMassvsRapidity0080  = new TH2D("hMassvsRapidity0080" ,"hMassvsRapidity0080" ,400,0,20 ,10000,-10,10 );
  //Event mixing
  fEventColl = new AliAnalysisMuMuEventCollection **[fzVertexBins]; 
  
  for (unsigned short i=0; i<fzVertexBins; i++) {
    fEventColl[i] = new AliAnalysisMuMuEventCollection *[fnCentBins];
    for (unsigned short j=0; j<fnCentBins; j++) {
      fEventColl[i][j] = new AliAnalysisMuMuEventCollection(fnEventsToMix+1, fMaxFirstMult, fMaxSecondMult);
    }
  }
  

}


//__________________________________________________
Bool_t AnalysisDMTaskDiMuons::Init(Int_t runNumber) {
  
  //TFile *inputFile = new TFile(Form("/Volumes/Data/LHC15o/AliLMREvents_%i.root", runNumber));

  TFile *inputFile = new TFile(Form("/gpfs/alice/alice/lea/DARKMATTER/Data/AliLMREvents.%i.root", runNumber)); //Pb--Pb 2015
  //  TFile *inputFile = new TFile(Form("/gpfs/alice/alice/lea/DARKMATTER/Data2018/%i/AliLMREvents.%i.root", runNumber, runNumber)); //Pb--Pb 2018
  //  TFile *inputFile = new TFile(Form("/gpfs/alice/alice/lea/DARKMATTER/DatapPb/%i/AliLMREvents.%i.root", runNumber,runNumber)); //p--Pb

  //  TFile *inputFile = new TFile("/home/ramona/Desktop/DiMuPbPb_New/Data/AliLMREvents.root");
  
  if (!inputFile->IsOpen()) {
    AliError(Form("Cannot open file %s", inputFile->GetName()));
    return kFALSE;
  }

  fEvtTree  = (TTree*) inputFile->Get("Data");

  if (!fEvtTree) {
    AliError("Cannot find the input tree");
    return kFALSE;
  }

  fNEvents = fEvtTree->GetEntries();
  AliInfo(Form("Input tree has %d events\n", fNEvents));

  fEvtTree->SetBranchAddress("fAliLMREvent", &fEvt);

  fOutput = new TFile(Form("Output/DM_%i.root", runNumber), "RECREATE");
  //  fOutput = new TFile(Form("Output2018/DM_%i.root", runNumber), "RECREATE");
  //  fOutput = new TFile(Form("OutputpPb/DM_%i.root", runNumber), "RECREATE");
  if (!fOutput) {
    AliError("Cannot create the output file");
    return kFALSE;
  }

  //------------------------------------------------------------

  ftreePM02023  = new TTree("ftreePM02023" ,"ftreePM02023" );
  ftreePM02034  = new TTree("ftreePM02034" ,"ftreePM02034" );
  ftreePM02045  = new TTree("ftreePM02045" ,"ftreePM02045" );
  ftreePM02056  = new TTree("ftreePM02056" ,"ftreePM02056" );
  ftreePM02067  = new TTree("ftreePM02067" ,"ftreePM02067" );
  ftreePM02078  = new TTree("ftreePM02078" ,"ftreePM02078" );
  ftreePM02089  = new TTree("ftreePM02089" ,"ftreePM02089" );
  ftreePM02090  = new TTree("ftreePM02090" ,"ftreePM02090" );

  ftreePM205023 = new TTree("ftreePM205023","ftreePM205023");
  ftreePM205034 = new TTree("ftreePM205034","ftreePM205034");
  ftreePM205045 = new TTree("ftreePM205045","ftreePM205045");
  ftreePM205056 = new TTree("ftreePM205056","ftreePM205056");
  ftreePM205067 = new TTree("ftreePM205067","ftreePM205067");
  ftreePM205078 = new TTree("ftreePM205078","ftreePM205078");
  ftreePM205089 = new TTree("ftreePM205089","ftreePM205089");
  ftreePM205090 = new TTree("ftreePM205090","ftreePM205090");

  ftreePM509023 = new TTree("ftreePM509023","ftreePM509023");
  ftreePM509034 = new TTree("ftreePM509034","ftreePM509034");
  ftreePM509045 = new TTree("ftreePM509045","ftreePM509045");
  ftreePM509056 = new TTree("ftreePM509056","ftreePM509056");
  ftreePM509067 = new TTree("ftreePM509067","ftreePM509067");
  ftreePM509078 = new TTree("ftreePM509078","ftreePM509078");
  ftreePM509089 = new TTree("ftreePM509089","ftreePM509089");
  ftreePM509090 = new TTree("ftreePM509090","ftreePM509090");

  ftreePP02023  = new TTree("ftreePP02023" ,"ftreePP02023" );
  ftreePP02034  = new TTree("ftreePP02034" ,"ftreePP02034" );
  ftreePP02045  = new TTree("ftreePP02045" ,"ftreePP02045" );
  ftreePP02056  = new TTree("ftreePP02056" ,"ftreePP02056" );
  ftreePP02067  = new TTree("ftreePP02067" ,"ftreePP02067" );
  ftreePP02078  = new TTree("ftreePP02078" ,"ftreePP02078" );
  ftreePP02089  = new TTree("ftreePP02089" ,"ftreePP02089" );
  ftreePP02090  = new TTree("ftreePP02090" ,"ftreePP02090" );
  
  ftreePP205023 = new TTree("ftreePP205023","ftreePP205023");
  ftreePP205034 = new TTree("ftreePP205034","ftreePP205034");
  ftreePP205045 = new TTree("ftreePP205045","ftreePP205045");
  ftreePP205056 = new TTree("ftreePP205056","ftreePP205056");
  ftreePP205067 = new TTree("ftreePP205067","ftreePP205067");
  ftreePP205078 = new TTree("ftreePP205078","ftreePP205078");
  ftreePP205089 = new TTree("ftreePP205089","ftreePP205089");
  ftreePP205090 = new TTree("ftreePP205090","ftreePP205090");
  
  ftreePP509023 = new TTree("ftreePP509023","ftreePP509023");
  ftreePP509034 = new TTree("ftreePP509034","ftreePP509034");
  ftreePP509045 = new TTree("ftreePP509045","ftreePP509045");
  ftreePP509056 = new TTree("ftreePP509056","ftreePP509056");
  ftreePP509067 = new TTree("ftreePP509067","ftreePP509067");
  ftreePP509078 = new TTree("ftreePP509078","ftreePP509078");
  ftreePP509089 = new TTree("ftreePP509089","ftreePP509089");
  ftreePP509090 = new TTree("ftreePP509090","ftreePP509090");
  
  ftreeMM02023  = new TTree("ftreeMM02023" ,"ftreeMM02023" );
  ftreeMM02034  = new TTree("ftreeMM02034" ,"ftreeMM02034" );
  ftreeMM02045  = new TTree("ftreeMM02045" ,"ftreeMM02045" );
  ftreeMM02056  = new TTree("ftreeMM02056" ,"ftreeMM02056" );
  ftreeMM02067  = new TTree("ftreeMM02067" ,"ftreeMM02067" );
  ftreeMM02078  = new TTree("ftreeMM02078" ,"ftreeMM02078" );
  ftreeMM02089  = new TTree("ftreeMM02089" ,"ftreeMM02089" );
  ftreeMM02090  = new TTree("ftreeMM02090" ,"ftreeMM02090" );
  
  ftreeMM205023 = new TTree("ftreeMM205023","ftreeMM205023");
  ftreeMM205034 = new TTree("ftreeMM205034","ftreeMM205034");
  ftreeMM205045 = new TTree("ftreeMM205045","ftreeMM205045");
  ftreeMM205056 = new TTree("ftreeMM205056","ftreeMM205056");
  ftreeMM205067 = new TTree("ftreeMM205067","ftreeMM205067");
  ftreeMM205078 = new TTree("ftreeMM205078","ftreeMM205078");
  ftreeMM205089 = new TTree("ftreeMM205089","ftreeMM205089");
  ftreeMM205090 = new TTree("ftreeMM205090","ftreeMM205090");
  
  ftreeMM509023 = new TTree("ftreeMM509023","ftreeMM509023");
  ftreeMM509034 = new TTree("ftreeMM509034","ftreeMM509034");
  ftreeMM509045 = new TTree("ftreeMM509045","ftreeMM509045");
  ftreeMM509056 = new TTree("ftreeMM509056","ftreeMM509056");
  ftreeMM509067 = new TTree("ftreeMM509067","ftreeMM509067");
  ftreeMM509078 = new TTree("ftreeMM509078","ftreeMM509078");
  ftreeMM509089 = new TTree("ftreeMM509089","ftreeMM509089");
  ftreeMM509090 = new TTree("ftreeMM509090","ftreeMM509090");

  //--------------------
  //Mix
 
  ftreeMixPM02023  = new TTree("ftreeMixPM02023" ,"ftreeMixPM02023" );
  ftreeMixPM02034  = new TTree("ftreeMixPM02034" ,"ftreeMixPM02034" );
  ftreeMixPM02045  = new TTree("ftreeMixPM02045" ,"ftreeMixPM02045" );
  ftreeMixPM02056  = new TTree("ftreeMixPM02056" ,"ftreeMixPM02056" );
  ftreeMixPM02067  = new TTree("ftreeMixPM02067" ,"ftreeMixPM02067" );
  ftreeMixPM02078  = new TTree("ftreeMixPM02078" ,"ftreeMixPM02078" );
  ftreeMixPM02089  = new TTree("ftreeMixPM02089" ,"ftreeMixPM02089" );
  ftreeMixPM02090  = new TTree("ftreeMixPM02090" ,"ftreeMixPM02090" );

  ftreeMixPM205023 = new TTree("ftreeMixPM205023","ftreeMixPM205023");
  ftreeMixPM205034 = new TTree("ftreeMixPM205034","ftreeMixPM205034");
  ftreeMixPM205045 = new TTree("ftreeMixPM205045","ftreeMixPM205045");
  ftreeMixPM205056 = new TTree("ftreeMixPM205056","ftreeMixPM205056");
  ftreeMixPM205067 = new TTree("ftreeMixPM205067","ftreeMixPM205067");
  ftreeMixPM205078 = new TTree("ftreeMixPM205078","ftreeMixPM205078");
  ftreeMixPM205089 = new TTree("ftreeMixPM205089","ftreeMixPM205089");
  ftreeMixPM205090 = new TTree("ftreeMixPM205090","ftreeMixPM205090");

  ftreeMixPM509023 = new TTree("ftreeMixPM509023","ftreeMixPM509023");
  ftreeMixPM509034 = new TTree("ftreeMixPM509034","ftreeMixPM509034");
  ftreeMixPM509045 = new TTree("ftreeMixPM509045","ftreeMixPM509045");
  ftreeMixPM509056 = new TTree("ftreeMixPM509056","ftreeMixPM509056");
  ftreeMixPM509067 = new TTree("ftreeMixPM509067","ftreeMixPM509067");
  ftreeMixPM509078 = new TTree("ftreeMixPM509078","ftreeMixPM509078");
  ftreeMixPM509089 = new TTree("ftreeMixPM509089","ftreeMixPM509089");
  ftreeMixPM509090 = new TTree("ftreeMixPM509090","ftreeMixPM509090");

  ftreeMixPP02023  = new TTree("ftreeMixPP02023" ,"ftreeMixPP02023" );
  ftreeMixPP02034  = new TTree("ftreeMixPP02034" ,"ftreeMixPP02034" );
  ftreeMixPP02045  = new TTree("ftreeMixPP02045" ,"ftreeMixPP02045" );
  ftreeMixPP02056  = new TTree("ftreeMixPP02056" ,"ftreeMixPP02056" );
  ftreeMixPP02067  = new TTree("ftreeMixPP02067" ,"ftreeMixPP02067" );
  ftreeMixPP02078  = new TTree("ftreeMixPP02078" ,"ftreeMixPP02078" );
  ftreeMixPP02089  = new TTree("ftreeMixPP02089" ,"ftreeMixPP02089" );
  ftreeMixPP02090  = new TTree("ftreeMixPP02090" ,"ftreeMixPP02090" );
  
  ftreeMixPP205023 = new TTree("ftreeMixPP205023","ftreeMixPP205023");
  ftreeMixPP205034 = new TTree("ftreeMixPP205034","ftreeMixPP205034");
  ftreeMixPP205045 = new TTree("ftreeMixPP205045","ftreeMixPP205045");
  ftreeMixPP205056 = new TTree("ftreeMixPP205056","ftreeMixPP205056");
  ftreeMixPP205067 = new TTree("ftreeMixPP205067","ftreeMixPP205067");
  ftreeMixPP205078 = new TTree("ftreeMixPP205078","ftreeMixPP205078");
  ftreeMixPP205089 = new TTree("ftreeMixPP205089","ftreeMixPP205089");
  ftreeMixPP205090 = new TTree("ftreeMixPP205090","ftreeMixPP205090");
  
  ftreeMixPP509023 = new TTree("ftreeMixPP509023","ftreeMixPP509023");
  ftreeMixPP509034 = new TTree("ftreeMixPP509034","ftreeMixPP509034");
  ftreeMixPP509045 = new TTree("ftreeMixPP509045","ftreeMixPP509045");
  ftreeMixPP509056 = new TTree("ftreeMixPP509056","ftreeMixPP509056");
  ftreeMixPP509067 = new TTree("ftreeMixPP509067","ftreeMixPP509067");
  ftreeMixPP509078 = new TTree("ftreeMixPP509078","ftreeMixPP509078");
  ftreeMixPP509089 = new TTree("ftreeMixPP509089","ftreeMixPP509089");
  ftreeMixPP509090 = new TTree("ftreeMixPP509090","ftreeMixPP509090");
  
  ftreeMixMM02023  = new TTree("ftreeMixMM02023" ,"ftreeMixMM02023" );
  ftreeMixMM02034  = new TTree("ftreeMixMM02034" ,"ftreeMixMM02034" );
  ftreeMixMM02045  = new TTree("ftreeMixMM02045" ,"ftreeMixMM02045" );
  ftreeMixMM02056  = new TTree("ftreeMixMM02056" ,"ftreeMixMM02056" );
  ftreeMixMM02067  = new TTree("ftreeMixMM02067" ,"ftreeMixMM02067" );
  ftreeMixMM02078  = new TTree("ftreeMixMM02078" ,"ftreeMixMM02078" );
  ftreeMixMM02089  = new TTree("ftreeMixMM02089" ,"ftreeMixMM02089" );
  ftreeMixMM02090  = new TTree("ftreeMixMM02090" ,"ftreeMixMM02090" );
  
  ftreeMixMM205023 = new TTree("ftreeMixMM205023","ftreeMixMM205023");
  ftreeMixMM205034 = new TTree("ftreeMixMM205034","ftreeMixMM205034");
  ftreeMixMM205045 = new TTree("ftreeMixMM205045","ftreeMixMM205045");
  ftreeMixMM205056 = new TTree("ftreeMixMM205056","ftreeMixMM205056");
  ftreeMixMM205067 = new TTree("ftreeMixMM205067","ftreeMixMM205067");
  ftreeMixMM205078 = new TTree("ftreeMixMM205078","ftreeMixMM205078");
  ftreeMixMM205089 = new TTree("ftreeMixMM205089","ftreeMixMM205089");
  ftreeMixMM205090 = new TTree("ftreeMixMM205090","ftreeMixMM205090");
  
  ftreeMixMM509023 = new TTree("ftreeMixMM509023","ftreeMixMM509023");
  ftreeMixMM509034 = new TTree("ftreeMixMM509034","ftreeMixMM509034");
  ftreeMixMM509045 = new TTree("ftreeMixMM509045","ftreeMixMM509045");
  ftreeMixMM509056 = new TTree("ftreeMixMM509056","ftreeMixMM509056");
  ftreeMixMM509067 = new TTree("ftreeMixMM509067","ftreeMixMM509067");
  ftreeMixMM509078 = new TTree("ftreeMixMM509078","ftreeMixMM509078");
  ftreeMixMM509089 = new TTree("ftreeMixMM509089","ftreeMixMM509089");
  ftreeMixMM509090 = new TTree("ftreeMixMM509090","ftreeMixMM509090");

  //----------
  // TREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

  ftreePM02023  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02034  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02045  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02056  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02067  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02078  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02089  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM02090  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
		
  ftreePM205023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM205090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
		
  ftreePM509023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePM509090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 

  //additional info --------------------------------------------------------

  ftreePM02023  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02034  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02045  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02056  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02067  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02078  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02089  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM02090  ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 

  ftreePM02023  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02034  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02045  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02056  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02067  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02078  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02089  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM02090  ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 

  ftreePM02023  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02034  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02045  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02056  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02067  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02078  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02089  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM02090  ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 

  //2
  ftreePM02023  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02034  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02045  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02056  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02067  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02078  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02089  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM02090  ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 

  ftreePM02023  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02034  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02045  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02056  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02067  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02078  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02089  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM02090  ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 

  ftreePM02023  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02034  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02045  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02056  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02067  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02078  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02089  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM02090  ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 

  
  //-----------
  
  ftreePM205023 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205034 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205045 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205056 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205067 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205078 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205089 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM205090 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 

  ftreePM205023 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205034 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205045 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205056 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205067 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205078 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205089 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM205090 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 

  ftreePM205023 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205034 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205045 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205056 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205067 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205078 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205089 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM205090 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 

  //---
  
  ftreePM205023 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205034 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205045 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205056 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205067 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205078 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205089 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM205090 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 

  ftreePM205023 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205034 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205045 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205056 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205067 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205078 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205089 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM205090 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 

  ftreePM205023 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205034 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205045 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205056 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205067 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205078 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205089 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM205090 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 


  //--------------
		
  ftreePM509023 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509034 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509045 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509056 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509067 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509078 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509089 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
  ftreePM509090 ->Branch("tp1x"      ,&tp1x      ,"tp1x/D"    ); 
		
  ftreePM509023 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509034 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509045 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509056 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509067 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509078 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509089 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 
  ftreePM509090 ->Branch("tp1y"      ,&tp1y      ,"tp1y/D"    ); 

		
  ftreePM509023 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509034 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509045 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509056 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509067 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509078 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509089 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 
  ftreePM509090 ->Branch("tp1z"      ,&tp1z      ,"tp1z/D"    ); 

  //2
		
  ftreePM509023 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509034 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509045 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509056 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509067 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509078 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509089 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
  ftreePM509090 ->Branch("tp2x"      ,&tp2x      ,"tp2x/D"    ); 
		
  ftreePM509023 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509034 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509045 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509056 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509067 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509078 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509089 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 
  ftreePM509090 ->Branch("tp2y"      ,&tp2y      ,"tp2y/D"    ); 

		
  ftreePM509023 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509034 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509045 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509056 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509067 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509078 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509089 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 
  ftreePM509090 ->Branch("tp2z"      ,&tp2z      ,"tp2z/D"    ); 


  //--------------------------------------------------------------------
		
  ftreePP02023  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02034  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02045  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02056  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02067  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02078  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02089  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP02090  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		
  ftreePP205023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP205090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		
  ftreePP509023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreePP509090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		
  ftreeMM02023  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02034  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02045  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02056  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02067  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02078  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02089  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM02090  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		
  ftreeMM205023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM205090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		
  ftreeMM509023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMM509090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 

  //--------------------
  //Mix
 
  ftreeMixPM02023  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02034  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02045  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02056  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02067  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02078  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02089  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM02090  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
		   
  ftreeMixPM205023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM205090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
		   
  ftreeMixPM509023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPM509090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
		   
  ftreeMixPP02023  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02034  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02045  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02056  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02067  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02078  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02089  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP02090  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		   
  ftreeMixPP205023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP205090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		   
  ftreeMixPP509023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixPP509090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		   
  ftreeMixMM02023  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02034  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02045  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02056  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02067  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02078  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02089  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM02090  ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		   
  ftreeMixMM205023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM205090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  		   
  ftreeMixMM509023 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509034 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509045 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509056 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509067 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509078 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509089 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 
  ftreeMixMM509090 ->Branch("tMass"      ,&tMass      ,"tMass/D"    ); 

  //------

  // // 1: number of events - 2: centrality - 3: index on number of muons in the event
  // Int_t mixingPoolBins[3]    = {    fNEvents,    1};
  // Double_t mixingPoolMins[3] = {        -0.5,   0.};
  // Double_t mixingPoolMaxs[3] = {fNEvents-0.5, 100.};
  // fMixingPool = new THnSparseI("fMixingPool", "fMixingPool", 2, mixingPoolBins, mixingPoolMins, mixingPoolMaxs);

  return kTRUE;
}


//__________________________________________________
AnalysisDMTaskDiMuons::~AnalysisDMTaskDiMuons() {
  delete fEvtTree;
  delete fEvt;
  delete fEvent;
  delete fMu1;
  delete fMu2;
  delete fUSMass;
  delete fUSMassC;
  delete fUSMassLS;
  delete fUSMassMM;
  delete fUSMassPP;
  delete fUSMassMix;
  delete fUSMassMixLS;
  delete fOutput;
  delete fCentralityBins;
  //  delete fMixingPool;

  for (unsigned short i=0; i<fzVertexBins; i++) {
    for (unsigned short j=0; j<fnCentBins; j++) {
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }
  
}


//__________________________________________________
void AnalysisDMTaskDiMuons::Analize() {
  auto nAnalizedEvents = 0;
 
  int centralityBin=0;
  
  for (auto ievt = 0; ievt < fNEvents; ++ievt) {
    
    //  for (auto ievt = 9; ievt < 10; ++ievt) {
    // cout<<"\n\n\n\n\n\n\nEvent:   "<<ievt<<endl;
    
    fEvtTree->GetEntry(ievt);
    
    // Trigger selection
    if (!(HasTriggerFired(kDirect))) continue;

    // Event cuts
    if (!HasPassedEventCuts(fEvt)) continue;
  
    Int_t nMuons = fEvt->GetNMuons();
    if (nMuons < 2) continue;
    
    Double_t centrality = fEvt->GetMultiplicity("V0M")          ;//TMath::Min(fEvt->GetMultiplicity(fActivitySelection), 99.9f);
    if(centrality>100)continue;
    
    hCentrality->Fill(centrality);
    
    int fCount = 0;
    int sCount = 0;
  
    bool isPositive = kFALSE;
    
    int zBin=0;
    Double_t zvtx = fEvt->GetZVertex();
    hZvertex->Fill(zvtx);
     
    
    double zStep=2*10/double(fzVertexBins), zStart=-10.;
    
    for (int i=0; i<fzVertexBins; i++) {
      if ((zvtx > zStart+i*zStep) && (zvtx < zStart+(i+1)*zStep)) {
	zBin=i;
	break;
      }
    }
     
    //---------------------------------------------------------
    //CENTRALITY
    // ... and centrality    //FIXME : find out how centrality wokrs in AOD in pp and pPb
    if(centrality < 5.) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
    else if(centrality < 10.) centralityBin=18;
    else if(centrality < 15.) centralityBin=17;
    else if(centrality < 20.) centralityBin=16;
    else if(centrality < 25.) centralityBin=15;
    else if(centrality < 30.) centralityBin=14;
    else if(centrality < 35.) centralityBin=13;
    else if(centrality < 40.) centralityBin=12;
    else if(centrality < 45.) centralityBin=11;
    else if(centrality < 50.) centralityBin=10;
    else if(centrality < 55.) centralityBin=9; 
    else if(centrality < 60.) centralityBin=8;
    else if(centrality < 65.) centralityBin=7;
    else if(centrality < 70.) centralityBin=6;
    else if(centrality < 75.) centralityBin=5;
    else if(centrality < 80.) centralityBin=4;
    else if(centrality < 85.) centralityBin=3;
    else if(centrality < 90.) centralityBin=2;   // from here on, wont be filled because the range selected in AddTask is 0-90
    else if(centrality < 95.) centralityBin=1;
    else if(centrality <= 100.) centralityBin=0;

    //centralityBin=0;
   
    fEventColl[zBin][centralityBin]->FifoShift();
    fEvent = fEventColl[zBin][centralityBin]->fEvt;
    
    // cout<<"centrality: "<<centrality<<" lbin: "<<centralityBin<<endl;

    //----------------------------------------------------------
    
    // Int_t centralityBin = fCentralityBins->FindBin(centrality)-1;
    // centralityBin = fCentralityBins->FindBin(centrality)-1;
    //    if (centralityBin < 0) continue;
    
    for (Int_t iMu = 0; iMu < nMuons; ++iMu) {
      
      //    // cout<<"i: "<<iMu<<endl;

      fMu1 = (AliLMRMuon*) fEvt->GetMuon(iMu);
      if (!HasPassedSingleMuCuts(fMu1)) continue;

      if (!(fMu1->GetRabs()>= 17.6  && fMu1->GetRabs() <= 89.5))continue;
      if (fMu1->GetChi2()>4) continue; 
      if (fMu1->GetpDCA()>500) continue; 

      double pt1 = TMath::Sqrt(fMu1->Px()*fMu1->Px()+fMu1->Py()*fMu1->Py());

      hRabs1->Fill(fMu1->GetRabs());
      hPt1->Fill(fMu1->Charge()*pt1);
      hpDCA1->Fill(fMu1->GetpDCA());
      hChi21->Fill(fMu1->GetChi2());
      hChi2Match1->Fill(fMu1->GetChi2Match());

      TLorentzVector vMu1 = fMu1->P4();

      fEvent->fReconstructedFirst[fCount].fVectMu = vMu1;
      // if(fMu1->Charge() > 0) isPositive = kTRUE;
      fEvent->fReconstructedFirst[fCount].charge = fMu1->Charge();

      fCount++;
      
      for (Int_t jMu = iMu+1; jMu < nMuons; ++jMu) {
	
	fMu2 = (AliLMRMuon*) fEvt->GetMuon(jMu);
      
	if (!HasPassedSingleMuCuts(fMu2))  continue;

	if (!(fMu2->GetRabs()  >= 17.6  && fMu2->GetRabs() <= 89.5))continue;
	if (fMu2->GetChi2()>4) continue; 
	if (fMu2->GetpDCA()>500) continue; 

	double pt2 = TMath::Sqrt(fMu2->Px()*fMu2->Px()+fMu2->Py()*fMu2->Py());

	//	if(pt2 < (-pt1+2))continue; //pt1pt2cut

	hRabs2->Fill(fMu2->GetRabs());
	hPt2->Fill(pt2);
	hpDCA2->Fill(fMu2->GetpDCA());
	hChi22->Fill(fMu2->GetChi2());
	hChi2Match2->Fill(fMu2->GetChi2Match());

        TLorentzVector vMu2 = fMu2->P4();
	
	/*
	  fEvent->fReconstructedSecond[sCount].sVectMu = vMu2;
	  if(fMu2->Charge() > 0) isPositive = kTRUE;
	  fEvent->fReconstructedSecond[sCount].isPos = isPositive;
	  sCount++;
	*/
	
	// opposite sign
	// if (fMu1->Charge() == fMu2->Charge()) continue;
	
	TLorentzVector dimu = vMu1+vMu2;
	if(!SelectDimu(&dimu))continue;
	
	if (fMu1->Charge() == fMu2->Charge()){ //LIKE-SIGN
	  fUSMassLS->Fill(dimu.M(),centrality);
	  
	  if(centrality>=0 && centrality<90){ 
	    hMassvsRapidity0080 ->Fill(dimu.M(),dimu.Rapidity());
	  }

	  if(fMu1->Charge()>0){
	    fUSMassPP->Fill(dimu.M(),centrality);
	  }
	  
	  else
	    fUSMassMM->Fill(dimu.M(),centrality);
	  
	}
	// // cout<<"second count: "<<sCount<<endl;
	else{ // opposite size

	  fUSMassC->Fill(dimu.M(),centrality);
	  
	  if(centrality>=0 && centrality<10){
	    hpDCA1vsMass0010 ->Fill(fMu1->GetpDCA(),dimu.M());
	    hpDCA2vsMass0010 ->Fill(fMu2->GetpDCA(),dimu.M());
	    hpT1vsMass0010   ->Fill(pt1,dimu.M());
	    hpT2vsMass0010   ->Fill(pt2,dimu.M());
	    hpDCAvsPt10010   ->Fill(fMu1->GetpDCA(),pt1);
	    hpDCAvsPt20010   ->Fill(fMu2->GetpDCA(),pt2);    
	    
	    hpT1vspT20010    ->Fill(pt1,pt2); 
	    hMassvsPt0010    ->Fill(dimu.M(),dimu.Pt());
	  }

	  if(centrality>=10 && centrality<40){

	    hpDCA1vsMass1040 ->Fill(fMu1->GetpDCA(),dimu.M());	  
	    hpDCA2vsMass1040 ->Fill(fMu2->GetpDCA(),dimu.M());	  
	    hpT1vsMass1040   ->Fill(pt1,dimu.M());		  
	    hpT2vsMass1040   ->Fill(pt2,dimu.M());		  
	    hpDCAvsPt11040   ->Fill(fMu1->GetpDCA(),pt1);	  
	    hpDCAvsPt21040   ->Fill(fMu2->GetpDCA(),pt2);   
	   
	    hpT1vspT21040    ->Fill(pt1,pt2); 
 
	    hMassvsPt1040    ->Fill(dimu.M(),dimu.Pt());       
	  }
	  
	  if(centrality>=40 && centrality<90){
	    
	    hpDCA1vsMass4080 ->Fill(fMu1->GetpDCA(),dimu.M());	  
	    hpDCA2vsMass4080 ->Fill(fMu2->GetpDCA(),dimu.M());	  
	    hpT1vsMass4080   ->Fill(pt1,dimu.M());		  
	    hpT2vsMass4080   ->Fill(pt2,dimu.M());		  
	    hpDCAvsPt14080   ->Fill(fMu1->GetpDCA(),pt1);	  
	    hpDCAvsPt24080   ->Fill(fMu2->GetpDCA(),pt2);      
	   
	    hpT1vspT24080  ->Fill(pt1,pt2); 
	    
	    hMassvsPt4080    ->Fill(dimu.M(),dimu.Pt());     

	  }

	  if(centrality>=0 && centrality<90){
	    
	    hpDCA1vsMass0080 ->Fill(fMu1->GetpDCA(),dimu.M());	  
	    hpDCA2vsMass0080 ->Fill(fMu2->GetpDCA(),dimu.M());	  
	    hpT1vsMass0080   ->Fill(pt1,dimu.M());		  
	    hpT2vsMass0080   ->Fill(pt2,dimu.M());		  
	    hpDCAvsPt10080   ->Fill(fMu1->GetpDCA(),pt1);	  
	    hpDCAvsPt20080   ->Fill(fMu2->GetpDCA(),pt2);      
	   
	    hpT1vspT20080  ->Fill(pt1,pt2); 
	    
	    hMassvsPt0080    ->Fill(dimu.M(),dimu.Pt());     

	    
	    hPhivsMass0080  ->Fill(dimu.Phi()  ,dimu.M());
	    hThetavsMass0080->Fill(dimu.Theta(),dimu.M());
	    hPhi1vsMass0080 ->Fill(vMu1.Phi(),dimu.M());
	    hEta1vsMass0080 ->Fill(vMu2.Eta(),dimu.M());
	    hPhi2vsMass0080 ->Fill(vMu1.Phi(),dimu.M());
	    hEta2vsMass0080 ->Fill(vMu2.Eta(),dimu.M());

	    
	  }

	}
	  // cout<<"---------> Massa 1: "<<(vMu1+vMu2).M()<<endl;
	
      }
    }
    
    fEvent->fNumberCandidateFirst  = fCount;
    //    fEvent->fNumberCandidateSecond = sCount;
    // // cout<<"candiate first:  "<<fCount<<" "<<fEvent->fNumberCandidateFirst<<endl;
    // // cout<<"candiate second: "<<sCount<<" "<<fEvent->fNumberCandidateSecond<<endl;
    // cout<<"---------------------------------------\nDo pairs"<<endl;
    
    DoPairs(centrality); 
    nAnalizedEvents++;
    
  }//Event loop
   
  // // cout<<"Do pairs??"<<endl;

 
  cout << "# good events: " << nAnalizedEvents << " -- "<<fNEvents<< endl;
  Terminate();
}

//__________________________________________________
void AnalysisDMTaskDiMuons::DoPairs( const Double_t centrality) {

  TLorentzVector vect1;
  TLorentzVector vect2;
  // bool isPos1 = kFALSE;
  // bool isPos2 = kFALSE;
  Short_t        charge1;
  Short_t        charge2;
  double mass = -999.;
  double pt   = -999.;
  
  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
  
  // // cout<<"centrality "<<lcentrality<<endl;
  // cout<<"-------------------------------------------------------------"<<endl;
  // cout<<"candidates 1: "<<fEvent->fNumberCandidateFirst<<endl;
  // cout<<"candidates 2: "<<fEvent->fNumberCandidateSecond<<endl;
    
  for (int i=0; i<fEvent->fNumberCandidateFirst; ++i) {
    
    vect1   = fEvent->fReconstructedFirst[i].fVectMu;
    charge1 = fEvent->fReconstructedFirst[i].charge;
   
    // cout<<i<<" Mu1 pT: "<<vect1.Pt()<<endl;

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; ++eventNumber) {  
      // For same event pairs
      if (!multmixedcounted && eventNumber!=0 && ((fEvent+eventNumber)->fNumberCandidateFirst)!=0.) evmultmixed++;
      
      //for (int j=0; j<(fEvent+eventNumber)->fNumberCandidateSecond; j++) {
      for (int j=i+1; j<(fEvent+eventNumber)->fNumberCandidateFirst; ++j) {
	//for (Int_t jMu = iMu+1; jMu < nMuons; ++jMu) {  

	// vect2  = (fEvent+eventNumber)->fReconstructedSecond[j].sVectMu;
	// isPos2 = (fEvent+eventNumber)->fReconstructedSecond[j].isPos;

	vect2   = (fEvent+eventNumber)->fReconstructedFirst[j].fVectMu;
	//	isPos2 = (fEvent+eventNumber)->fReconstructedFirst[j].isPos;
	charge2 = (fEvent+eventNumber)->fReconstructedFirst[j].charge;
	// cout<<j<<" Mu2 pT: "<<vect2.Pt()<<endl;

	//if(vect2.M() == 0)continue;
	
	//if (charge1 == charge2) continue;

	TLorentzVector dimu = vect1+vect2;
	if(!SelectDimu(&dimu))continue;
	
	//	mass = ((vect1+vect2).M());
	mass = (dimu.M());
	pt   = (dimu.Pt());

	//	tCentrality      = centrality;
	tMass            = mass;
	
	tp1x = vect1.Px();
	tp1y = vect1.Py();
	tp1z = vect1.Pz();

	tp2x = vect2.Px();
	tp2y = vect2.Py();
	tp2z = vect2.Pz();
	

	//	cout<<isPos1<<" "<<isPos2<<" ls: "<<isPos1*isPos2<<endl;
	
	//	if((isPos1 && isPos2) || (!isPos1 && !isPos2))continue; //exclude ls pair
	// if((isPos1  && isPos2))continue;
	// if((!isPos1 && !isPos2))continue; //exclude ls pair
	// cout<<isPos1<<" "<<isPos2<<endl;
	//	if((isPos1 * isPos2)>0)continue; //exclude ls pair
	//	// cout<<"---------> pt1 dentro : "<<vect1.Pt()<<endl;

	//// cout<<"---------> Massa 2: "<<mass<<endl;
	
	if (eventNumber==0) {//Same event pair histogramming
	  //  // cout<<"vect 2 i: "<<vect2.Pt()<<endl;
	  // cout<<"---------> Massa 2: "<<mass<<endl;

	  //	  // cout<<"Massa otside: "<<mass<<endl;
	  
	  if (charge1 != charge2){ 
	    fUSMass->Fill(mass,centrality);
	    //	    ftreeSC->Fill();
	    if(fWriteTree == kTRUE){
	      if(centrality>=0 && centrality<20){
		if(pt>=2 && pt < 3)
		  ftreePM02023 ->Fill();
		if(pt>=3 && pt < 4)
		  ftreePM02034 ->Fill();
		if(pt>=4 && pt < 5)
		  ftreePM02045 ->Fill();
		if(pt>=5 && pt < 6)
		  ftreePM02056 ->Fill();
		if(pt>=6 && pt < 7)
		  ftreePM02067 ->Fill();
		if(pt>=7 && pt < 8)
		  ftreePM02078 ->Fill();
		if(pt>=8 && pt < 9)
		  ftreePM02089 ->Fill();
		if(pt>=9 )
		  ftreePM02090 ->Fill();
	      }

	      if(centrality>=20 && centrality<50){
		
		if(pt>=2 && pt < 3)
		  ftreePM205023 ->Fill();
		if(pt>=3 && pt < 4)
		  ftreePM205034 ->Fill();
		if(pt>=4 && pt < 5)
		  ftreePM205045 ->Fill();
		if(pt>=5 && pt < 6)
		  ftreePM205056 ->Fill();
		if(pt>=6 && pt < 7)
		  ftreePM205067 ->Fill();
		if(pt>=7 && pt < 8)
		  ftreePM205078 ->Fill();
		if(pt>=8 && pt < 9)
		  ftreePM205089 ->Fill();
		if(pt>=9 )
		  ftreePM205090 ->Fill();
	      }
	      
	      
	      if(centrality>=50 && centrality<90){
		if(pt>=2 && pt < 3)
		  ftreePM509023 ->Fill();
		if(pt>=3 && pt < 4)
		  ftreePM509034 ->Fill();
		if(pt>=4 && pt < 5)
		  ftreePM509045 ->Fill();
		if(pt>=5 && pt < 6)
		  ftreePM509056 ->Fill();
		if(pt>=6 && pt < 7)
		  ftreePM509067 ->Fill();
		if(pt>=7 && pt < 8)
		  ftreePM509078 ->Fill();
		if(pt>=8 && pt < 9)
		  ftreePM509089 ->Fill();
		if(pt>=9 )
		  ftreePM509090 ->Fill();
	      }
	      
	    }
	  }
	  
	  else{ //like-sign
	    
	    if (charge1 > 0){ 
	      
	      if(fWriteTree == kTRUE){
		if(centrality>=0 && centrality<20){
		  if(pt>=2 && pt < 3)
		    ftreePP02023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreePP02034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreePP02045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreePP02056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreePP02067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreePP02078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreePP02089 ->Fill();
		  if(pt>=9 )
		    ftreePP02090 ->Fill();
		}

		if(centrality>=20 && centrality<50){
		
		  if(pt>=2 && pt < 3)
		    ftreePP205023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreePP205034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreePP205045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreePP205056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreePP205067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreePP205078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreePP205089 ->Fill();
		  if(pt>=9 )
		    ftreePP205090 ->Fill();
		}
	      
	      
		if(centrality>=50 && centrality<90){
		  if(pt>=2 && pt < 3)
		    ftreePP509023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreePP509034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreePP509045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreePP509056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreePP509067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreePP509078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreePP509089 ->Fill();
		  if(pt>=9 )
		    ftreePP509090 ->Fill();
		}
	      
	      }
	    }
	    
	    else{ 
	      
	      if(fWriteTree == kTRUE){
		if(centrality>=0 && centrality<20){
		  if(pt>=2 && pt < 3)
		    ftreeMM02023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMM02034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMM02045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMM02056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMM02067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMM02078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMM02089 ->Fill();
		  if(pt>=9 )
		    ftreeMM02090 ->Fill();
		}

		if(centrality>=20 && centrality<50){
		
		  if(pt>=2 && pt < 3)
		    ftreeMM205023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMM205034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMM205045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMM205056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMM205067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMM205078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMM205089 ->Fill();
		  if(pt>=9 )
		    ftreeMM205090 ->Fill();
		}
	      
	      
		if(centrality>=50 && centrality<90){
		  if(pt>=2 && pt < 3)
		    ftreeMM509023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMM509034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMM509045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMM509056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMM509067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMM509078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMM509089 ->Fill();
		  if(pt>=9 )
		    ftreeMM509090 ->Fill();
		}
	      
	      }
	    }
	    
	  }//fine like-sign	 
	  
	  // if(mass<0.2){
	  //   // cout<<mass<<endl;
	  //   // cout<<"Pt1: "<<vect1.Pt()<<" Pt2: "<<vect2.Pt()<<" M2: "<<vect2.M()<<endl;
	  // }
	}else {//Mixed-event pair histogramming
	  if (charge1 != charge2){ 
	  
	    fUSMassMix->Fill(mass,centrality);

	    if(centrality>=0 && centrality<10)
	      hLSMassvsPt0010 ->Fill(mass,pt);
	    if(centrality>=10 && centrality<40)	  
	      hLSMassvsPt1040 ->Fill(mass,pt);
	    if(centrality>=40 /*&& centrality<80*/)	    
	      hLSMassvsPt4080 ->Fill(mass,pt);

	    if(fWriteTree == kTRUE){
	      if(centrality>=0 && centrality<20){
		if(pt>=2 && pt < 3)
		  ftreeMixPM02023 ->Fill();
		if(pt>=3 && pt < 4)
		  ftreeMixPM02034 ->Fill();
		if(pt>=4 && pt < 5)
		  ftreeMixPM02045 ->Fill();
		if(pt>=5 && pt < 6)
		  ftreeMixPM02056 ->Fill();
		if(pt>=6 && pt < 7)
		  ftreeMixPM02067 ->Fill();
		if(pt>=7 && pt < 8)
		  ftreeMixPM02078 ->Fill();
		if(pt>=8 && pt < 9)
		  ftreeMixPM02089 ->Fill();
		if(pt>=9 )
		  ftreeMixPM02090 ->Fill();
	      }

	      if(centrality>=20 && centrality<50){
		
		if(pt>=2 && pt < 3)
		  ftreeMixPM205023 ->Fill();
		if(pt>=3 && pt < 4)
		  ftreeMixPM205034 ->Fill();
		if(pt>=4 && pt < 5)
		  ftreeMixPM205045 ->Fill();
		if(pt>=5 && pt < 6)
		  ftreeMixPM205056 ->Fill();
		if(pt>=6 && pt < 7)
		  ftreeMixPM205067 ->Fill();
		if(pt>=7 && pt < 8)
		  ftreeMixPM205078 ->Fill();
		if(pt>=8 && pt < 9)
		  ftreeMixPM205089 ->Fill();
		if(pt>=9 )
		  ftreeMixPM205090 ->Fill();
	      }
	      
	      
	      if(centrality>=50 && centrality<90){
		if(pt>=2 && pt < 3)
		  ftreeMixPM509023 ->Fill();
		if(pt>=3 && pt < 4)
		  ftreeMixPM509034 ->Fill();
		if(pt>=4 && pt < 5)
		  ftreeMixPM509045 ->Fill();
		if(pt>=5 && pt < 6)
		  ftreeMixPM509056 ->Fill();
		if(pt>=6 && pt < 7)
		  ftreeMixPM509067 ->Fill();
		if(pt>=7 && pt < 8)
		  ftreeMixPM509078 ->Fill();
		if(pt>=8 && pt < 9)
		  ftreeMixPM509089 ->Fill();
		if(pt>=9 )
		  ftreeMixPM509090 ->Fill();
	      }
	      
	    }
	  }
	  else{//like-sign
	    
	    if (charge1 > 0){ 
	      
	      if(fWriteTree == kTRUE){
		if(centrality>=0 && centrality<20){
		  if(pt>=2 && pt < 3)
		    ftreeMixPP02023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMixPP02034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMixPP02045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMixPP02056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMixPP02067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMixPP02078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMixPP02089 ->Fill();
		  if(pt>=9 )
		    ftreeMixPP02090 ->Fill();
		}

		if(centrality>=20 && centrality<50){
		
		  if(pt>=2 && pt < 3)
		    ftreeMixPP205023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMixPP205034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMixPP205045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMixPP205056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMixPP205067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMixPP205078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMixPP205089 ->Fill();
		  if(pt>=9 )
		    ftreeMixPP205090 ->Fill();
		}
	      
	      
		if(centrality>=50 && centrality<90){
		  if(pt>=2 && pt < 3)
		    ftreeMixPP509023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMixPP509034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMixPP509045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMixPP509056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMixPP509067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMixPP509078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMixPP509089 ->Fill();
		  if(pt>=9 )
		    ftreeMixPP509090 ->Fill();
		}
	      
	      }
	    }
	    
	    else{ //MM 
	      
	      if(fWriteTree == kTRUE){
		if(centrality>=0 && centrality<20){
		  if(pt>=2 && pt < 3)
		    ftreeMixMM02023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMixMM02034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMixMM02045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMixMM02056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMixMM02067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMixMM02078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMixMM02089 ->Fill();
		  if(pt>=9 )
		    ftreeMixMM02090 ->Fill();
		}

		if(centrality>=20 && centrality<50){
		
		  if(pt>=2 && pt < 3)
		    ftreeMixMM205023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMixMM205034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMixMM205045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMixMM205056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMixMM205067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMixMM205078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMixMM205089 ->Fill();
		  if(pt>=9 )
		    ftreeMixMM205090 ->Fill();
		}
	      
	      
		if(centrality>=50 && centrality<90){
		  if(pt>=2 && pt < 3)
		    ftreeMixMM509023 ->Fill();
		  if(pt>=3 && pt < 4)
		    ftreeMixMM509034 ->Fill();
		  if(pt>=4 && pt < 5)
		    ftreeMixMM509045 ->Fill();
		  if(pt>=5 && pt < 6)
		    ftreeMixMM509056 ->Fill();
		  if(pt>=6 && pt < 7)
		    ftreeMixMM509067 ->Fill();
		  if(pt>=7 && pt < 8)
		    ftreeMixMM509078 ->Fill();
		  if(pt>=8 && pt < 9)
		    ftreeMixMM509089 ->Fill();
		  if(pt>=9 )
		    ftreeMixMM509090 ->Fill();
		}
	      
	      }
	    }
	  }
	}//mix
	
      }//second ptc
    }
    if (evmultmixed!=0) multmixedcounted = kTRUE;
  } //first particle 
  
  
  
  if  (multmixedcounted) fHistMultiplicityOfMixedEvent->Fill(evmultmixed);//TBF
  
}


//__________________________________________________
void AnalysisDMTaskDiMuons::AnalizeMixing() {

}


//__________________________________________________
void AnalysisDMTaskDiMuons::Terminate() {
  fOutput->Write();
 
  fUSMass->Write();
  fUSMassC->Write();
  fUSMassLS->Write();
  fUSMassPP->Write();
  fUSMassMM->Write();
  fUSMassMix->Write();
  fUSMassMixLS->Write();
  fHistMultiplicityOfMixedEvent->Write();
  hZvertex->Write();
  hCentrality->Write();
  
  if(fWriteHisto == kTRUE){

    hRabs1->Write();
    hPt1->Write();
    hpDCA1->Write();
    hChi21->Write();
    hChi2Match1->Write();
    hRabs2->Write();
    hPt2->Write();
    hpDCA2->Write();
    hChi22->Write();
    hChi2Match2->Write();
    
    hpDCA1vsMass0080 ->Write();
    hpDCA2vsMass0080 ->Write();
    hpDCA1vsMass0010 ->Write();
    hpDCA2vsMass0010 ->Write();
    hpDCA1vsMass1040 ->Write();
    hpDCA2vsMass1040 ->Write();
    hpDCA1vsMass4080 ->Write();
    hpDCA2vsMass4080 ->Write();
    
    hpT1vsMass0080   ->Write();
    hpT2vsMass0080   ->Write();
    hpT1vsMass0010   ->Write();
    hpT2vsMass0010   ->Write();
    hpT1vsMass1040   ->Write();
    hpT2vsMass1040   ->Write();
    hpT1vsMass4080   ->Write();
    hpT2vsMass4080   ->Write();

    hpT1vspT20080->Write();
    hpT1vspT20010->Write();
    hpT1vspT21040->Write();
    hpT1vspT24080->Write();

    hpDCAvsPt10080   ->Write();
    hpDCAvsPt20080   ->Write();
    hpDCAvsPt10010   ->Write();
    hpDCAvsPt20010   ->Write();
    hpDCAvsPt11040   ->Write();
    hpDCAvsPt21040   ->Write();
    hpDCAvsPt14080   ->Write();
    hpDCAvsPt24080   ->Write();
    
    hMassvsPt0080  ->Write();
    hMassvsPt0010  ->Write();
    hMassvsPt1040  ->Write();
    hMassvsPt4080  ->Write();
    
    hLSMassvsPt0080 ->Write();
    hLSMassvsPt0010 ->Write();
    hLSMassvsPt1040 ->Write();
    hLSMassvsPt4080 ->Write();

    hPhivsMass0080  ->Write();
    hThetavsMass0080->Write();
    
    hPhi1vsMass0080 ->Write();
    hEta1vsMass0080 ->Write();
    hPhi2vsMass0080 ->Write();
    hEta2vsMass0080 ->Write();
    hMassvsRapidity0080->Write();
  }

  if(fWriteTree == kTRUE){

    //---------------------------

 ftreePM02023->Write();
 ftreePM02034->Write();
 ftreePM02045->Write();
 ftreePM02056->Write();
 ftreePM02067->Write();
 ftreePM02078->Write();
 ftreePM02089->Write();
 ftreePM02090->Write();
 
 ftreePM205023->Write();
 ftreePM205034->Write();
 ftreePM205045->Write();
 ftreePM205056->Write();
 ftreePM205067->Write();
 ftreePM205078->Write();
 ftreePM205089->Write();
 ftreePM205090->Write();
 
 ftreePM509023->Write();
 ftreePM509034->Write();
 ftreePM509045->Write();
 ftreePM509056->Write();
 ftreePM509067->Write();
 ftreePM509078->Write();
 ftreePM509089->Write();
 ftreePM509090->Write();
 
 //--//--------------

 ftreePP02023->Write();
 ftreePP02034->Write();
 ftreePP02045->Write();
 ftreePP02056->Write();
 ftreePP02067->Write();
 ftreePP02078->Write();
 ftreePP02089->Write();
 ftreePP02090->Write();
 
 ftreePP205023->Write();
 ftreePP205034->Write();
 ftreePP205045->Write();
 ftreePP205056->Write();
 ftreePP205067->Write();
 ftreePP205078->Write();
 ftreePP205089->Write();
 ftreePP205090->Write();
 
 ftreePP509023->Write();
 ftreePP509034->Write();
 ftreePP509045->Write();
 ftreePP509056->Write();
 ftreePP509067->Write();
 ftreePP509078->Write();
 ftreePP509089->Write();
 ftreePP509090->Write();
 
 //---------------------

 ftreeMM02023->Write();
 ftreeMM02034->Write();
 ftreeMM02045->Write();
 ftreeMM02056->Write();
 ftreeMM02067->Write();
 ftreeMM02078->Write();
 ftreeMM02089->Write();
 ftreeMM02090->Write();
 
 ftreeMM205023->Write();
 ftreeMM205034->Write();
 ftreeMM205045->Write();
 ftreeMM205056->Write();
 ftreeMM205067->Write();
 ftreeMM205078->Write();
 ftreeMM205089->Write();
 ftreeMM205090->Write();
 
 ftreeMM509023->Write();
 ftreeMM509034->Write();
 ftreeMM509045->Write();
 ftreeMM509056->Write();
 ftreeMM509067->Write();
 ftreeMM509078->Write();
 ftreeMM509089->Write();
 ftreeMM509090->Write();

 //------------------
 

 ftreeMixPM02023->Write();
 ftreeMixPM02034->Write();
 ftreeMixPM02045->Write();
 ftreeMixPM02056->Write();
 ftreeMixPM02067->Write();
 ftreeMixPM02078->Write();
 ftreeMixPM02089->Write();
 ftreeMixPM02090->Write();
 
 ftreeMixPM205023->Write();
 ftreeMixPM205034->Write();
 ftreeMixPM205045->Write();
 ftreeMixPM205056->Write();
 ftreeMixPM205067->Write();
 ftreeMixPM205078->Write();
 ftreeMixPM205089->Write();
 ftreeMixPM205090->Write();
 
 ftreeMixPM509023->Write();
 ftreeMixPM509034->Write();
 ftreeMixPM509045->Write();
 ftreeMixPM509056->Write();
 ftreeMixPM509067->Write();
 ftreeMixPM509078->Write();
 ftreeMixPM509089->Write();
 ftreeMixPM509090->Write();
 
 //----------------

 ftreeMixPP02023->Write();
 ftreeMixPP02034->Write();
 ftreeMixPP02045->Write();
 ftreeMixPP02056->Write();
 ftreeMixPP02067->Write();
 ftreeMixPP02078->Write();
 ftreeMixPP02089->Write();
 ftreeMixPP02090->Write();
 
 ftreeMixPP205023->Write();
 ftreeMixPP205034->Write();
 ftreeMixPP205045->Write();
 ftreeMixPP205056->Write();
 ftreeMixPP205067->Write();
 ftreeMixPP205078->Write();
 ftreeMixPP205089->Write();
 ftreeMixPP205090->Write();
 
 ftreeMixPP509023->Write();
 ftreeMixPP509034->Write();
 ftreeMixPP509045->Write();
 ftreeMixPP509056->Write();
 ftreeMixPP509067->Write();
 ftreeMixPP509078->Write();
 ftreeMixPP509089->Write();
 ftreeMixPP509090->Write();

 //---------------------

 ftreeMixMM02023->Write();
 ftreeMixMM02034->Write();
 ftreeMixMM02045->Write();
 ftreeMixMM02056->Write();
 ftreeMixMM02067->Write();
 ftreeMixMM02078->Write();
 ftreeMixMM02089->Write();
 ftreeMixMM02090->Write();
 
 ftreeMixMM205023->Write();
 ftreeMixMM205034->Write();
 ftreeMixMM205045->Write();
 ftreeMixMM205056->Write();
 ftreeMixMM205067->Write();
 ftreeMixMM205078->Write();
 ftreeMixMM205089->Write();
 ftreeMixMM205090->Write();
 
 ftreeMixMM509023->Write();
 ftreeMixMM509034->Write();
 ftreeMixMM509045->Write();
 ftreeMixMM509056->Write();
 ftreeMixMM509067->Write();
 ftreeMixMM509078->Write();
 ftreeMixMM509089->Write();
 ftreeMixMM509090->Write();
 
}
  
  // for (auto i = 0; i < fNCentralityBins; ++i) {
  //   fUSMassCentrality[i]->Write();
  //   fUSMassCentralityLS[i]->Write();
  // }
  //  fMixingPool->Write();
  fOutput->Close();
}


//__________________________________________________
Bool_t AnalysisDMTaskDiMuons::HasTriggerFired(AnalysisOption opt) {
  Short_t triggers = GetTriggerMask(fEvt->GetTriggerString());
  
  // printf("%d\t%d\t%d\t%d\n", triggers, fTriggerIncludedData, fTriggerExcludedData, (triggers&fTriggerIncludedData) &&! (triggers&fTriggerExcludedData));
  if (opt == kDirect) {
    return (triggers&fTriggerIncludedData) && !(triggers&fTriggerExcludedData);
  } else if (opt == kMixing) {
    return (triggers&fTriggerIncludedMixing) && !(triggers&fTriggerExcludedMixing);
  }
  return kFALSE;
}

//__________________________________________________

Bool_t AnalysisDMTaskDiMuons::SelectDimu(TLorentzVector *dimu) { 
  Bool_t isAccepted = kFALSE; 
  // insert here possible cuts on dimu ....
  if (dimu->Rapidity() < -2.5 && dimu->Rapidity() > -4) isAccepted = kTRUE; 
  return isAccepted; 
}

//__________________________________________________
Bool_t AnalysisDMTaskDiMuons::HasPassedEventCuts(AliLMREvent *event) {
  if (!(event->GetPhysicsSelectionMask()&fPS)) return kFALSE;
  if (!(event->GetZVertex() > fZvxEventMin  && event->GetZVertex() < fZvxEventMax)) return kFALSE;
  if (event->GetNVtxContributors() < 1) return kFALSE;
  if (event->IsPileupFromSPD()) return kFALSE;
  if (!(event->GetMultiplicity(fActivitySelection) < fActivityMax)) return kFALSE;
  return kTRUE;
}


//__________________________________________________
Bool_t AnalysisDMTaskDiMuons::HasPassedSingleMuCuts(AliLMRMuon *mu) {
  if (!fMu1) return kFALSE;
  TLorentzVector vMu = mu->P4();
  if (!(mu->GetTriggerMatch() >= fTrigLevel)) return kFALSE;
  if (!(vMu.Eta() > fEtaSingleMuMin && vMu.Eta() < fEtaSingleMuMax)) return kFALSE;
  if (!(vMu.Pt()  > fPtSingleMuMin  && vMu.Pt()  < fPtSingleMuMax))  return kFALSE;

 
  // hRabs1->Fill(fMu1->GetRabs());
  // hPt1->Fill(pt1);
  // hpDCA1->Fill(fMu1->GetpDCA());
  // hChi21->Fill(fMu1->GetChi2());
  // hChi2Match1->Fill(fMu1->GetChi2Match());
      
  return kTRUE;
}


// //__________________________________________________ DEFAULT trigger
// void AnalysisDMTaskDiMuons::SetTrigger(TriggerOptions optData, TriggerOptions optMixing) {
//   fTriggerIncludedData   = optData;
//   fTriggerIncludedMixing = optMixing;

//   switch (fTriggerIncludedData) {
//     case kMSH:
//     case kMSL:
//     case kMUL:
//     case kMLL:
//     case kMB:
//     case kSingleMuon:
//     case kDimuon:
//       fIsTriggerSet = kTRUE;
//       break;

//     default:
//     AliError("No valid trigger option for Data!!!");
//     fIsTriggerSet = kFALSE;
//     return;
//   }

//   switch (fTriggerIncludedMixing) {
//     case kMSH:
//     case kMSL:
//     case kMUL:
//     case kMLL:
//     case kSingleMuon:
//     case kDimuon:
//       fIsTriggerSet = kTRUE;
//       break;

//     default:
//       AliError("No valid trigger option for Mixing!!!");
//       fIsTriggerSet = kFALSE;
//       return;
//   }
// }

//__________________________________________________
void AnalysisDMTaskDiMuons::SetTrigger(TriggerOptions optData, TriggerOptions optMixing) {
  fTriggerIncludedData   = optData;
  fTriggerIncludedMixing = optMixing;

  switch (fTriggerIncludedData) {
  case kMUL:
    fIsTriggerSet = kTRUE;
    break;
  case kDimuon:
    fIsTriggerSet = kTRUE;
    break;
    
    default:
      AliError("No valid trigger option for Data!!!");
      fIsTriggerSet = kFALSE;
      return;
  }
  
  switch (fTriggerIncludedMixing) {
  case kMUL:
    fIsTriggerSet = kTRUE;
      break;
      
  default:
    AliError("No valid trigger option for Mixing!!!");
    fIsTriggerSet = kFALSE;
    return;
  }
}


//__________________________________________________
void AnalysisDMTaskDiMuons::SetNOTTrigger(TriggerOptions optData, TriggerOptions optMixing) {
  fTriggerExcludedData   = optData;
  fTriggerExcludedMixing = optMixing;
}


//__________________________________________________
Short_t AnalysisDMTaskDiMuons::GetTriggerMask(TString word) {
  Short_t mask = 0;
  if (word.Contains("CMSH7")) mask|=kMSH;
  if (word.Contains("CMSL7")) mask|=kMSL;
  if (word.Contains("CMUL7")) mask|=kMUL;
  if (word.Contains("CMLL7")) mask|=kMLL;
  return mask;
}


//__________________________________________________
void AnalysisDMTaskDiMuons::BuildPoolForMixing() {
  for (auto ievt = 0; ievt < fNEvents; ++ievt) {
    fEvtTree->GetEntry(ievt);

    // Trigger selection
    if (!(HasTriggerFired(kMixing))) continue;

    // Event cuts
    if (!HasPassedEventCuts(fEvt)) continue;

    Int_t nMuons = fEvt->GetNMuons();
    if (nMuons < 1) continue;

    Double_t centrality = TMath::Min(fEvt->GetMultiplicity(fActivitySelection), 99.9f);

    for (Int_t iMu = 0; iMu < nMuons; ++iMu) {
      fMu1 = (AliLMRMuon*) fEvt->GetMuon(iMu);

    }

  }
}


//__________________________________________________
void AnalysisDMTaskDiMuons::SetCentralityBins(Int_t nBins, Double_t bins[]) {
  /*
  fNCentralityBins = nBins-1;
  fCentralityBins = new TAxis(fNCentralityBins, bins);
  cout << "Centrality bins" << endl;
  cout << "---------------" << endl;
  for (auto i = 1; i <= fCentralityBins->GetNbins(); ++i) {
    cout << "    " << setw(2) << fCentralityBins->GetBinLowEdge(i) << " - " << fCentralityBins->GetBinUpEdge(i) << endl;
  }
  cout << endl;

  // US pairs plots
  fUSMassCentrality = new TH1D*[fNCentralityBins];
  fUSMassCentralityLS = new TH1D*[fNCentralityBins];
  string name, crange, name2;
  for (auto i = 0; i < fNCentralityBins; ++i) {
    crange = to_string((int) fCentralityBins->GetBinLowEdge(i+1)) + "_" + to_string((int) fCentralityBins->GetBinUpEdge(i+1));
    name   = "fUSMass_"   + crange;
    name2  = "fUSMassLS_" + crange;
    fUSMassCentrality[i]= new TH1D(name.c_str(), name.c_str(), 20000, 0, 20);
    fUSMassCentralityLS[i]= new TH1D(name2.c_str(), name2.c_str(), 20000, 0, 20);
  }

  // 1: number of events - 2: centrality
  Int_t mixingPoolBins[2]    = {    fNEvents, fNCentralityBins};
  Double_t mixingPoolMins[2] = {        -0.5,   0.};
  Double_t mixingPoolMaxs[2] = {fNEvents-0.5, 100.};
  fMixingPool = new THnSparseI("fMixingPool", "fMixingPool", 2, mixingPoolBins, mixingPoolMins, mixingPoolMaxs);
  fMixingPool->GetAxis(1)->Set(fNCentralityBins, bins);
  */
  
}
