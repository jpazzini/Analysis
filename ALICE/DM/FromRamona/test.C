void test(Int_t nrun = 296065 /*245068*/) {
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  // gROOT->LoadMacro("AliAnalysisMuMuEventCollection.cxx++g");
  // gROOT->LoadMacro("AnalysisDMTaskDiMuons.cxx++g");

  gROOT->LoadMacro("AliAnalysisMuMuEventCollection.cxx+g");
  gROOT->LoadMacro("AnalysisDMTaskDiMuons.cxx+g");

  TString methodActivity = "V0M";
  Double_t activityMax   = 250.0; //data

  enum PhysicsSelection{kNoPS=1<<0, kINT7=1<<1, kINT8=1<<2, kINT7MUON=1<<3, kINT8MUON=1<<4};

  AnalysisDMTaskDiMuons* task = new AnalysisDMTaskDiMuons();
  //  if (!task->Init(246994)) return;
  task->Init(nrun);
  //task->SetTrigger(TriggerOptions::kDimuon, TriggerOptions::kMSL); //uncomment
  task->SetTrigger(TriggerOptions::kDimuon, TriggerOptions::kMUL); 
  task->SetNOTTrigger(0, TriggerOptions::kDimuon);
  task->SetPhysicsSelection(kINT7MUON);
  //  task->SetZvxEvent(-2000., 2000.);
  //  task->SetZvxEvent(-20., 20.);
  task->SetZvxEvent(-10., 10.);
  task->SetMatchTriggerLevel(2);
  task->SetEtaSingleMuCut(-4., -2.5);
  //  task->SetPtSingleMuCut(0.25, 99999999.);
  //  task->SetPtSingleMuCut(0.25, 100.);
  task->SetPtSingleMuCut(.25, 100.);
  // task->SetPtSingleMuCut(1, 100.);
  task->SetActivitySelection(methodActivity);
  task->SetActivityMax(activityMax);
 
  // task->SetWriteTree(kTRUE);
  // task->SetWriteHisto(kFALSE);

  task->SetWriteTree(kTRUE);
  task->SetWriteHisto(kTRUE);



  Double_t centralityBins[] = {0., 10., 40., 100.};
  task->SetCentralityBins(sizeof(centralityBins)/sizeof(double), centralityBins);

  task->Analize();
}
