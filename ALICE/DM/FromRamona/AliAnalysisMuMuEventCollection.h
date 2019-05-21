#ifndef AliAnalysisMuMuEventCollection_cxx
#define AliAnalysisMuMuEventCollection_cxx
#include <iostream>

#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "AliLMRMuon.h"

//From FEMTO
 
using namespace std;

class AliReconstructedFirst {
  
 public:

  AliReconstructedFirst();
  
  ~AliReconstructedFirst();
  
  /* AliReconstructedFirst(const AliReconstructedFirst&); */
  /* AliReconstructedFirst & operator=(const AliReconstructedFirst&); */
  
  TLorentzVector fVectMu;
  Short_t        charge;

  ClassDef(AliReconstructedFirst, 1);   

};

class AliReconstructedSecond {
  
 public:

  AliReconstructedSecond();
  
  ~AliReconstructedSecond();
  
  /* AliReconstructedSecond(const AliReconstructedSecond&); */
  /* AliReconstructedSecond & operator=(const AliReconstructedSecond&); */
  
  TLorentzVector sVectMu;
  Short_t        charge;

  ClassDef(AliReconstructedSecond, 1);   

};


class AliAnalysisMuMuEvent {
  
 public:
  
  int fNumberCandidateFirst;
  int fNumberCandidateSecond;
  
  AliReconstructedFirst *fReconstructedFirst;
  AliReconstructedSecond *fReconstructedSecond;
  
  ClassDef(AliAnalysisMuMuEvent, 1);
  
};



class AliAnalysisMuMuEventCollection  {

   public:

    AliAnalysisMuMuEventCollection();

    AliAnalysisMuMuEventCollection(short eventBuffSize, int maxFirstMult, int maxSecondMult);
    
    ~AliAnalysisMuMuEventCollection();
    
    void FifoShift();
    
    AliAnalysisMuMuEvent *fEvt;
    
 private:
    
    short fifo; //Size of the Event Storage buffer
    
    void SetBuffSize(short eventBuffSize){fifo = eventBuffSize;}
    
    ClassDef(AliAnalysisMuMuEventCollection, 1);
    
};




#endif
