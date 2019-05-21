////////////////////////////////////////////////////////////////////////////////
//                       To be implemented
////////////////////////////////////////////////////////////////////////////////



#include "AliAnalysisMuMuEventCollection.h"
#include "TLorentzVector.h"


AliAnalysisMuMuEventCollection::~AliAnalysisMuMuEventCollection(){

  for(int i = 0; i < fifo; i++){
     if((fEvt + i)->fReconstructedFirst != NULL){
       delete [] (fEvt + i)->fReconstructedFirst;
     }
     if((fEvt + i)->fReconstructedSecond != NULL){
       delete [] (fEvt + i)->fReconstructedSecond;
     }
   }
   delete [] fEvt;
 }
 //_____________________________________________________________________________


AliAnalysisMuMuEventCollection::AliAnalysisMuMuEventCollection() : fEvt(0x0), fifo(0) {


}

//_____________________________________________________________________________

AliAnalysisMuMuEventCollection::AliAnalysisMuMuEventCollection(short eventBuffSize, int maxFirstMult, int maxSecondMult) : fEvt(0x0), fifo(0) {

  SetBuffSize(eventBuffSize);

  fEvt = new AliAnalysisMuMuEvent[fifo];  //allocate pointer array of AliAnalysisMuMuEvents

  for(int ii = 0; ii < fifo; ii++){ //Initialize particle table pointers to NULL

    (fEvt + ii)->fReconstructedFirst = NULL;

    (fEvt + ii)->fNumberCandidateFirst = 0;

    (fEvt + ii)->fReconstructedFirst = new AliReconstructedFirst[maxFirstMult];


    (fEvt + ii)->fReconstructedSecond = NULL;

    (fEvt + ii)->fNumberCandidateSecond = 0;

    (fEvt + ii)->fReconstructedSecond = new AliReconstructedSecond[maxSecondMult];


  }

}


//_____________________________________________________________________________
void AliAnalysisMuMuEventCollection::FifoShift() { //Shift elements in FIFO by one and clear last element in FIFO 

  for(unsigned short i=fifo-1 ; i > 0; i--) {

    for(int j=0; j<(fEvt + i-1)->fNumberCandidateFirst; j++){

      (fEvt + i)->fReconstructedFirst[j] = (fEvt + i-1)->fReconstructedFirst[j];

    }

    (fEvt + i)->fNumberCandidateFirst = (fEvt + i-1)->fNumberCandidateFirst;

    for(int j=0; j<(fEvt + i-1)->fNumberCandidateSecond; j++){

      (fEvt + i)->fReconstructedSecond[j] = (fEvt + i-1)->fReconstructedSecond[j];

    }

    (fEvt + i)->fNumberCandidateSecond = (fEvt + i-1)->fNumberCandidateSecond;

    // for(int j=0; j<3; j++){

    //   (fEvt + i)->fPrimaryVertex[j] = (fEvt + i-1)->fPrimaryVertex[j];

    // }

  }

   (fEvt)->fNumberCandidateFirst=0;
   (fEvt)->fNumberCandidateSecond=0;

   // for(int j=0; j<3; j++) {

   //   (fEvt)->fPrimaryVertex[j] = 0.;

   // }

}



//_____________________________________________________________________________

AliReconstructedFirst::AliReconstructedFirst() :
  
  fVectMu(),
  charge(0)
  
{
  
}
//_____________________________________________________________________________

AliReconstructedSecond::AliReconstructedSecond() :

  sVectMu(),
  charge(0)
  
  
{

}
//_____________________________________________________________________________

AliReconstructedFirst::~AliReconstructedFirst()

{
  
}

//_____________________________________________________________________________

AliReconstructedSecond::~AliReconstructedSecond()

{
  
}
