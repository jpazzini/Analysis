/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was generated by  Zijun Xu                           * 
 *****************************************************************************/ 

#include "Riostream.h" 

#include "HWWLVJRooPdfs.h"
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooExponential.h" 
#include <math.h> 
#include "TMath.h" 

#include <algorithm>
#include <vector>
#include <string>

#include "RooPlot.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TIterator.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooCurve.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
using namespace std;

void HWWLVJRooPdfs(){}

//// Erf*Exp function implementation 
Double_t ErfExp(Double_t x, Double_t c, Double_t offset, Double_t width){
    if(width<1e-2)width=1e-2;
    if (c==0)c=-1e-7;
	return TMath::Exp(c*x)*(1.+TMath::Erf((x-offset)/width))/2. ;
}

Double_t ErfExp(Double_t x, Double_t x_min, Double_t x_max, Double_t c, Double_t offset, Double_t width){
    if(width<1e-2)width=1e-2;
    if (c==0)c=1e-7;
    double minTerm = (TMath::Exp(c*c*width*width/4+c*offset) * 
					TMath::Erf((2*x_min-c*width*width-
							2*offset)/2/width) - 
					TMath::Exp(c*x_min) * 
					TMath::Erf((x_min-offset)/width) - 
					TMath::Exp(c*x_min))/-2/c;
	double maxTerm = (TMath::Exp(c*c*width*width/4+c*offset) * 
					TMath::Erf((2*x_max-c*width*width-
							2*offset)/2/width) - 
					TMath::Exp(c*x_max) * 
					TMath::Erf((x_max-offset)/width) - 
					TMath::Exp(c*x_max))/-2/c;
	Double_t integral=(maxTerm-minTerm) ;
	return TMath::Exp(c*x)*(1.+TMath::Erf((x-offset)/width))/2./integral ;
}

//// Single Exp function 
Double_t Exp(Double_t x, Double_t c){
	return TMath::Exp(c*x);
}

Double_t Exp(Double_t x, Double_t x_min, Double_t x_max, Double_t c){
	Double_t integral ;
    if(c==0.){
        integral=x_max-x_min;
    }else{
        integral= ( TMath::Exp(c*x_max)-TMath::Exp(c*x_min) ) / c;
    }
	return TMath::Exp(c*x)/integral ;
}


//// Erf*Exp pdf 
ClassImp(RooErfExpPdf) 

RooErfExpPdf::RooErfExpPdf(const char *name, const char *title, 
					RooAbsReal& _x,
					RooAbsReal& _c,
					RooAbsReal& _offset,
					RooAbsReal& _width) :
    			                RooAbsPdf(name,title), 
 			                x("x","x",this,_x),
			                c("c","c",this,_c),
			                offset("offset","offset",this,_offset),
			                width("width","width",this,_width){ } 


RooErfExpPdf::RooErfExpPdf(const RooErfExpPdf& other, const char* name) :  
	RooAbsPdf(other,name), 
	x("x",this,other.x),
	c("c",this,other.c),
	offset("offset",this,other.offset),
	width("width",this,other.width){}



Double_t RooErfExpPdf::evaluate() const { 

    Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
    return ErfExp(x,c,offset,width_tmp) ; 
} 

Int_t RooErfExpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  { 

	if (matchArgs(allVars,analVars,x)) return 1 ; 
	return 0 ; 
} 

Double_t RooErfExpPdf::analyticalIntegral(Int_t code, const char* rangeName) const  { 

    Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
	if (code==1) { 
        Double_t minTerm=0;
        Double_t maxTerm=0;
        if(c==0){ 
            Double_t delta=-1e-7;
            minTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
					TMath::Erf((2*x.min(rangeName)-delta*width_tmp*width_tmp-
							2*offset)/2/width_tmp) - 
					TMath::Exp(delta*x.min(rangeName)) * 
					TMath::Erf((x.min(rangeName)-offset)/width_tmp) - 
					TMath::Exp(delta*x.min(rangeName)))/-2/delta;
		    maxTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
					TMath::Erf((2*x.max(rangeName)-delta*width_tmp*width_tmp-
							2*offset)/2/width_tmp) - 
					TMath::Exp(delta*x.max(rangeName)) * 
					TMath::Erf((x.max(rangeName)-offset)/width_tmp) - 
					TMath::Exp(delta*x.max(rangeName)))/-2/delta;
        
        }else{
            minTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
					TMath::Erf((2*x.min(rangeName)-c*width_tmp*width_tmp-
							2*offset)/2/width_tmp) - 
					TMath::Exp(c*x.min(rangeName)) * 
					TMath::Erf((x.min(rangeName)-offset)/width_tmp) - 
					TMath::Exp(c*x.min(rangeName)))/-2/c;
		    maxTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
					TMath::Erf((2*x.max(rangeName)-c*width_tmp*width_tmp-
							2*offset)/2/width_tmp) - 
					TMath::Exp(c*x.max(rangeName)) * 
					TMath::Erf((x.max(rangeName)-offset)/width_tmp) - 
					TMath::Exp(c*x.max(rangeName)))/-2/c;
        }
		return (maxTerm-minTerm) ;
	} 
	return 0 ; 
} 

/// RooAlpha pdf as ratio of two Erf*Exp

ClassImp(RooAlpha)

RooAlpha::RooAlpha(){}

RooAlpha::RooAlpha(const char *name, const char *title,
		   RooAbsReal& _x,
		   RooAbsReal& _c,
		   RooAbsReal& _offset,
		   RooAbsReal& _width,
		   RooAbsReal& _ca,
		   RooAbsReal& _offseta,
		   RooAbsReal& _widtha,
                   Double_t _xmin,
                   Double_t _xmax) :
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  c("c","c",this,_c),
  offset("offset","offset",this,_offset),
  width("width","width",this,_width),
  ca("ca","ca",this,_ca),
  offseta("offseta","offseta",this,_offseta),
  widtha("widtha","widtha",this,_widtha){
        xmin=_xmin;
        xmax=_xmax;
}

RooAlpha::RooAlpha(const RooAlpha& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  c("c",this,other.c),
  offset("offset",this,other.offset),
  width("width",this,other.width),
  ca("ca",this,other.ca),
  offseta("offseta",this,other.offseta),
  widtha("widtha",this,other.widtha){
        xmin=other.xmin;
        xmax=other.xmax;
}

double RooAlpha::evaluate() const{
    Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
    Double_t widtha_tmp=widtha; if(widtha<1e-2){ widtha_tmp=1e-2;}
    return ErfExp(x,xmin,xmax,c,offset,width_tmp)/ErfExp(x,xmin,xmax,ca,offseta,widtha_tmp);
}


/// Alpha function given by the ratio of two exponential functions
ClassImp(RooAlpha)

RooAlphaExp::RooAlphaExp(){}

RooAlphaExp::RooAlphaExp(const char *name, const char *title,
		   RooAbsReal& _x,
		   RooAbsReal& _c,
		   RooAbsReal& _ca,
                   Double_t _xmin,
                   Double_t _xmax) :
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  c("c","c",this,_c),
  ca("ca","ca",this,_ca){
        xmin=_xmin;
        xmax=_xmax;
}

RooAlphaExp::RooAlphaExp(const RooAlphaExp& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  c("c",this,other.c),
  ca("ca",this,other.ca){
        xmin=other.xmin;
        xmax=other.xmax;
}

double RooAlphaExp::evaluate() const{
  return Exp(x,xmin,xmax,c)/Exp(x,xmin,xmax,ca);
}


///// Relativistic BW
ClassImp(RooBWRunPdf) 

 RooBWRunPdf::RooBWRunPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mean,
                        RooAbsReal& _width) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width){ } 


 RooBWRunPdf::RooBWRunPdf(const RooBWRunPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width){ } 

 Double_t RooBWRunPdf::evaluate() const { 
   return (x*x*width/mean) / ( (x*x-mean*mean)*(x*x-mean*mean) + (x*x*width/mean)*(x*x*width/mean) );
 } 

///// Erf*Pow2 pdf 
ClassImp(RooErfPow2Pdf) 

Double_t  ErfPow2(Double_t x,Double_t c0,Double_t c1, Double_t offset, Double_t width){

   if(width<1e-2)width=1e-2;
   Double_t sqrt_s=2000.;
   return TMath::Power(x/sqrt_s ,-1*(c0+c1*TMath::Log(x/sqrt_s)) )*(1+ TMath::Erf((x-offset)/width)) /2. ; 
 }

RooErfPow2Pdf::RooErfPow2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _c1,
                        RooAbsReal& _offset,
                        RooAbsReal& _width) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   c1("c1","c1",this,_c1),
   offset("offset","offset",this,_offset),
   width("width","width",this,_width){} 


RooErfPow2Pdf::RooErfPow2Pdf(const RooErfPow2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   c1("c1",this,other.c1),
   offset("offset",this,other.offset),
   width("width",this,other.width){ } 

Double_t RooErfPow2Pdf::evaluate() const { 
   Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
   return ErfPow2(x,c0,c1,offset,width_tmp);
 } 


///// Alpha  function given by the ratio of two Erf*Pow2 pdf

ClassImp(RooAlpha4ErfPow2Pdf) 

RooAlpha4ErfPow2Pdf::RooAlpha4ErfPow2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _c1,
                        RooAbsReal& _offset,
                        RooAbsReal& _width,
                        RooAbsReal& _c0a,
                        RooAbsReal& _c1a,
                        RooAbsReal& _offseta,
                        RooAbsReal& _widtha) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   c1("c1","c1",this,_c1),
   offset("offset","offset",this,_offset),
   width("width","width",this,_width),
   c0a("c0a","c0a",this,_c0a),
   c1a("c1a","c1a",this,_c1a),
   offseta("offseta","offseta",this,_offseta),
   widtha("widtha","widtha",this,_widtha){} 


RooAlpha4ErfPow2Pdf::RooAlpha4ErfPow2Pdf(const RooAlpha4ErfPow2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   c1("c1",this,other.c1),
   offset("offset",this,other.offset),
   width("width",this,other.width),
   c0a("c0a",this,other.c0a),
   c1a("c1a",this,other.c1a),
   offseta("offseta",this,other.offseta),
   widtha("widtha",this,other.widtha){} 



Double_t RooAlpha4ErfPow2Pdf::evaluate() const { 

    Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
    Double_t widtha_tmp=widtha; if(widtha<1e-2){ widtha_tmp=1e-2;}
    return ErfPow2(x,c0,c1,offset,width_tmp)/ErfPow2(x,c0a,c1a,offseta,widtha_tmp);
 } 


//////// Erf Pow Exp pdf 
ClassImp(RooErfPowExpPdf) 

Double_t  ErfPowExp(Double_t x,Double_t c0,Double_t c1, Double_t offset, Double_t width){
   if(width<1e-2)width=1e-2;
   Double_t sqrt_s=2000.;
   return TMath::Power(x/sqrt_s ,-1*(c1*TMath::Log(x/sqrt_s)) )*TMath::Exp(-1*x/sqrt_s*c0)*(1+ TMath::Erf((x-offset)/width)) /2. ; 
}

RooErfPowExpPdf::RooErfPowExpPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _c1,
                        RooAbsReal& _offset,
                        RooAbsReal& _width) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   c1("c1","c1",this,_c1),
   offset("offset","offset",this,_offset),
   width("width","width",this,_width){} 


RooErfPowExpPdf::RooErfPowExpPdf(const RooErfPowExpPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   c1("c1",this,other.c1),
   offset("offset",this,other.offset),
   width("width",this,other.width){} 


Double_t RooErfPowExpPdf::evaluate() const { 
   Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
   return ErfPowExp(x,c0,c1,offset,width_tmp);
} 


////// Alpha by the ratio of two ErfPowExp Pdf

ClassImp(RooAlpha4ErfPowExpPdf) 

RooAlpha4ErfPowExpPdf::RooAlpha4ErfPowExpPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _c1,
                        RooAbsReal& _offset,
                        RooAbsReal& _width,
                        RooAbsReal& _c0a,
                        RooAbsReal& _c1a,
                        RooAbsReal& _offseta,
                        RooAbsReal& _widtha) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   c1("c1","c1",this,_c1),
   offset("offset","offset",this,_offset),
   width("width","width",this,_width),
   c0a("c0a","c0a",this,_c0a),
   c1a("c1a","c1a",this,_c1a),
   offseta("offseta","offseta",this,_offseta),
   widtha("widtha","widtha",this,_widtha){} 


RooAlpha4ErfPowExpPdf::RooAlpha4ErfPowExpPdf(const RooAlpha4ErfPowExpPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   c1("c1",this,other.c1),
   offset("offset",this,other.offset),
   width("width",this,other.width),
   c0a("c0a",this,other.c0a),
   c1a("c1a",this,other.c1a),
   offseta("offseta",this,other.offseta),
   widtha("widtha",this,other.widtha){} 

Double_t RooAlpha4ErfPowExpPdf::evaluate() const { 
   Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
   Double_t widtha_tmp=widtha; if(widtha<1e-2){ widtha_tmp=1e-2;}
   return ErfPowExp(x,c0,c1,offset,width_tmp)/ErfPowExp(x,c0a,c1a,offseta,widtha_tmp);
} 

////// Gaus Exp Pdf 
ClassImp(RooGausExpPdf) 

Double_t  GausExp(Double_t x,Double_t c,Double_t mean, Double_t sigma){
        if(sigma<1e-2)sigma=1e-2;
	return TMath::Exp(c*x)+TMath::Exp(-(x-mean)*(x-mean)/(2*sigma*sigma)) ; 
}

RooGausExpPdf::RooGausExpPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _mean,
                        RooAbsReal& _sigma) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   mean("mean","mean",this,_mean),
   sigma("sigma","sigma",this,_sigma){ } 


RooGausExpPdf::RooGausExpPdf(const RooGausExpPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   mean("mean",this,other.mean),
   sigma("sigma",this,other.sigma){} 


Double_t RooGausExpPdf::evaluate() const { 
   Double_t width_tmp=sigma; if(sigma<1e-2){ width_tmp=1e-2;}
   return GausExp(x,c,mean,width_tmp);
} 


/////////////// Alpha for Gaus Exp Function

ClassImp(RooGausExpPdf) 

RooAlpha4GausExpPdf::RooAlpha4GausExpPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _mean,
                        RooAbsReal& _sigma,
                        RooAbsReal& _ca,
                        RooAbsReal& _meana,
		        RooAbsReal& _sigmaa):
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   mean("mean","mean",this,_mean),
   sigma("sigma","sigma",this,_sigma),
   ca("ca","ca",this,_ca),
   meana("meana","meana",this,_meana),
   sigmaa("sigmaa","sigmaa",this,_sigmaa){} 


RooAlpha4GausExpPdf::RooAlpha4GausExpPdf(const RooAlpha4GausExpPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   mean("mean",this,other.mean),
   sigma("sigma",this,other.sigma),
   ca("ca",this,other.ca),
   meana("meana",this,other.meana),
   sigmaa("sigmaa",this,other.sigmaa){} 

Double_t RooAlpha4GausExpPdf::evaluate() const { 
    Double_t width_tmp=sigma; if(sigma<1e-2){ width_tmp=1e-2;}
    Double_t widtha_tmp=sigmaa; if(sigmaa<1e-2){ widtha_tmp=1e-2;}
    return GausExp(x,c,mean,width_tmp)/GausExp(x,ca,meana,widtha_tmp);} 


////// Erf*Pow Pdf 
ClassImp(RooErfPowPdf) 

Double_t  ErfPow(Double_t x,Double_t c, Double_t offset, Double_t width){
   if(width<1e-2)width=1e-2;
   Double_t sqrt_s=2000.;
   return TMath::Power(x/sqrt_s ,c)*(1+ TMath::Erf((x-offset)/width)) /2. ; 
}

RooErfPowPdf::RooErfPowPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _offset,
                        RooAbsReal& _width) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   offset("offset","offset",this,_offset),
   width("width","width",this,_width){ } 


RooErfPowPdf::RooErfPowPdf(const RooErfPowPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   offset("offset",this,other.offset),
   width("width",this,other.width){} 


Double_t RooErfPowPdf::evaluate() const { 
   Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
   return ErfPow(x,c,offset,width_tmp);} 


//////// Alpha given by the ratio of two Erf*Pow

ClassImp(RooAlpha4ErfPowPdf) 

RooAlpha4ErfPowPdf::RooAlpha4ErfPowPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _offset,
                        RooAbsReal& _width,
                        RooAbsReal& _ca,
                        RooAbsReal& _offseta,
                        RooAbsReal& _widtha) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   offset("offset","offset",this,_offset),
   width("width","width",this,_width),
   ca("ca","ca",this,_ca),
   offseta("offseta","offseta",this,_offseta),
   widtha("widtha","widtha",this,_widtha){} 


RooAlpha4ErfPowPdf::RooAlpha4ErfPowPdf(const RooAlpha4ErfPowPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   offset("offset",this,other.offset),
   width("width",this,other.width),
   ca("ca",this,other.ca),
   offseta("offseta",this,other.offseta),
   widtha("widtha",this,other.widtha){} 



Double_t RooAlpha4ErfPowPdf::evaluate() const { 
   Double_t width_tmp=width; if(width<1e-2){ width_tmp=1e-2;}
   Double_t widtha_tmp=widtha; if(widtha<1e-2){ widtha_tmp=1e-2;}
   return ErfPow(x,c,offset,width_tmp)/ErfPow(x,ca,offseta,widtha_tmp);
} 

//////////////////////////////////////////RooPow2Pdf.cxx
ClassImp(RooPow2Pdf) 

RooPow2Pdf::RooPow2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _p0,
                        RooAbsReal& _p1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   p0("p0","p0",this,_p0),
   p1("p1","p1",this,_p1){} 


RooPow2Pdf::RooPow2Pdf(const RooPow2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   p0("p0",this,other.p0),
   p1("p1",this,other.p1){} 

Double_t RooPow2Pdf::evaluate() const { 
   Double_t sqrt_s=2000.;
   return TMath::Power( x/sqrt_s,-1*( p0+p1*TMath::Log(x/sqrt_s) ) )  ; 
} 

////////////////////////////RooPowPdf.cxx
ClassImp(RooPowPdf) 

RooPowPdf::RooPowPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _p0) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   p0("p0","p0",this,_p0){} 


RooPowPdf::RooPowPdf(const RooPowPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   p0("p0",this,other.p0){} 



Double_t RooPowPdf::evaluate() const { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   Double_t sqrt_s=2000.;
   return TMath::Power( x/sqrt_s, p0 )  ; 
} 

/////////////////////////////  RooQCDPdf.cxx
ClassImp(RooQCDPdf) 

RooQCDPdf::RooQCDPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _p0,
                        RooAbsReal& _p1,
                        RooAbsReal& _p2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   p0("p0","p0",this,_p0),
   p1("p1","p1",this,_p1),
   p2("p2","p2",this,_p2){} 


RooQCDPdf::RooQCDPdf(const RooQCDPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   p0("p0",this,other.p0),
   p1("p1",this,other.p1),
   p2("p2",this,other.p2){} 



Double_t RooQCDPdf::evaluate() const { 

   Double_t sqrt_s=2000.;
   return TMath::Power(1-x/sqrt_s ,p0)/TMath::Power(x/sqrt_s, p1+p2*TMath::Log(x/sqrt_s))  ; 
} 

//////////////////////////////////////////RooUser1Pdf.cxx
ClassImp(RooUser1Pdf) 

RooUser1Pdf::RooUser1Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _p0,
                        RooAbsReal& _p1
                        ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   p0("p0","p0",this,_p0),
   p1("p1","p1",this,_p1){ } 


RooUser1Pdf::RooUser1Pdf(const RooUser1Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   p0("p0",this,other.p0),
   p1("p1",this,other.p1){} 



Double_t RooUser1Pdf::evaluate() const { 
   Double_t sqrt_s=500.;
   return TMath::Power(1-x/sqrt_s ,p0)/TMath::Power(x/sqrt_s, p1)  ; 
} 



///////////////////////////////////////////////RooExpNPdf.cxx
Double_t ExpN(Double_t x, Double_t c, Double_t n){
    return TMath::Exp( c*x+n/x ); 
}

///////////////////////////////////////////////RooExpN2Pdf.cxx
Double_t ExpN2(Double_t x, Double_t c, Double_t n, Double_t w){
    return TMath::Exp( c*x+n/x+w/TMath::Power(x,2.) ); 
}

ClassImp(RooExpNPdf) 

RooExpNPdf::RooExpNPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _n) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   n("n","n",this,_n){} 


RooExpNPdf::RooExpNPdf(const RooExpNPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   n("n",this,other.n){} 


Double_t RooExpNPdf::evaluate() const { 
   return ExpN(x,c,n); 
 } 


ClassImp(RooAlpha4ExpNPdf) 

RooAlpha4ExpNPdf::RooAlpha4ExpNPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _n0,
                        RooAbsReal& _c1,
                        RooAbsReal& _n1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   n0("n0","n0",this,_n0),
   c1("c1","c1",this,_c1),
   n1("n1","n1",this,_n1){} 


RooAlpha4ExpNPdf::RooAlpha4ExpNPdf(const RooAlpha4ExpNPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   n0("n0",this,other.n0),
   c1("c1",this,other.c1),
   n1("n1",this,other.n1){} 


Double_t RooAlpha4ExpNPdf::evaluate() const { 
   return ExpN(x, c0-c1, n0-n1); 
} 


///////////////////////////////////////////////RooErfExpNPdf.cxx
Double_t ErfExpN(Double_t x, Double_t c, Double_t n, Double_t o, Double_t w ){
  return ExpN(x,c,n)*((1.+TMath::Erf((x-o)/w))/2);
}

ClassImp(RooErfExpNPdf) 

RooErfExpNPdf::RooErfExpNPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _n,
                        RooAbsReal& _o,
                        RooAbsReal& _w ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   n("n","n",this,_n),
   o("o","o",this,_o),
   w("w","w",this,_w) {} 

RooErfExpNPdf::RooErfExpNPdf(const RooErfExpNPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   n("n",this,other.n),
   o("o",this,other.o),
   w("w",this,other.w){} 

Double_t RooErfExpNPdf::evaluate() const { 
   return ErfExpN(x, c, n, o, w);
 } 


ClassImp(RooAlpha4ErfExpNPdf) 

RooAlpha4ErfExpNPdf::RooAlpha4ErfExpNPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _n0,
                        RooAbsReal& _o0,
                        RooAbsReal& _w0,
                        RooAbsReal& _c1,
                        RooAbsReal& _n1,
                        RooAbsReal& _o1,
                        RooAbsReal& _w1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   n0("n0","n0",this,_n0),
   o0("o0","o0",this,_o0),
   w0("w0","w0",this,_w0),
   c1("c1","c1",this,_c1),
   n1("n1","n1",this,_n1),
   o1("o1","o1",this,_o1),
   w1("w1","w1",this,_w1){} 


RooAlpha4ErfExpNPdf::RooAlpha4ErfExpNPdf(const RooAlpha4ErfExpNPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   n0("n0",this,other.n0),
   o0("o0",this,other.o0),
   w0("w0",this,other.w0),
   c1("c1",this,other.c1),
   n1("n1",this,other.n1),
   o1("o1",this,other.o1),
   w1("w1",this,other.w1){} 


Double_t RooAlpha4ErfExpNPdf::evaluate() const { 
   return ErfExpN(x, c0, n0, o0, w0)/ErfExpN(x, c1, n1, o1, w1) ;
} 


///////////////////////////////////////////////RooExpTailPdf.cxx
Double_t ExpTail2(Double_t x, Double_t s, Double_t a, Double_t w){
    return TMath::Exp( -x/(s+a*x) +w/x ); 
}
 
ClassImp(RooExpTail2Pdf) 

RooExpTail2Pdf::RooExpTail2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s,
                        RooAbsReal& _a,
                        RooAbsReal& _w) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s("s","s",this,_s),
   a("a","a",this,_a),
   w("w","w",this,_w){} 


RooExpTail2Pdf::RooExpTail2Pdf(const RooExpTail2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s("s",this,other.s),
   a("a",this,other.a),
   w("w",this,other.w){} 


Double_t RooExpTail2Pdf::evaluate() const { 
   return ExpTail2(x, s, a, w) ; 
 } 


ClassImp(RooAlpha4ExpTail2Pdf) 

 RooAlpha4ExpTail2Pdf::RooAlpha4ExpTail2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s0,
                        RooAbsReal& _a0,
                        RooAbsReal& _w0,
                        RooAbsReal& _s1,
                        RooAbsReal& _a1,
                        RooAbsReal& _w1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s0("s0","s0",this,_s0),
   a0("a0","a0",this,_a0),
   w0("w0","w0",this,_w0),
   s1("s1","s1",this,_s1),
   a1("a1","a1",this,_a1),
   w1("w1","w1",this,_w1){} 


RooAlpha4ExpTail2Pdf::RooAlpha4ExpTail2Pdf(const RooAlpha4ExpTail2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s0("s0",this,other.s0),
   a0("a0",this,other.a0),
   w0("w0",this,other.w0),
   s1("s1",this,other.s1),
   a1("a1",this,other.a1),
   w1("w1",this,other.w1){} 



Double_t RooAlpha4ExpTail2Pdf::evaluate() const { 
   return ExpTail2(x, s0, a0, w0)/ExpTail2(x, s1, a1, w1) ; 
} 



///////////////////////////////////////////////RooExpTailPdf.cxx
Double_t ExpTail(Double_t x, Double_t s, Double_t a){
    return TMath::Exp( -x/(s+a*x) ); 
}
 
ClassImp(RooExpTailPdf) 

RooExpTailPdf::RooExpTailPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s,
                        RooAbsReal& _a) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s("s","s",this,_s),
   a("a","a",this,_a){} 


RooExpTailPdf::RooExpTailPdf(const RooExpTailPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s("s",this,other.s),
   a("a",this,other.a){} 


Double_t RooExpTailPdf::evaluate() const { 
   return ExpTail(x, s, a) ; 
 } 

ClassImp(RooAlpha4ExpTailPdf) 

 RooAlpha4ExpTailPdf::RooAlpha4ExpTailPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s0,
                        RooAbsReal& _a0,
                        RooAbsReal& _s1,
                        RooAbsReal& _a1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s0("s0","s0",this,_s0),
   a0("a0","a0",this,_a0),
   s1("s1","s1",this,_s1),
   a1("a1","a1",this,_a1){} 


RooAlpha4ExpTailPdf::RooAlpha4ExpTailPdf(const RooAlpha4ExpTailPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s0("s0",this,other.s0),
   a0("a0",this,other.a0),
   s1("s1",this,other.s1),
   a1("a1",this,other.a1){} 



Double_t RooAlpha4ExpTailPdf::evaluate() const { 
   return ExpTail(x, s0, a0)/ExpTail(x, s1, a1) ; 
} 


///////////////////////////////////////////////RooErfExpTailPdfSimple.cxx
ClassImp(RooExpN2Pdf) 

RooExpN2Pdf::RooExpN2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s,
                        RooAbsReal& _a,
                        RooAbsReal& _w) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s("s","s",this,_s),
   a("a","a",this,_a),
   w("w","w",this,_w) {} 


RooExpN2Pdf::RooExpN2Pdf(const RooExpN2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s("s",this,other.s),
   a("a",this,other.a),
   w("w",this,other.w) {} 


Double_t RooExpN2Pdf::evaluate() const { 
   return ExpN2(x, s, a, w) ; 
 } 

ClassImp(RooAlpha4ExpN2Pdf) 

 RooAlpha4ExpN2Pdf::RooAlpha4ExpN2Pdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s0,
                        RooAbsReal& _a0,
                        RooAbsReal& _w0,
                        RooAbsReal& _s1,
                        RooAbsReal& _a1,
                        RooAbsReal& _w1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s0("s0","s0",this,_s0),
   a0("a0","a0",this,_a0),
   w0("w0","w0",this,_w0),
   s1("s1","s1",this,_s1),
   a1("a1","a1",this,_a1),
   w1("w1","w1",this,_w1){} 


RooAlpha4ExpN2Pdf::RooAlpha4ExpN2Pdf(const RooAlpha4ExpN2Pdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s0("s0",this,other.s0),
   a0("a0",this,other.a0),
   w0("w0",this,other.w0),
   s1("s1",this,other.s1),
   a1("a1",this,other.a1),
   w1("w1",this,other.w1){} 



Double_t RooAlpha4ExpN2Pdf::evaluate() const { 
   return ExpN2(x, s0, a0, w0)/ExpN2(x, s1, a1, w1) ; 
} 



///////////////////////////////////////////////RooErfExpTailPdf.cxx
Double_t ErfExpTail(Double_t x, Double_t s, Double_t a, Double_t o, Double_t w){
  Double_t val = ExpTail(x,s,a)*((1.+TMath::Erf((x-o)/w))/2);
  return val ;
}
 
ClassImp(RooErfExpTailPdf) 

RooErfExpTailPdf::RooErfExpTailPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s,
                        RooAbsReal& _a,
                        RooAbsReal& _o,
                        RooAbsReal& _w ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s("s","s",this,_s),
   a("a","a",this,_a),
   o("o","o",this,_o),
   w("w","w",this,_w) {} 


RooErfExpTailPdf::RooErfExpTailPdf(const RooErfExpTailPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s("s",this,other.s),
   a("a",this,other.a),
   o("o",this,other.o),
   w("w",this,other.w){} 


Double_t RooErfExpTailPdf::evaluate() const { 
   return ErfExpTail(x, s, a, o, w) ; 
 } 

ClassImp(RooAlpha4ErfExpTailPdf) 

 RooAlpha4ErfExpTailPdf::RooAlpha4ErfExpTailPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _s0,
                        RooAbsReal& _a0,
                        RooAbsReal& _o0,
                        RooAbsReal& _w0,
                        RooAbsReal& _s1,
                        RooAbsReal& _a1,
                        RooAbsReal& _o1,
                        RooAbsReal& _w1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   s0("s0","s0",this,_s0),
   a0("a0","a0",this,_a0),
   o0("o0","o0",this,_o0),
   w0("w0","w0",this,_w0),
   s1("s1","s1",this,_s1),
   a1("a1","a1",this,_a1),
   o1("o1","o1",this,_o1),
   w1("w1","w1",this,_w1){} 


RooAlpha4ErfExpTailPdf::RooAlpha4ErfExpTailPdf(const RooAlpha4ErfExpTailPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   s0("s0",this,other.s0),
   a0("a0",this,other.a0),
   o0("o0",this,other.o0),
   w0("w0",this,other.w0),
   s1("s1",this,other.s1),
   a1("a1",this,other.a1),
   o1("o1",this,other.o1),
   w1("w1",this,other.w1){} 



Double_t RooAlpha4ErfExpTailPdf::evaluate() const { 
   return ErfExpTail(x, s0, a0, o0, w0)/ErfExpTail(x, s1, a1, o1, w1) ; 
} 





///////////////////////////////////////////////Roo2ExpPdf.cxx
Double_t TwoExp(Double_t x, Double_t c0, Double_t c1, Double_t frac){
	if(frac<0){frac=0.;}
	if(frac>1){frac=1.;}
	return TMath::Exp(x*c0)+frac*TMath::Exp(x*c1);
}

ClassImp(Roo2ExpPdf) 

Roo2ExpPdf::Roo2ExpPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c0,
                        RooAbsReal& _c1,
                        RooAbsReal& _frac) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c0("c0","c0",this,_c0),
   c1("c1","c1",this,_c1),
   frac("frac","frac",this,_frac){} 


Roo2ExpPdf::Roo2ExpPdf(const Roo2ExpPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c0("c0",this,other.c0),
   c1("c1",this,other.c1),
   frac("frac",this,other.frac){} 

Double_t Roo2ExpPdf::evaluate() const { 
   return TwoExp(x,c0,c1,frac);
} 

ClassImp(RooAlpha42ExpPdf) 

RooAlpha42ExpPdf::RooAlpha42ExpPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c00,
                        RooAbsReal& _c01,
                        RooAbsReal& _frac0,
                        RooAbsReal& _c10,
                        RooAbsReal& _c11,
                        RooAbsReal& _frac1) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c00("c00","c00",this,_c00),
   c01("c01","c01",this,_c01),
   frac0("frac0","frac0",this,_frac0),
   c10("c10","c10",this,_c10),
   c11("c11","c11",this,_c11),
   frac1("frac1","frac1",this,_frac1){} 


RooAlpha42ExpPdf::RooAlpha42ExpPdf(const RooAlpha42ExpPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c00("c00",this,other.c00),
   c01("c01",this,other.c01),
   frac0("frac0",this,other.frac0),
   c10("c10",this,other.c10),
   c11("c11",this,other.c11),
   frac1("frac1",this,other.frac1){} 



Double_t RooAlpha42ExpPdf::evaluate() const { 
   //return 1.0 ; 
   return TwoExp(x,c00,c01,frac0)/TwoExp(x,c10,c11,frac1);
} 



// RooAnaExpNPdf.cxx
ClassImp(RooAnaExpNPdf) 

 RooAnaExpNPdf::RooAnaExpNPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                        RooAbsReal& _n) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   n("n","n",this,_n){} 


 RooAnaExpNPdf::RooAnaExpNPdf(const RooAnaExpNPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   n("n",this,other.n){} 

 Double_t RooAnaExpNPdf::evaluate() const { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return ExpN(x,c,n) ; 
 } 

Int_t RooAnaExpNPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const { 

   if (matchArgs(allVars,analVars,x)){ return 1; } 
   return 0 ; 
 } 

Double_t gamma_in_mathematica(Double_t a, Double_t z){
	return (1-TMath::Gamma(a,z))*TMath::Gamma(a);
}

Double_t integral_ExpN(Double_t x, Double_t c, Double_t n){
	return -1*( TMath::Power(-1*c , -1/n ) * gamma_in_mathematica( 1/n, -1*c*TMath::Power(x, n) )   )/n   ;
}

Double_t RooAnaExpNPdf::analyticalIntegral(Int_t code, const char* rangeName) const  { 

	if (code==1) { 
		Double_t x_min=x.min(rangeName);
		Double_t x_max=x.max(rangeName);

        Double_t minTerm=integral_ExpN(x_min,c,n);
        Double_t maxTerm=integral_ExpN(x_max,c,n);
		cout<<"maxTerm-minTerm="<<maxTerm-minTerm<<endl;
		return (maxTerm-minTerm) ;
	} 
	return 0 ; 
} 


RooDoubleCrystalBall::RooDoubleCrystalBall(){}

RooDoubleCrystalBall::RooDoubleCrystalBall(const char *name, const char *title, 
			 RooAbsReal& _x,
			 RooAbsReal& _mean,
			 RooAbsReal& _width,
			 RooAbsReal& _alpha1,
			 RooAbsReal& _n1,
			 RooAbsReal& _alpha2,
			     RooAbsReal& _n2
			 ) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  width("width","width",this,_width),
  alpha1("alpha1","alpha1",this,_alpha1),
  n1("n1","n1",this,_n1),
  alpha2("alpha2","alpha2",this,_alpha2),
  n2("n2","n2",this,_n2){} 


RooDoubleCrystalBall::RooDoubleCrystalBall(const RooDoubleCrystalBall& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mean("mean",this,other.mean),
  width("width",this,other.width),
  alpha1("alpha1",this,other.alpha1),
  n1("n1",this,other.n1),
  alpha2("alpha2",this,other.alpha2),
  n2("n2",this,other.n2){} 

double RooDoubleCrystalBall::evaluate() const { 
  double t = (x-mean)/width;
  if(t>-alpha1 && t<alpha2){
    return exp(-0.5*t*t);
  }else if(t<-alpha1){
    double A1 = pow(n1/fabs(alpha1),n1)*exp(-alpha1*alpha1/2);
    double B1 = n1/fabs(alpha1)-fabs(alpha1);
    return A1*pow(B1-t,-n1);
  }else if(t>alpha2){
    double A2 = pow(n2/fabs(alpha2),n2)*exp(-alpha2*alpha2/2);
    double B2 = n2/fabs(alpha2)-fabs(alpha2);
    return A2*pow(B2+t,-n2);
  }else{
    cout << "ERROR evaluating range..." << endl;
    return 99;
  }
   
} 

Int_t RooDoubleCrystalBall::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const {
  if (matchArgs(allVars,analVars,x)) return 1;
  return 0;
}

Double_t RooDoubleCrystalBall::analyticalIntegral(Int_t code, const char* rangeName) const {
  assert(code==1) ;
 
  double central=0;
  double left=0;
  double right=0;
 
  static const Double_t root2 = sqrt(2) ;
  static const Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
  Double_t xscale = root2*width;
 
  //compute gaussian contribution
  double central_low =max(x.min(rangeName),mean - alpha1*width );
  double central_high=min(x.max(rangeName),mean + alpha2*width );
  if(central_low < central_high) // is the gaussian part in range?
    central = rootPiBy2*width*(TMath::Erf((central_high-mean)/xscale)-TMath::Erf((central_low-mean)/xscale));
 
  //compute left tail;
  double A1 = pow(n1/fabs(alpha1),n1)*exp(-alpha1*alpha1/2);
  double B1 = n1/fabs(alpha1)-fabs(alpha1);
 
  double left_low=x.min(rangeName);
  double left_high=min(x.max(rangeName),mean - alpha1*width);
  if(left_low < left_high){ //is the left tail in range?
    if(fabs(n1-1.0)>1.e-5)
      left = A1/(-n1+1.0)*width*(pow(B1-(left_low-mean)/width,-n1+1.)-pow(B1-(left_high-mean)/width,-n1+1.));
    else
      left = A1*width*(log(B1-(left_low-mean)/width) - log(B1-(left_high-mean)/width) );
  }
 
  //compute right tail;
  double A2 = pow(n2/fabs(alpha2),n2)*exp(-alpha2*alpha2/2);
  double B2 = n2/fabs(alpha2)-fabs(alpha2);
 
  double right_low=max(x.min(rangeName),mean + alpha2*width);
  double right_high=x.max(rangeName);
  if(right_low < right_high){ //is the right tail in range?
    if(fabs(n2-1.0)>1.e-5)
      right = A2/(-n2+1.0)*width*(pow(B2+(right_high-mean)/width,-n2+1.)-pow(B2+(right_low-mean)/width,-n2+1.));
    else
      right = A2*width*(log(B2+(right_high-mean)/width) - log(B2+(right_low-mean)/width) );
  }
     
  return left+central+right;
 
}




