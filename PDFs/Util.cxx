/*
 * =====================================================================================
 *
 *       Filename:  Util.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/26/2012 03:27:25 PM CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   Zijun Xu (xuzijun123@gmail.com), 
 *        Company:  
 *
 * =====================================================================================
 */


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

void draw_error_band(RooAbsData &rdata, RooAbsPdf &rpdf, RooRealVar &rrv_number_events , RooFitResult *rfres, RooPlot *mplot, Int_t kcolor=6,char* opt="F", Int_t number_point=100, const Int_t number_errorband=2000)
{
	TRandom3 rand(1234);

	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
    rpdf.Print("v");
    rpdf.getParameters(RooArgSet(*rrv_x))->Print("v");
	rrv_number_events.Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	//Double_t width_x=rrv_x->getBinWidth(1);
	Double_t width_x=mplot->getFitRangeBinW();
	//rdata.plotOn(mplot);
	//rpdf.plotOn(mplot,RooFit::VisualizeError(*rfres,1),RooFit::FillColor(kOrange));
	//rpdf.plotOn(mplot);

	Double_t number_events_mean = rrv_number_events.getVal();
	Double_t number_events_sigma= rrv_number_events.getError();

	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rrv_number_events.getVal()*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	//par_pdf->Print("v");

	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		syst[j]=new TGraph(number_point+1);
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;
		Double_t number_events_tmp = rand.Gaus(number_events_mean,number_events_sigma);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , number_events_tmp*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	RooArgList par_tmp = rfres->floatParsFinal();
	*par_pdf = par_tmp;

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	//errorband->SetFillColor(kcolor);
	errorband->SetFillColor(kBlack);
    errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ){ mplot->addObject(errorband,"E3"); //mplot->addObject(bkgpred);
            }
	if( TString(opt).Contains("L") ){ mplot->addObject(am); mplot->addObject(ap); }
	//mplot->addObject(errorband,"E3"); 
	//mplot->addObject(am); 
	//mplot->addObject(ap); 
}
void draw_error_band( RooAbsPdf &rpdf, char* xaxis_name, RooRealVar &rrv_number_events , RooArgList &paras, RooWorkspace &ws, RooPlot *mplot, Int_t kcolor=6,char* opt="F", Int_t number_point=100, const Int_t number_errorband=2000)
{
	TRandom3 rand(1234);

	RooRealVar *rrv_x=ws.var(xaxis_name);
    rpdf.Print("v");
    rpdf.getParameters(RooArgSet(*rrv_x))->Print("v");
	rrv_number_events.Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	//Double_t width_x=rrv_x->getBinWidth(1);
	Double_t width_x=mplot->getFitRangeBinW();

	Double_t number_events_mean = rrv_number_events.getVal();
	Double_t number_events_sigma= rrv_number_events.getError();

	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rrv_number_events.getVal()*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
            ws.var(paras[ipara].GetName())->setConstant(0);
	        ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,ws.var(paras[ipara].GetName())->getError()) );
            //ws.var(paras[ipara].GetName())->Print();
		}
    //{double tmpb;cout<<"tmpb";cin>>tmpb;}

		Double_t number_events_tmp = rand.Gaus(number_events_mean,number_events_sigma);
		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , number_events_tmp*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	//errorband->SetFillColor(kcolor);
	errorband->SetFillColor(kBlack);
    errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ){
            mplot->addObject(errorband,"E3");
            //mplot->addObject(bkgpred);
            }
	if( TString(opt).Contains("L") ){ mplot->addObject(am); mplot->addObject(ap); }
}


void draw_error_band_extendPdf(RooAbsData &rdata, RooExtendPdf &rpdf, RooFitResult *rfres, RooPlot *mplot, Int_t kcolor=6,char* opt="F", Int_t number_point=100, const Int_t number_errorband=2000)
{
	TRandom3 rand(1234);

	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	//Double_t width_x=rrv_x->getBinWidth(1);
	Double_t width_x=mplot->getFitRangeBinW();

	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , rpdf.expectedEvents(*rrv_x)*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	par_pdf->Print("v");

	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		syst[j]=new TGraph(number_point+1);
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , rpdf.expectedEvents(*rrv_x)*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	RooArgList par_tmp = rfres->floatParsFinal();
	*par_pdf = par_tmp;

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	//errorband->SetFillColor(kcolor);
	errorband->SetFillColor(kBlack);
    errorband->SetFillStyle(3013);

	if( TString(opt).Contains("F") ){ mplot->addObject(errorband,"E3"); }
	if( TString(opt).Contains("L") ){ mplot->addObject(am); mplot->addObject(ap); }
	//mplot->addObject(errorband,"E3"); 
	//mplot->addObject(am); 
	//mplot->addObject(ap); 
}

void draw_error_band2(RooAbsData &rdata, RooAbsPdf &rpdf, RooRealVar &rrv_number_events , RooFitResult *rfres, RooPlot *mplot, Int_t kcolor=9,char* opt="F", Int_t number_point=100, const Int_t number_errorband=2000)//don't think of covariance-matrix
{
	TRandom3 rand(1234);

	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
	rrv_number_events.Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	//Double_t width_x=rrv_x->getBinWidth(1);
	Double_t width_x=mplot->getFitRangeBinW();

	Double_t number_events_mean = rrv_number_events.getVal();
	Double_t number_events_sigma= rrv_number_events.getError();

	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		//cout<<rpdf.getVal(*rrv_x)<<endl;
		bkgpred->SetPoint( i , x_min+delta_x*i , rrv_number_events.getVal()*rpdf.getVal(*rrv_x)*width_x );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor);

	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	par_pdf->Print("v");

	TGraph* syst[number_errorband];
	RooArgList par_ParsFinal = rfres->floatParsFinal();
	for(int j=0;j<number_errorband;j++){
		syst[j]=new TGraph(number_point+1);
		//RooArgList par_tmp = rfres->randomizePars();
		//*par_pdf = par_tmp;
        TIterator *par=par_pdf->createIterator();
        TIterator *par_final=par_ParsFinal.createIterator();
        par->Reset();
        par_final->Reset();
        RooRealVar* param=(RooRealVar*)par->Next();
        RooRealVar* param_final=(RooRealVar*)par_final->Next();
        while(param){
            param->setVal(rand.Gaus(param_final->getVal(),param_final->getError() ) );
            param=(RooRealVar*)par->Next();
            param_final=(RooRealVar*)par_final->Next();
             
        }
		Double_t number_events_tmp = rand.Gaus(number_events_mean,number_events_sigma);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i , number_events_tmp*rpdf.getVal(*rrv_x)*width_x);
		}
	}

	*par_pdf = par_ParsFinal;

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	//errorband->SetFillColor(kcolor);
	errorband->SetFillColor(kBlack);
    errorband->SetFillStyle(3013);
    errorband->SetName("Uncertainty");

	if( TString(opt).Contains("F") ){ mplot->addObject(errorband,"E3"); }
	if( TString(opt).Contains("L") ){ mplot->addObject(am); mplot->addObject(ap); }
	//mplot->addObject(errorband,"E3"); 
	//mplot->addObject(am); 
	//mplot->addObject(ap); 
}
void draw_error_band_Decor( char* pdf_name, char* xaxis_name, RooArgList &paras, RooWorkspace &ws,RooRealVar &rrv_shape_scale , RooPlot *mplot, Int_t kcolor=6,char* opt="F", Int_t number_point=100, const Int_t number_errorband=2000)
{
	TRandom3 rand(1234);

	RooRealVar *rrv_x=ws.var(xaxis_name);
	rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	//Double_t width_x=rrv_x->getBinWidth(1);
	Double_t width_x=mplot->getFitRangeBinW();

	Double_t shape_scale = rrv_shape_scale.getVal();
	Double_t shape_scale_error = rrv_shape_scale.getError();

	TGraph *bkgpred=new TGraph(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , shape_scale*ws.pdf(pdf_name)->getVal() );
	}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor+3);

	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
	        ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,1.) );
		}

        Double_t shape_scale_tmp=rand.Gaus(shape_scale,shape_scale_error);
		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i ,shape_scale_tmp*ws.pdf(pdf_name)->getVal());
		}
	}

	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	//errorband->SetFillColor(kcolor);
	errorband->SetFillColor(kBlack);
    errorband->SetFillStyle(3013);
    errorband->SetName("Uncertainty");

	//mplot->addObject(errorband,"E3"); 
	//mplot->addObject(bkgpred); 
	//mplot->addObject(am); 
	//mplot->addObject(ap); 

    if( TString(opt).Contains("F") ){
            mplot->addObject(errorband,"E3"); 
            mplot->addObject(bkgpred);
            }
	if( TString(opt).Contains("L") ){
            mplot->addObject(am); mplot->addObject(ap); 
            mplot->addObject(bkgpred);
            }

	for(Int_t ipara=0;ipara<paras.getSize();ipara++){
		ws.var(paras[ipara].GetName())->setVal(0.);
	}

} 


void draw_error_band_shape_Decor( char* pdf_name, char* xaxis_name, RooArgList &paras, RooWorkspace &ws,Double_t sigma , RooPlot *mplot, Int_t kcolor=6,char* opt="F", Int_t fillstyle=3013,char* uncertainty_title="", Int_t number_point=100, const Int_t number_errorband=2000)
{
	TRandom3 rand(1234);

	RooRealVar *rrv_x=ws.var(xaxis_name);
	rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	Double_t delta_x=(x_max-x_min)/number_point;
	//Double_t width_x=rrv_x->getBinWidth(1);
	Double_t width_x=mplot->getFitRangeBinW();

	TGraph *bkgpred=new TGraph(number_point+1);
    //{double tmpb;cout<<"tmpb";cin>>tmpb;}
	for(int i =0 ; i<= number_point ; i++){
		rrv_x->setVal(x_min+delta_x*i); 
		bkgpred->SetPoint( i , x_min+delta_x*i , ws.pdf(pdf_name)->getVal(*rrv_x)*width_x );
	}
    //{double tmpb;cout<<"tmpb";cin>>tmpb;}
	bkgpred->SetLineWidth(2);
	bkgpred->SetLineColor(kcolor+3);

	TGraph* syst[number_errorband];
	for(int j=0;j<number_errorband;j++){
        //cout<<"random "<<j<<"---------------------------------------------------"<<endl;
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
	        ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,sigma) );
            //cout<<paras[ipara].GetName()<<endl;
	        //ws.var(paras[ipara].GetName())->Print("");
		}

		syst[j]=new TGraph(number_point+1);
		for(int i =0 ; i<=number_point ; i++){
			rrv_x->setVal(x_min+delta_x*i); 
			syst[j]->SetPoint( i , x_min+delta_x*i ,ws.pdf(pdf_name)->getVal(*rrv_x)*width_x);
		}
	}
	std::vector<double> val;
	val.resize(number_errorband);
	TGraph *ap=new TGraph(number_point+1);
	TGraph *am=new TGraph(number_point+1);
	TGraphAsymmErrors* errorband=new TGraphAsymmErrors(number_point+1);
	for(int i =0 ; i<= number_point ; i++){
		for(int j=0;j<number_errorband;j++){
			val[j]=(syst[j])->GetY()[i];
		}
		std::sort(val.begin(),val.end());
		ap->SetPoint(i, x_min+delta_x*i,val[Int_t(0.16*number_errorband)]);
		am->SetPoint(i, x_min+delta_x*i,val[Int_t(0.84*number_errorband)]);
		errorband->SetPoint(i, x_min+delta_x*i,bkgpred->GetY()[i] );
		errorband->SetPointError(i, 0.,0., bkgpred->GetY()[i]-val[Int_t(0.84*number_errorband)],val[Int_t(0.16*number_errorband)]-bkgpred->GetY()[i]);
        if (i==55)cout<<val[Int_t(0.16*number_errorband)]<<"  "<<bkgpred->GetY()[i] <<"  "<<val[Int_t(0.84*number_errorband)] <<endl;
	}
	ap->SetLineWidth(2);
	ap->SetLineColor(kcolor);
	am->SetLineWidth(2);
	am->SetLineColor(kcolor);
	//errorband->SetFillColor(kBlack+7*(sigma-1));
	errorband->SetFillColor(kcolor);
    //errorband->SetFillStyle(3013+sigma-1);
    errorband->SetFillStyle(fillstyle);
    //errorband->SetName(Form("Uncertainty of %g#sigma",sigma));
    errorband->SetName(Form("%s %g#sigma",uncertainty_title,sigma));

    if( TString(opt).Contains("F") ){ mplot->addObject(errorband,"E3"); }
	if( TString(opt).Contains("L") ){ mplot->addObject(am); mplot->addObject(ap); }
	for(Int_t ipara=0;ipara<paras.getSize();ipara++){
		ws.var(paras[ipara].GetName())->setVal(0.);
	}
} 



double Calc_error_extendPdf(RooAbsData &rdata, RooExtendPdf &rpdf, RooFitResult *rfres, char* range, const Int_t calc_times=2000)
{
	TRandom3 rand(1234);

	RooArgSet* argset_obs=rpdf.getObservables(rdata);
	TIterator *par=argset_obs->createIterator();
	par->Reset();
	RooRealVar *rrv_x=(RooRealVar*)par->Next();
	rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();

    RooAbsReal* fullInt = rpdf.createIntegral(*rrv_x,*rrv_x);
    RooAbsReal* signalInt = rpdf.createIntegral(*rrv_x,*rrv_x,range);
    double fullInt_var=fullInt->getVal();
    double signalInt_var=signalInt->getVal()/fullInt_var;

    double signal_number_media=signalInt_var*rpdf.expectedEvents(*rrv_x);
    cout<<"signal_number_media="<<signal_number_media<<endl;


	RooArgSet* par_pdf  = rpdf.getParameters(RooArgSet(*rrv_x)) ;
	par_pdf->Print("v");

	std::vector<double> val;
	val.resize(calc_times);
	for(int j=0;j<calc_times;j++){
		RooArgList par_tmp = rfres->randomizePars();
		*par_pdf = par_tmp;

        double fullInt_var=fullInt->getVal();
        double signalInt_var=signalInt->getVal()/fullInt_var;
        //cout<<"signalInt_var="<<signalInt_var<<endl;

        double signal_number_media=signalInt_var*rpdf.expectedEvents(*rrv_x);
        val[j]=signal_number_media;
	}

	std::sort(val.begin(),val.end());
    //cout<<"signal_number_low_error="<<signal_number_media-val[Int_t(0.16*calc_times)];
    //cout<<"signal_number_hig_error="<<val[Int_t(0.84*calc_times)]-signal_number_media;
    return (val[Int_t(0.84*calc_times)]-val[Int_t(0.16*calc_times)])/2.;
}

double Calc_error( char* rpdfname, char* xaxis_name , RooArgList &paras, RooWorkspace &ws,char*range, const Int_t calc_times=2000){
	TRandom3 rand(1234);

	RooRealVar *rrv_x=ws.var(xaxis_name); rrv_x->Print();
	Double_t x_min=rrv_x->getMin();
	Double_t x_max=rrv_x->getMax();
	RooAbsPdf*rpdf=ws.pdf(rpdfname);

    RooAbsReal* fullInt = rpdf->createIntegral(*rrv_x,*rrv_x);
    RooAbsReal* signalInt = rpdf->createIntegral(*rrv_x,*rrv_x,range);
    double fullInt_var=fullInt->getVal();
    double signalInt_var=signalInt->getVal()/fullInt_var;

    double signal_number_media=signalInt_var;
    //cout<<"signal_number_media="<<signal_number_media<<endl;


	std::vector<double> val;
	val.resize(calc_times);
	for(int j=0;j<calc_times;j++){
		for(Int_t ipara=0;ipara<paras.getSize();ipara++){
            ws.var(paras[ipara].GetName())->setConstant(0);
	        ws.var(paras[ipara].GetName())->setVal( rand.Gaus(0.,ws.var(paras[ipara].GetName())->getError()) );
		}

        double fullInt_var=fullInt->getVal();
        double signalInt_var=signalInt->getVal()/fullInt_var;

        double signal_number_media=signalInt_var;
        val[j]=signal_number_media;
        //cout<<"signal_number_tmp="<<signal_number_media<<endl;
	}
	for(Int_t ipara=0;ipara<paras.getSize();ipara++){ ws.var(paras[ipara].GetName())->setVal(0.); }

	std::sort(val.begin(),val.end());
	double number_error=(val[Int_t(0.84*calc_times)]-val[Int_t(0.16*calc_times)])/2./signal_number_media;
    return number_error;
}


TH1* makeTH1() 
{
	// Create ROOT TH1 filled with a Gaussian distribution

	TH1D* hh = new TH1D("hh","hh",25,-10,10) ;
	for (int i=0 ; i<100 ; i++) {
		hh->Fill(gRandom->Gaus(0,3)) ;
	}
	return hh ;
}


TTree* makeTTree() 
{
	// Create ROOT TTree filled with a Gaussian distribution in x and a uniform distribution in y

	TTree* tree = new TTree("tree","tree") ;
	Double_t* px = new Double_t ;
	tree->Branch("x",px,"x/D") ;
	for (int i=0 ; i<100 ; i++) {
		*px = gRandom->Gaus(0,3) ;
		tree->Fill() ;
	}
	return tree ;
}


