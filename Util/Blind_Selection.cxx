/*
 * =====================================================================================
 *
 *       Filename:  Blind_Selection.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/25/2013 12:20:31 PM CST
 *       Revision:  none
 *       Compiler:  gcc, root
 *
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
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


void Blind_Selection(char* channel){
	Double_t lumi_ful=0;
	Double_t lumi_sub=0.;
	if (TString(channel)=="el"){ lumi_ful=19.2;lumi_sub=5.1;}
	if (TString(channel)=="mu"){ lumi_ful=19.3;lumi_sub=5.3;}
	Double_t selection_scale=lumi_sub/lumi_ful;


	TString datafile = Form("ofile_data.root",channel);
	//TString datafile = Form("../trainingtrees_19_Jan24/trainingtrees_%s/ofile_data.root",channel);
   //Get old file, old tree and set top branch address
   TFile *oldfile = new TFile(datafile);
   TTree *oldtree = (TTree*)oldfile->Get("otree");
   Long64_t nentries = oldtree->GetEntries();
   Double_t number_newtree=0;

   //Create a new file + a clone of old tree in new file
   //TFile *newfile = new TFile(Form("ofile_data_%s_%g.root", lumi_sub,channel),"recreate");
   TFile *newfile = new TFile(Form("ofile_data_sub.root"),"recreate");
   TTree *newtree = oldtree->CloneTree(0);

   TRandom3 rand(1234);

   for (Long64_t i=0;i<nentries; i++){
      oldtree->GetEntry(i);
      if ( rand.Uniform(1)<=selection_scale){
		  newtree->Fill();
		  number_newtree=number_newtree+1;
	  }
   }
   cout<<"number_newtree="<<number_newtree<<" lumi_sub="<<number_newtree/nentries*lumi_ful<<endl;
   newtree->Print();
   newtree->AutoSave();
   delete oldfile;
   delete newfile;
}

