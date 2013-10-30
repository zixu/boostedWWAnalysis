#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import subprocess
from array import array
from subprocess import Popen

from ROOT import gROOT, gStyle, gSystem, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooSimultaneous, RooGenericPdf, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C")

from ROOT import setTDRStyle
from ROOT import RooTrace

ROOT.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPalette(1);

ROOT.gSystem.Load("PDFs/RooRelBWRunningWidth_cxx.so")

############################################################

def FitMassPoint(filename, massin, massmin, massmax, nbins=50):
    
    inputFile = ROOT.TFile(filename);
    tree = inputFile.Get("WJet");
    
    print "n entries: ", tree.GetEntries();
    
    # RooFitting
    
    rrv_mass   = ROOT.RooRealVar("rrv_mass","rrv_mass",massmin,massmax)
    rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 

    rrv_mH2    = ROOT.RooRealVar("rrv_mH2","rrv_mH2", massin, massmin, massmax )
    rrv_gamma2 = ROOT.RooRealVar("rrv_gamma2","rrv_gamma2",20.,massmax)

    rds_raw = ROOT.RooDataSet("rds_raw","rds_raw",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight)); ## created the raw dataset
    rds_cps = ROOT.RooDataSet("rds_cps","rds_cps",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight)); ## create the cps dataset
    rds_cps_intf = ROOT.RooDataSet("rds_cps_intf","rds_cps_intf",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight)); ## create the cps + interference dataset

    model2_pdf = ROOT.RooRelBWRunningWidth("model2_pdf","model2_pdf",rrv_mass,rrv_mH2,rrv_gamma2); ## BWRunningWidth Pdf
    
    ### loop on Higgs signal events 
    for i in range(tree.GetEntries()):
        
        if i % 10000 == 0: print "reweighter, i: ", i;
        tree.GetEntry(i);
        
        curmass = getattr(tree,"W_H_mass_gen"); ## take the higgs generated mass inside the window
        if curmass < massmax and curmass > massmin:
            
            rrv_mass.setVal( curmass );            
            tmpweight_cps = getattr(tree,"complexpolewtggH"+str(massin))/getattr(tree,"avecomplexpolewtggH"+str(massin));
            tmpweight_cps_intf = getattr(tree,"complexpolewtggH"+str(massin))*getattr(tree,"interferencewtggH"+str(massin))/getattr(tree,"avecomplexpolewtggH"+str(massin));    

            rds_raw.add( RooArgSet( rrv_mass ), 1. );
            rds_cps.add( RooArgSet( rrv_mass ), tmpweight_cps );
            rds_cps_intf.add( RooArgSet( rrv_mass ), tmpweight_cps_intf );


    print ">>>>"
    RooTrace.dump(ROOT.cout,ROOT.kTRUE);
    RooTrace.mark();                
    print "<<<<"

    model2_pdf.fitTo(rds_cps,RooFit.Save(0), RooFit.SumW2Error(kTRUE) );

#    mplot = rrv_mass.frame(RooFit.Title("mass plot"));
#    rds_raw.plotOn(mplot, RooFit.MarkerColor(kBlack), RooFit.LineColor(kBlack), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
#    rds_cps_intf.plotOn(mplot, RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
#    rds_cps.plotOn(mplot, RooFit.MarkerColor(kRed), RooFit.LineColor(kRed), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
#    model_pdf.plotOn(mplot, RooFit.LineColor(kRed));
#    model2_pdf.plotOn(mplot, RooFit.LineColor(kRed), RooFit.LineStyle(2) );
#    model3_pdf.plotOn(mplot, RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(2) );
#    model4_pdf.plotOn(mplot, RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(3) );

    print "rds_cps_intf.sumEntries() = ", rds_cps_intf.sumEntries()
    print "model2_pdf: mH = ", rrv_mH2.getVal(), ", gamma = ", rrv_gamma2.getVal();    
#
#    dummy_h1 = ROOT.TH1F("dummy_h1","dummy_h1",1,0,1); 
#    dummy_h1.SetMarkerColor( ROOT.kBlack );
#    dummy_h2 = ROOT.TH1F("dummy_h2","dummy_h2",1,0,1); 
#    dummy_h2.SetMarkerColor( ROOT.kRed );
#    dummy_h3 = ROOT.TH1F("dummy_h3","dummy_h3",1,0,1); 
#    dummy_h3.SetMarkerColor( ROOT.kBlue );
#    dummy_h5 = ROOT.TH1F("dummy_h5","dummy_h5",1,0,1); 
#    dummy_h5.SetLineColor( ROOT.kRed ); dummy_h5.SetLineStyle( 2 );
#
#    L = TLegend(0.65,0.60,0.93,0.85);
#    L.SetFillStyle(0);
#    L.AddEntry(dummy_h1,"Powheg","p");
#    L.AddEntry(dummy_h2,"w/CPS weight","p");
#    L.AddEntry(dummy_h3,"w/CPS,Intf weight","p");
#    L.AddEntry(dummy_h5,"Fit, BW (running)","l");
#
#    can2 = ROOT.TCanvas("can2","can2",800,800);
#    mplot.Draw();
#    L.Draw();
#    can2.SaveAs("massFits/mass_rf_"+str(massin)+".eps");
#    can2.SaveAs("massFits/mass_rf_"+str(massin)+".png");

    outputpar = [];
    outputpar.append( rrv_mH2.getVal() );
    outputpar.append( rrv_gamma2.getVal() );


    model2_pdf.Delete();
    rds_raw.Delete(),
    rds_cps.Delete(),
    rds_cps_intf.Delete();
    rrv_mH2.Delete(),
    rrv_gamma2.Delete();
    rrv_mass.Delete();

    return outputpar

### definition of the BW relativistic
def MyRunningWidthBW(mww, mH, gamma):

    numerator = (mww*mww*gamma/mH) ;
    denominator = (mww*mww-mH*mH)*(mww*mww-mH*mH) + (mww*mww*gamma/mH)*(mww*mww*gamma/mH);
    pdf = numerator/denominator;
    return pdf;


def lineshapeWidthReweight(point, meanSM, gammaSM, Cprime, massmin, massmax):

    
    weight = MyRunningWidthBW(point,meanSM,gammaSM*Cprime)/MyRunningWidthBW(point,meanSM,gammaSM);
    return weight;

def IntfRescale(curIntfRw,cPrime,BRnew):

    curIoverS = curIntfRw - 1;
    newWeight = 1 + ((1-BRnew)*curIoverS)/cPrime;
    if newWeight < 0: newWeight = 0.01;
    ratio = newWeight/curIntfRw;
        
    return ratio;
