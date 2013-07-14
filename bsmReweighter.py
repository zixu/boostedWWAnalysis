#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooSimultaneous, RooGenericPdf, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan
from array import array
from ROOT import gROOT, gStyle, gSystem, TLatex
import subprocess
from subprocess import Popen

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C")
from ROOT import setTDRStyle
from ROOT import RooTrace
ROOT.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPalette(1);

ROOT.gSystem.Load("PDFs/RooRelBWRunningWidth_cxx.so")

############################################################
############################################################

def FitMassPoint(filename, massin, massmin, massmax, nbins=50):

#    massin = 800;
#    massmin = 400;
#    massmax = 1400;
#    nbins = 50;
#    massin = 600;
#    massmin = 200;
#    massmax = 1000;
#    nbins = 40;
    
    inputFile = ROOT.TFile(filename);
    tree = inputFile.Get("WJet");
#    tree.Print();
    
    print "n entries: ", tree.GetEntries();
    
    ############################################################
    # RooFitting
    
    rrv_mass = ROOT.RooRealVar("rrv_mass","rrv_mass",massmin,massmax)
    rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 

#    rrv_mH = ROOT.RooRealVar("rrv_mH","rrv_mH", massin, massmin, massmax )
#    rrv_gamma = ROOT.RooRealVar("rrv_gamma","rrv_gamma",20.,massmax)
    rrv_mH2 = ROOT.RooRealVar("rrv_mH2","rrv_mH2", massin, massmin, massmax )
    rrv_gamma2 = ROOT.RooRealVar("rrv_gamma2","rrv_gamma2",20.,massmax)
#    rrv_mH3 = ROOT.RooRealVar("rrv_mH3","rrv_mH3", massin, massmin, massmax )
#    rrv_gamma3 = ROOT.RooRealVar("rrv_gamma3","rrv_gamma3",20.,massmax)
#    rrv_mH4 = ROOT.RooRealVar("rrv_mH4","rrv_mH4", massin, massmin, massmax )
#    rrv_gamma4 = ROOT.RooRealVar("rrv_gamma4","rrv_gamma4",20.,massmax)

    rds_raw = ROOT.RooDataSet("rds_raw","rds_raw",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_cps = ROOT.RooDataSet("rds_cps","rds_cps",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_cps_intf = ROOT.RooDataSet("rds_cps_intf","rds_cps_intf",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight));

    model2_pdf = ROOT.RooRelBWRunningWidth("model2_pdf","model2_pdf",rrv_mass,rrv_mH2,rrv_gamma2);

    for i in range(tree.GetEntries()):
        
        if i % 10000 == 0: print "reweighter, i: ", i;
        tree.GetEntry(i);
        
        curmass = getattr(tree,"W_H_mass_gen");
        if curmass < massmax and curmass > massmin:
            
            rrv_mass.setVal( curmass );
            
            tmpweight_cps = getattr(tree,"complexpolewtggH"+str(massin))/getattr(tree,"avecomplexpolewtggH"+str(massin));
            tmpweight_cps_intf = getattr(tree,"complexpolewtggH"+str(massin))*getattr(tree,"interferencewtggH"+str(massin))/getattr(tree,"avecomplexpolewtggH"+str(massin));    

            rds_raw.add( RooArgSet( rrv_mass ), 1. )
            rds_cps.add( RooArgSet( rrv_mass ), tmpweight_cps )
            rds_cps_intf.add( RooArgSet( rrv_mass ), tmpweight_cps_intf )

    ############################################################
    print ">>>>"
    RooTrace.dump(ROOT.cout,ROOT.kTRUE);
    RooTrace.mark();                
    print "<<<<"

#    rfresult2 = model2_pdf.fitTo(rds_cps,RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
    model2_pdf.fitTo(rds_cps,RooFit.Save(0), RooFit.SumW2Error(kTRUE) );

    print ">>>>2"
    RooTrace.dump(ROOT.cout,ROOT.kTRUE);
    RooTrace.mark();                
    print "<<<<2"

#    mplot = rrv_mass.frame(RooFit.Title("mass plot"));
#    rds_raw.plotOn(mplot, RooFit.MarkerColor(kBlack), RooFit.LineColor(kBlack), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
#    rds_cps_intf.plotOn(mplot, RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
#    rds_cps.plotOn(mplot, RooFit.MarkerColor(kRed), RooFit.LineColor(kRed), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
##    model_pdf.plotOn(mplot, RooFit.LineColor(kRed));
#    model2_pdf.plotOn(mplot, RooFit.LineColor(kRed), RooFit.LineStyle(2) );
##    model3_pdf.plotOn(mplot, RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(2) );
##    model4_pdf.plotOn(mplot, RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(3) );
#
##    print "rds_cps_intf.sumEntries() = ", rds_cps_intf.sumEntries()
##    print "model_pdf: mH = ", rrv_mH.getVal(), ", gamma = ", rrv_gamma.getVal();
##    print "model2_pdf: mH = ", rrv_mH2.getVal(), ", gamma = ", rrv_gamma2.getVal();    
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

    print ">>>>2.5"
    RooTrace.dump(ROOT.cout,ROOT.kTRUE);
    RooTrace.mark();                
    print "<<<<2.5"    

#    mplot.Delete();
#    rfresult2.Delete();
    model2_pdf.Delete();
    rds_raw.Delete(), rds_cps.Delete(), rds_cps_intf.Delete();
    rrv_mH2.Delete(), rrv_gamma2.Delete();
    rrv_mass.Delete();
#    del model2_pdf;
#    del rds_raw, rds_cps, rds_cps_intf;
#    del rrv_mH2, rrv_gamma2;

    print ">>>>2.9"
    RooTrace.dump(ROOT.cout,ROOT.kTRUE);
    RooTrace.mark();                
    print "<<<<2.9"    

    return outputpar

def MyRunningWidthBW(mww, mH, gamma):

    numerator = (mww*mww*gamma/mH) ;
    denominator = (mww*mww-mH*mH)*(mww*mww-mH*mH) + (mww*mww*gamma/mH)*(mww*mww*gamma/mH);
    pdf = numerator/denominator;
    return pdf;


def lineshapeWidthReweight(point, meanSM, gammaSM, Cprime, massmin, massmax):

#    rrv2_mass = ROOT.RooRealVar("rrv2_mass","rrv2_mass",massmin,massmax)
#
#    rrv2_mH1 = ROOT.RooRealVar("rrv2_mH1","rrv2_mH1", meanSM )
#    rrv2_gamma1 = ROOT.RooRealVar("rrv2_gamma1","rrv2_gamma1", gammaSM )
#    
#    rrv2_mH2 = ROOT.RooRealVar("rrv2_mH2","rrv2_mH2", meanSM )
#    rrv2_gamma2 = ROOT.RooRealVar("rrv2_gamma2","rrv2_gamma2", gammaSM*Cprime)
#
#    model1_pdf2 = ROOT.RooRelBWRunningWidth("model1_pdf2","model1_pdf2",rrv2_mass,rrv2_mH1,rrv2_gamma1);
#    model2_pdf2 = ROOT.RooRelBWRunningWidth("model2_pdf2","model2_pdf2",rrv2_mass,rrv2_mH2,rrv2_gamma2);
#
#    rrv2_mass.setVal(point);
#    weight = model2_pdf2.getVal()/model1_pdf2.getVal();    
#    
##    del model1_pdf2, model2_pdf2
##    del rrv2_mH1, rrv2_mH2, rrv2_gamma1, rrv2_gamma2;
#    
#    model1_pdf2.Delete(); model2_pdf2.Delete();
#    rrv2_mH1.Delete(); rrv2_mH2.Delete(); rrv2_gamma1.Delete(); rrv2_gamma2.Delete();
#    rrv2_mass.Delete();
#     
#    return weight;
    
    weight2 = MyRunningWidthBW(point,meanSM,gammaSM*Cprime)/MyRunningWidthBW(point,meanSM,gammaSM);
#    print "weight: ", weight,", weight2: ", weight2
    return weight2;

def IntfRescale(curIntfRw,cPrime,BRnew):

    curIoverS = curIntfRw - 1;
    newWeight = 1 + ((1-BRnew)*curIoverS)/cPrime;
    if newWeight < 0: newWeight = 0.01;
    ratio = newWeight/curIntfRw;
        
    return ratio;


#if __name__ == '__main__':
#
#
##    FitMassPoint( 900, 400, 1600 );
#
#    masses_in = [500, 550, 600, 700, 800, 900, 1000]
#    masses_min = [200, 200, 200, 250, 400, 400, 400]
#    masses_max = [1000, 1000, 1000, 1200, 1600, 1600, 1600]
#    gammaSM = [68.0,93.0,123.0,199.0,304.0,449.0,647.0]
#    
#    outputpars = [];
#    
#    for i in range(len(masses_in)):
#        outputpars.append( FitMassPoint( masses_in[i], masses_min[i], masses_max[i] ) ); 
#
#    for i in range(len(outputpars)):
#        print "mass: ", masses_in[i], ", gammaSM: ", gammaSM[i], ", gammaFit: ", outputpars[i][1], ", meanFit: ", outputpars[i][0] ;


