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
ROOT.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPalette(1);

ROOT.gSystem.Load("PDFs/RooRelBWHighMass_cxx.so")
ROOT.gSystem.Load("PDFs/RooRelBWRunningWidth_cxx.so")

############################################################
############################################################

def FitMassPoint(massin, massmin, massmax, nbins=50):

#    massin = 800;
#    massmin = 400;
#    massmax = 1400;
#    nbins = 50;
#    massin = 600;
#    massmin = 200;
#    massmax = 1000;
#    nbins = 40;
    
    inputFile = ROOT.TFile("/uscms_data/d2/andersj/Wjj/2012/data/Moriond2013/ReducedTrees/RD_mu_HWWMH"+str(massin)+"_CMSSW532_private.root");
    tree = inputFile.Get("WJet");
#    tree.Print();
    
    print "n entries: ", tree.GetEntries();
    
    ############################################################
    # RooFitting
    
    rrv_mass = ROOT.RooRealVar("rrv_mass","rrv_mass",massmin,massmax)
    rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 

    rrv_mH = ROOT.RooRealVar("rrv_mH","rrv_mH", massin, massmin, massmax )
    rrv_gamma = ROOT.RooRealVar("rrv_gamma","rrv_gamma",20.,massmax)
    rrv_mH2 = ROOT.RooRealVar("rrv_mH2","rrv_mH2", massin, massmin, massmax )
    rrv_gamma2 = ROOT.RooRealVar("rrv_gamma2","rrv_gamma2",20.,massmax)
#    rrv_gamma2 = ROOT.RooRealVar("rrv_gamma2","rrv_gamma2",350.,600.)    
    rrv_mH3 = ROOT.RooRealVar("rrv_mH3","rrv_mH3", massin, massmin, massmax )
    rrv_gamma3 = ROOT.RooRealVar("rrv_gamma3","rrv_gamma3",20.,massmax)
    rrv_mH4 = ROOT.RooRealVar("rrv_mH4","rrv_mH4", massin, massmin, massmax )
    rrv_gamma4 = ROOT.RooRealVar("rrv_gamma4","rrv_gamma4",20.,massmax)

    rds_raw = ROOT.RooDataSet("rds_raw","rds_raw",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_cps = ROOT.RooDataSet("rds_cps","rds_cps",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_cps_intf = ROOT.RooDataSet("rds_cps_intf","rds_cps_intf",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight));

    model_pdf = ROOT.RooRelBWHighMass("model_pdf","model_pdf",rrv_mass,rrv_mH,rrv_gamma);
    model2_pdf = ROOT.RooRelBWRunningWidth("model2_pdf","model2_pdf",rrv_mass,rrv_mH2,rrv_gamma2);

    model3_pdf = ROOT.RooRelBWRunningWidth("model3_pdf","model3_pdf",rrv_mass,rrv_mH3,rrv_gamma3);
    model4_pdf = ROOT.RooRelBWRunningWidth("model4_pdf","model4_pdf",rrv_mass,rrv_mH4,rrv_gamma4);

    for i in range(tree.GetEntries()):
        
        if i % 10000 == 0: print "i: ", i;
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

    rfresult = model_pdf.fitTo(rds_cps,RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
    rfresult2 = model2_pdf.fitTo(rds_cps,RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
    rrv_mH3.setVal(rrv_mH2.getVal());
    rrv_gamma3.setVal(0.5*rrv_gamma2.getVal());
    rrv_mH4.setVal(rrv_mH2.getVal());
    rrv_gamma4.setVal(0.2*rrv_gamma2.getVal());

    mplot = rrv_mass.frame(RooFit.Title("mass plot"));
    rds_raw.plotOn(mplot, RooFit.MarkerColor(kBlack), RooFit.LineColor(kBlack), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    rds_cps_intf.plotOn(mplot, RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    rds_cps.plotOn(mplot, RooFit.MarkerColor(kRed), RooFit.LineColor(kRed), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
#    model_pdf.plotOn(mplot, RooFit.LineColor(kRed));
    model2_pdf.plotOn(mplot, RooFit.LineColor(kRed), RooFit.LineStyle(2) );
    model3_pdf.plotOn(mplot, RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(2) );
    model4_pdf.plotOn(mplot, RooFit.LineColor(ROOT.kGreen+2), RooFit.LineStyle(3) );

    print "rds_cps_intf.sumEntries() = ", rds_cps_intf.sumEntries()
    print "model_pdf: mH = ", rrv_mH.getVal(), ", gamma = ", rrv_gamma.getVal();
    print "model2_pdf: mH = ", rrv_mH2.getVal(), ", gamma = ", rrv_gamma2.getVal();    

    dummy_h1 = ROOT.TH1F("dummy_h1","dummy_h1",1,0,1); 
    dummy_h1.SetMarkerColor( ROOT.kBlack );
    dummy_h2 = ROOT.TH1F("dummy_h2","dummy_h2",1,0,1); 
    dummy_h2.SetMarkerColor( ROOT.kRed );
    dummy_h3 = ROOT.TH1F("dummy_h3","dummy_h3",1,0,1); 
    dummy_h3.SetMarkerColor( ROOT.kBlue );
    dummy_h4 = ROOT.TH1F("dummy_h4","dummy_h4",1,0,1); 
    dummy_h4.SetLineColor( ROOT.kRed );
    dummy_h5 = ROOT.TH1F("dummy_h5","dummy_h5",1,0,1); 
    dummy_h5.SetLineColor( ROOT.kRed ); dummy_h5.SetLineStyle( 2 );
    dummy_h6 = ROOT.TH1F("dummy_h6","dummy_h6",1,0,1); 
    dummy_h6.SetLineColor( ROOT.kGreen+2 ); dummy_h6.SetLineStyle( 2 );
    dummy_h7 = ROOT.TH1F("dummy_h7","dummy_h7",1,0,1); 
    dummy_h7.SetLineColor( ROOT.kGreen+2 ); dummy_h6.SetLineStyle( 3 );

    L = TLegend(0.65,0.60,0.93,0.85);
    L.SetFillStyle(0);
    L.AddEntry(dummy_h1,"Powheg","p");
    L.AddEntry(dummy_h2,"w/CPS weight","p");
    L.AddEntry(dummy_h3,"w/CPS,Intf weight","p");
#    L.AddEntry(dummy_h4,"Fit, BW (Mario)","l");
    L.AddEntry(dummy_h5,"Fit, BW (running)","l");
    L.AddEntry(dummy_h6,"BW (running), width*0.5","l");
    L.AddEntry(dummy_h7,"Fit, BW (running), width*0.2","l");


    can2 = ROOT.TCanvas("can2","can2",800,800);
    mplot.Draw();
    L.Draw();
    ROOT.gPad.SetLogy();
    can2.SaveAs("massFits/mass_rf_"+str(massin)+".eps");
    can2.SaveAs("massFits/mass_rf_"+str(massin)+".png");

    outputpar = [];
    outputpar.append( rrv_mH2.getVal() );
    outputpar.append( rrv_gamma2.getVal() );
    outputpar.append( rrv_mH.getVal() );
    outputpar.append( rrv_gamma.getVal() );
    return outputpar

if __name__ == '__main__':


#    FitMassPoint( 900, 400, 1600 );

    masses_in = [500, 550, 600, 700, 800, 900, 1000]
    masses_min = [200, 200, 200, 250, 400, 400, 400]
    masses_max = [1000, 1000, 1000, 1200, 1600, 1600, 1600]
    gammaSM = [68.0,93.0,123.0,199.0,304.0,449.0,647.0]
    
    outputpars = [];
    
    for i in range(len(masses_in)):
        outputpars.append( FitMassPoint( masses_in[i], masses_min[i], masses_max[i] ) ); 

    for i in range(len(outputpars)):
        print "mass: ", masses_in[i], ", gammaSM: ", gammaSM[i], ", gammaFit: ", outputpars[i][1], ", meanFit: ", outputpars[i][0] ;

#    ############################################################
#    # simple histograms...
#    # histograms
#    if runSimple:
#        m_raw = ROOT.TH1F("m_raw","m_raw",50, 200, 1000);
#        m_cps = ROOT.TH1F("m_cps","m_cps",50, 200, 1000);    
#        m_cps_intf = ROOT.TH1F("m_cps_intf","m_cps_intf",50, 200, 1000);        
#        
#        for i in range(tree.GetEntries()):
#            if i % 10000 == 0: print "i: ", i;
#            
#            tree.GetEntry(i);
#            m_raw.Fill( getattr(tree,"W_H_mass_gen") );
#            m_cps.Fill( getattr(tree,"W_H_mass_gen"), getattr(tree,"complexpolewtggH600")/getattr(tree,"avecomplexpolewtggH600") );
#            m_cps_intf.Fill( getattr(tree,"W_H_mass_gen"), getattr(tree,"complexpolewtggH600")*getattr(tree,"interferencewtggH600")/getattr(tree,"avecomplexpolewtggH600") );
#        
#        
#        # draw
#        can1 = ROOT.TCanvas("can1","can1",800,800);
#        
#        m_raw.SetMaximum( 1.2*max(m_raw.GetMaximum(),m_cps.GetMaximum(),m_cps_intf.GetMaximum()) );
#        
#        m_raw.Draw();
#        m_cps.SetLineColor( ROOT.kRed );
#        m_cps.Draw("sames");
#        m_cps_intf.SetLineColor( ROOT.kBlue );
#        m_cps_intf.Draw("sames");
#        
#        can1.SaveAs("mass.eps");
#    ############################################################
