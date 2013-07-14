import os
import glob
import math

import ROOT
from ROOT import *
from ROOT import gROOT, gStyle, gSystem, TLatex
import subprocess

def makeCanvas(hmc, hmcsig, hd, name, legendOrientation, normalizeSignalToBkg, plotData):
    
    xmin = 0.;
    xmax = 0.;
    if legendOrientation == "left":
        xmin = 0.2;
        xmax = 0.4;
    elif legendOrientation == "right":
        xmin = 0.7;
        xmax = 0.9;
    elif legendOrientation == "center":
        xmin = 0.45;
        xmax = 0.65;
    else: 
        print "UNKNOWN ORIENTATION!! Default left";
        xmin = 0.2;
        xmax = 0.4;

    names = ["W+jets","WW","WZ","ZZ","ttbar"]
    
#    print hmc[0].GetTitle(), ", ",hmc[0].GetName(),", ",hmc[0].GetXaxis().GetName(),", ",hmc[0].GetXaxis().GetTitle();
    xtit = str(hmc[0].GetXaxis().GetTitle());
    ytit = str(hmc[0].GetYaxis().GetTitle());
    
    hstack = ROOT.THStack("hstack",hmc[0].GetTitle());
    for i in range(len(hmc)):
        hmc[i].SetFillColor( i+2 );
        hmc[i].SetFillStyle( 1001 );
        hstack.Add( hmc[i] );
        print "integral: ", hmc[i].Integral();
    
    if normalizeSignalToBkg: hmcsig.Scale( hstack.GetStack().Last().Integral()/hmcsig.Integral() );
    hmcsig.SetLineStyle( 2 );
    hmcsig.SetLineWidth( 3 );
    hmcsig.SetLineColor( 4 );

    maxmax = max(hstack.GetStack().Last().GetMaximum(),hmcsig.GetMaximum(),hd.GetMaximum());
    if not plotData: maxmax = max( hstack.GetStack().Last().GetMaximum(), hmcsig.GetMaximum() );

    leg = ROOT.TLegend(xmin,0.6,xmax,0.9)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)                
    for i in range(len(hmc)):
        leg.AddEntry( hmc[i], names[i], 'f' );
    leg.AddEntry( hmcsig, "MH600", 'l' );
    if plotData: leg.AddEntry( hd, "data", 'p' );                 
    
    c = ROOT.TCanvas("c","c",800,800);
    hstack.SetMaximum( maxmax * 1.2 );
    hstack.Draw("hist");
    hmcsig.Draw("histsames");
    if plotData: hd.Draw("pesames");

    hstack.GetXaxis().SetTitle( xtit );
    hstack.GetYaxis().SetTitle( ytit );
    leg.Draw();
#    ROOT.gPad.SetLogy();
    c.Update();
    
    c.SaveAs( name + ".eps");
    c.SaveAs( name + ".png");

##################################################################################################
##################################################################################################
##################################################################################################
class plotterClass:


    ### ------------------------------------------------
    def __init__(self, signalMC, bkgMCs, data):
        
        self.sigMCfile_ = ROOT.TFile( signalMC.getTrainingTreeName() );
        self.sigMClabel_ = signalMC.getSampleLabel();        
        self.bkgMCfile_ = [];
        self.bkgMClabel_ = [];        

        for i in range(len(bkgMCs)): 
            self.bkgMCfile_.append( ROOT.TFile(bkgMCs[i].getTrainingTreeName()) );
            self.bkgMClabel_.append( bkgMCs[i].getSampleLabel() );        
        
        self.datMCfile_ = ROOT.TFile(data.getTrainingTreeName());
        self.datMClabel_ = data.getSampleLabel();

        # make a list of the histograms you want to make and their parameters, bin,max,min,name    
        self.hvars = []; self.hvarsNBins = {}; self.hvarsMin = {}; self.hvarsMax = {}; self.hvarsXName = {}; self.hvarsLeg = {};
        self.hvars.append("v_pt"); 
        self.hvarsNBins[self.hvars[0]] = 25.; self.hvarsMin[self.hvars[0]] = 100.; self.hvarsMax[self.hvars[0]] = 600.; self.hvarsXName[self.hvars[0]] = "W pT"; self.hvarsLeg[self.hvars[0]] = "right";
        self.hvars.append("jet_tau2tau1"); 
        self.hvarsNBins[self.hvars[1]] = 30; self.hvarsMin[self.hvars[1]] = 0.; self.hvarsMax[self.hvars[1]] = 1.5; self.hvarsXName[self.hvars[1]] = "tau2/tau1"; self.hvarsLeg[self.hvars[1]] = "right";
        self.hvars.append("jet_pt_pr"); 
        self.hvarsNBins[self.hvars[2]] = 30; self.hvarsMin[self.hvars[2]] = 0.; self.hvarsMax[self.hvars[2]] = 600; self.hvarsXName[self.hvars[2]] = "pruned jet pT"; self.hvarsLeg[self.hvars[2]] = "right";
        self.hvars.append("jet_mass_pr"); 
        self.hvarsNBins[self.hvars[3]] = 30; self.hvarsMin[self.hvars[3]] = 0.; self.hvarsMax[self.hvars[3]] = 150; self.hvarsXName[self.hvars[3]] = "pruned mass"; self.hvarsLeg[self.hvars[3]] = "right";
        self.hvars.append("l_pt"); 
        self.hvarsNBins[self.hvars[4]] = 30; self.hvarsMin[self.hvars[4]] = 0.; self.hvarsMax[self.hvars[4]] = 600; self.hvarsXName[self.hvars[4]] = "lepton pT"; self.hvarsLeg[self.hvars[4]] = "right";
        self.hvars.append("l_eta"); 
        self.hvarsNBins[self.hvars[5]] = 25; self.hvarsMin[self.hvars[5]] = -2.5; self.hvarsMax[self.hvars[5]] = 2.5; self.hvarsXName[self.hvars[5]] = "lepton eta"; self.hvarsLeg[self.hvars[5]] = "right";
        self.hvars.append("mvaMET"); 
        self.hvarsNBins[self.hvars[6]] = 50; self.hvarsMin[self.hvars[6]] = 0.; self.hvarsMax[self.hvars[6]] = 500; self.hvarsXName[self.hvars[6]] = "mvaMET"; self.hvarsLeg[self.hvars[6]] = "right";
        self.hvars.append("njets"); 
        self.hvarsNBins[self.hvars[7]] = 5; self.hvarsMin[self.hvars[7]] = 0.; self.hvarsMax[self.hvars[7]] = 5; self.hvarsXName[self.hvars[7]] = "n extra jets"; self.hvarsLeg[self.hvars[7]] = "right";          
        self.hvars.append("nPV"); 
        self.hvarsNBins[self.hvars[8]] = 50; self.hvarsMin[self.hvars[8]] = 0.; self.hvarsMax[self.hvars[8]] = 50; self.hvarsXName[self.hvars[8]] = "n PV"; self.hvarsLeg[self.hvars[8]] = "right";   
        self.hvars.append("jet_grsens_ft"); 
        self.hvarsNBins[self.hvars[9]] = 30; self.hvarsMin[self.hvars[9]] = 0.; self.hvarsMax[self.hvars[9]] = 1.; self.hvarsXName[self.hvars[9]] = "filtered gr. sens."; self.hvarsLeg[self.hvars[9]] = "left";           
        self.hvars.append("jet_grsens_tr"); 
        self.hvarsNBins[self.hvars[10]] = 30; self.hvarsMin[self.hvars[10]] = 0.; self.hvarsMax[self.hvars[10]] = 1.; self.hvarsXName[self.hvars[10]] = "trimmed gr. sens."; self.hvarsLeg[self.hvars[10]] = "left";            
        self.hvars.append("jet_massdrop_pr"); 
        self.hvarsNBins[self.hvars[11]] = 30; self.hvarsMin[self.hvars[11]] = 0.; self.hvarsMax[self.hvars[11]] = 1.; self.hvarsXName[self.hvars[11]] = "mass drop"; self.hvarsLeg[self.hvars[11]] = "center";            
        self.hvars.append("jet_qjetvol"); 
        self.hvarsNBins[self.hvars[12]] = 30; self.hvarsMin[self.hvars[12]] = 0.; self.hvarsMax[self.hvars[12]] = 0.3; self.hvarsXName[self.hvars[12]] = "QJet volatility"; self.hvarsLeg[self.hvars[12]] = "right";            
        self.hvars.append("jet_jetconstituents"); 
        self.hvarsNBins[self.hvars[13]] = 50; self.hvarsMin[self.hvars[13]] = 0.; self.hvarsMax[self.hvars[13]] = 150; self.hvarsXName[self.hvars[13]] = "n constituents"; self.hvarsLeg[self.hvars[13]] = "right";            
        self.hvars.append("jet_rcore4"); 
        self.hvarsNBins[self.hvars[14]] = 30; self.hvarsMin[self.hvars[14]] = 0.; self.hvarsMax[self.hvars[14]] = 1.; self.hvarsXName[self.hvars[14]] = "R-core 0.4"; self.hvarsLeg[self.hvars[14]] = "right";            
        self.hvars.append("jet_rcore5"); 
        self.hvarsNBins[self.hvars[15]] = 30; self.hvarsMin[self.hvars[15]] = 0.; self.hvarsMax[self.hvars[15]] = 1.; self.hvarsXName[self.hvars[15]] = "R-core 0.5"; self.hvarsLeg[self.hvars[15]] = "right";            
        self.hvars.append("jet_rcore6"); 
        self.hvarsNBins[self.hvars[16]] = 30; self.hvarsMin[self.hvars[16]] = 0.; self.hvarsMax[self.hvars[16]] = 1.; self.hvarsXName[self.hvars[16]] = "R-core 0.6"; self.hvarsLeg[self.hvars[16]] = "right";            
        self.hvars.append("jet_rcore7"); 
        self.hvarsNBins[self.hvars[17]] = 30; self.hvarsMin[self.hvars[17]] = 0.; self.hvarsMax[self.hvars[17]] = 1.; self.hvarsXName[self.hvars[17]] = "R-core 0.7"; self.hvarsLeg[self.hvars[17]] = "left";            
        self.hvars.append("jet_planarflow04"); 
        self.hvarsNBins[self.hvars[18]] = 30; self.hvarsMin[self.hvars[18]] = 0.; self.hvarsMax[self.hvars[18]] = 1.; self.hvarsXName[self.hvars[18]] = "Planar Flow 0.4"; self.hvarsLeg[self.hvars[18]] = "left";            
        self.hvars.append("jet_planarflow05"); 
        self.hvarsNBins[self.hvars[19]] = 30; self.hvarsMin[self.hvars[19]] = 0.; self.hvarsMax[self.hvars[19]] = 1.; self.hvarsXName[self.hvars[19]] = "Planar Flow 0.5"; self.hvarsLeg[self.hvars[19]] = "left";            
        self.hvars.append("jet_planarflow06"); 
        self.hvarsNBins[self.hvars[20]] = 30; self.hvarsMin[self.hvars[20]] = 0.; self.hvarsMax[self.hvars[20]] = 1.; self.hvarsXName[self.hvars[20]] = "Planar Flow 0.6"; self.hvarsLeg[self.hvars[20]] = "left";            
        self.hvars.append("jet_planarflow07"); 
        self.hvarsNBins[self.hvars[21]] = 30; self.hvarsMin[self.hvars[21]] = 0.; self.hvarsMax[self.hvars[21]] = 1.; self.hvarsXName[self.hvars[21]] = "Planar Flow 0.7"; self.hvarsLeg[self.hvars[21]] = "left";            
        self.hvars.append("nbjets"); 
        self.hvarsNBins[self.hvars[22]] = 3; self.hvarsMin[self.hvars[22]] = 0.; self.hvarsMax[self.hvars[22]] = 3; self.hvarsXName[self.hvars[22]] = "n b-jets"; self.hvarsLeg[self.hvars[22]] = "right";            
        self.hvars.append("jet_pt1frac"); 
        self.hvarsNBins[self.hvars[23]] = 30; self.hvarsMin[self.hvars[23]] = 0.; self.hvarsMax[self.hvars[23]] = 1.; self.hvarsXName[self.hvars[23]] = "pt1/pt"; self.hvarsLeg[self.hvars[23]] = "left";            
        self.hvars.append("jet_pt2frac"); 
        self.hvarsNBins[self.hvars[24]] = 30; self.hvarsMin[self.hvars[24]] = 0.; self.hvarsMax[self.hvars[24]] = 1.; self.hvarsXName[self.hvars[24]] = "pt2/pt"; self.hvarsLeg[self.hvars[24]] = "right";           
        self.hvars.append("jet_sjdr"); 
        self.hvarsNBins[self.hvars[25]] = 30; self.hvarsMin[self.hvars[25]] = 0.; self.hvarsMax[self.hvars[25]] = 0.8; self.hvarsXName[self.hvars[25]] = "DR subjets"; self.hvarsLeg[self.hvars[25]] = "center";            
        self.hvars.append("deltaR_lca8jet"); 
        self.hvarsNBins[self.hvars[26]] = 30; self.hvarsMin[self.hvars[26]] = 0.; self.hvarsMax[self.hvars[26]] = 3.0; self.hvarsXName[self.hvars[26]] = "DR l,j"; self.hvarsLeg[self.hvars[26]] = "left";            
        self.hvars.append("deltaphi_METca8jet"); 
        self.hvarsNBins[self.hvars[27]] = 30; self.hvarsMin[self.hvars[27]] = 0.; self.hvarsMax[self.hvars[27]] = 3.15; self.hvarsXName[self.hvars[27]] = "DPhi MET,j"; self.hvarsLeg[self.hvars[27]] = "left";            
        self.hvars.append("deltaphi_Vca8jet"); 
        self.hvarsNBins[self.hvars[28]] = 30; self.hvarsMin[self.hvars[28]] = 0.; self.hvarsMax[self.hvars[28]] = 3.15; self.hvarsXName[self.hvars[28]] = "DPhi V,j"; self.hvarsLeg[self.hvars[28]] = "left";            
            
        self.cBTagMax = 10;
        self.cBTagMin = 0;
        self.cMassMin = 0;    
        self.cMassMax = 1000;    
        self.cPt1Frac = 1.0;
        self.cDRsubjets = 0.0;
        self.cPtMin = 0.0;
        self.cPtMax = 1000.0;
        self.cNJetsMax = 10.;    

    ## ------------------------
    ## ------------------------
    def SetDefaultCuts(self):

        self.cBTagMax = 10;
        self.cBTagMin = 0;
        self.cMassMin = 0;    
        self.cMassMax = 1000;    
        self.cPt1Frac = 1.0;
        self.cDRsubjets = 0.0;
        self.cPtMin = 0.0;
        self.cPtMax = 1000.0;
        self.cNJetsMax = 10.;

    ## ------------------------
    ## ------------------------
    def makeControlPlots(self, odir, tag):
    
        normalizeSignalToBkg = False;
        plotData = True;
        if tag == "ttbar":
            self.cBTagMax = 10;
            self.cBTagMin = 1;
            self.cDRsubjets = 0.3;
            self.cPtMin = 220.0;
            normalizeSignalToBkg = True;   
        elif tag == "signalregion":
            self.cBTagMax = 1;
            self.cBTagMin = 0;
            self.cMassMin = 60;    
            self.cMassMax = 100;    
            self.cDRsubjets = 0.;
            self.cPtMin = 220.0;            
            self.cNJetsMax = 1;                        
            normalizeSignalToBkg = True;   
            plotData = False;
        else: 
            self.SetDefaultCuts();

        ## make bkg histos
        bkgHistos = []
        for i in range(len(self.bkgMCfile_)): 
            curhisto = self.makeHistograms( self.bkgMCfile_[i], self.bkgMClabel_[i] );
            bkgHistos.append( curhisto );

        ## make dat histos
        datHistos = self.makeHistograms( self.datMCfile_, self.datMClabel_);

        ## make sig histos
        ## set default cuts here for ttbar to plot signal without b-tag cut
        if tag == "ttbar": 
            self.SetDefaultCuts();
            self.cDRsubjets = 0.3;
            self.cPtMin = 220.0;

        sigHistos = self.makeHistograms( self.sigMCfile_, self.sigMClabel_ );

        self.SetDefaultCuts();

        for j in range(len(self.hvars)):
            mcbkghistos = [];
            for k in range(len(self.bkgMCfile_)): mcbkghistos.append( bkgHistos[k][j] );
            makeCanvas(mcbkghistos, sigHistos[j], datHistos[j], odir+"/"+self.hvars[j], self.hvarsLeg[self.hvars[j]], normalizeSignalToBkg, plotData);            
            
    ## ------------------------
    ## ------------------------
    def makeHistograms( self, file, label ):
        
        tree = file.Get("otree");
        # define histograms
        histograms = [];
        for i in range(len(self.hvars)):
            name = "h_"+self.hvars[i]+"_"+label
            title = ";"+self.hvarsXName[self.hvars[i]]+";a.u.";
            bins = self.hvarsNBins[self.hvars[i]];
            min = self.hvarsMin[self.hvars[i]];
            max = self.hvarsMax[self.hvars[i]];
#            print name, ",",title,",",bins,",",min,",",max
            dummyhist = ROOT.TH1F( name, title, int(bins), min, max );
            histograms.append( dummyhist );
        # fill histograms
        print "n: ", tree.GetEntries();
        print "cuts: ", self.cBTagMax, ",",self.cBTagMin,",",self.cPt1Frac,",",self.cMassMin,",",self.cMassMax
        for i in range(tree.GetEntries()):
            if i % 10000 == 0: print "i: ", i
            
            tree.GetEntry(i);
            
            if tree.nbjets >= self.cBTagMin and tree.nbjets < self.cBTagMax and tree.jet_mass_pr > self.cMassMin and tree.jet_mass_pr < self.cMassMax and tree.jet_pt1frac < self.cPt1Frac and tree.jet_sjdr > self.cDRsubjets and tree.jet_pt_pr > self.cPtMin and tree.jet_pt_pr < self.cPtMax and tree.njets <         self.cNJetsMax:
                for j in range(len(self.hvars)):
                    histograms[j].Fill( getattr( tree, self.hvars[j] ), getattr( tree,"totalEventWeight" ) );

        return histograms

