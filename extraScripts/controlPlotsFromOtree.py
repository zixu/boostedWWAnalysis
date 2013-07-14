#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
from array import array
from ROOT import gROOT, gStyle, gSystem, TLatex
from ROOT import RooArgSet, RooFit

import subprocess
from subprocess import Popen

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C")
from ROOT import setTDRStyle
ROOT.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPalette(1);

from optparse import OptionParser

############################################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--efficiency', action='store_true', dest='efficiency', default=False, help='no X11 windows')
parser.add_option('--sigCompare', action='store_true', dest='sigCompare', default=False, help='no X11 windows')

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('-m','--mass', action="store", type="int", dest="mass", default=600, help='mass to check')

(options, args) = parser.parse_args()
############################################################

############################################################
# G l o b a l   I n p u t s

inputDir = "../trainingtrees_"+options.channel+"/";

# ---------------
# observables
# ---------------

obsAll = [];
#obsAll.append( ["mass_lvj","m_{l#nuJ}",20,400,2400] );
obsAll.append( ["jet_mass_pr","m_{J}",75,0,150] );
#obsAll.append( ["ungroomed_jet_pt","pT_{J}",20,100,600] );
#obsAll.append( ["jet_tau2tau1","#tau_{2}/#tau_{1}",30,0.,1.5] );
#obsAll.append( ["jet_massdrop_pr","#mu = m_{1}/m_{J}",20,0.,1.0] );

#obsAll.append( ["v_pt","V pT",50,0,500] );
#obsAll.append( ["l_pt","lepton pT",40,0,400] );
#obsAll.append( ["pfMET","pfMET",40,0,400] );
#obsAll.append( ["deltaphi_METca8jet","deltaphi_METca8jet",16,0,3.2] );
#obsAll.append( ["deltaphi_Vca8jet","deltaphi_Vca8jet",16,0,3.2] );
#
#obsAll.append( ["ttb_ca8_mass_pr","m_{J}",30,0,150] );
#obsAll.append( ["ttb_ht","HT - pT_{W} - pT_{ca8}",35,0,700] );
#obsAll.append( ["ttb_ca8_ungroomed_pt","pT_{J}",20,200,600] );
#obsAll.append( ["ttb_ca8_tau2tau1","#tau_{2}/#tau_{1}",30,0.,1.1] );
#obsAll.append( ["ttb_ca8_mu","#mu = m_{1}/m_{J}",20,0.,1.0] );
#obsAll.append( ["ttb_mlvj","m_{l#nuJ}",40,400,1400] );

obs=[]
obs_xname=[]
obs_lo=[]
obs_hi=[]
obs_bin=[]

for i in range(len(obsAll)):
    obs.append(obsAll[i][0]);
    obs_xname.append(obsAll[i][1]);
    obs_lo.append(obsAll[i][3]);
    obs_hi.append(obsAll[i][4]);
    obs_bin.append(obsAll[i][2]);   

print obs


# ---------------
# data file names
# ---------------

dataFileName = "ofile_data.root";
dataTag = "dat";
#mcFileNames = ["ofile_WJets_Herwig.root","ofile_TTbar.root","ofile_VV.root","ofile_STop.root"];
#mcTag = ["WJets","ttbar","VV","sTop"];
#mcColors = [ROOT.kGreen, ROOT.kRed, ROOT.kBlue, ROOT.kCyan];
mcFileNames = ["ofile_WJets_Herwig.root","ofile_VV.root","ofile_STop.root","ofile_TTbar.root"];
mcTag = ["WJets","VV","sTop","ttbar"];
mcColors = [ROOT.kGreen, ROOT.kBlue, ROOT.kCyan, ROOT.kRed];
#mcFileNames = ["ofile_VV.root","ofile_STop.root","ofile_TTbar.root","ofile_WJets_Herwig.root"];
#mcTag = ["VV","sTop","ttbar","WJets"];
#mcColors = [ROOT.kBlue, ROOT.kCyan, ROOT.kRed, ROOT.kGreen];

#signalFileName = "ofile_ggH600.root";
#sigTag = "ggH600";
signalFileName = "ofile_ggH1000.root";
sigTag = "200-300GeV";
signalFileName2 = "ofile_ggH1000.root";
sigTag2 = "450-550GeV";


# ---------------
# cuts
# ---------------


## GENERIC CUTS
cutVpt = 200;
cutLpt = 30;
cutMET = 50;
if options.channel == "el":
    cutLpt = 35;
    cutMET = 70;
cut_deltaphi_METca8jet = 1.57;

cutIssignal = 1;

cutNbjetsSSVHEveto = 0;
cut_t2t1 = 100.75;
cut_mu = 100.25;
cut_mJlo = 0;
cut_mJhi = 150;
cut_mlvJlo = 400;
cut_mlvJhi = 2400;
cut_jptLo = 200;
cut_jptHi = 2000;
##cut_deltaR!
##cut_deltaPhi!


if cutIssignal == 0:
    cutNbjetsSSVHEveto = 100;
    cut_t2t1 = 100;
    cut_mu = 100;
    cut_mJlo = -100;
    cut_mJhi = 10000;
    cut_mlvJlo = -100;
    cut_mlvJhi = 10000;
    cut_jptLo = -100;
    cut_jptHi = 1e9;


cutIsttbar = 0; # looser region
cutIsttbar_v2 = 0; 

ttb_massCutLo = 40;
ttb_massCutHi = 130;    
ttb_t2t1Cut = 100.53;
#ttb_massCutLo = 70;
#ttb_massCutHi = 100;    
#ttb_t2t1Cut = 0.53;
ttb_muCut = 100.;
ttb_ptCut = 200.;
ttb_ptCutMax = 1000.;
if cutIsttbar == 0 and cutIsttbar_v2 == 0:
    ttb_massCutLo = -100;
    ttb_massCutHi = 1000;    
    ttb_t2t1Cut = 100.;
    ttb_muCut = 100.;
    ttb_ptCut = 100.;        
    ttb_ptCutMax = 10000.;


# ---------------
# cuts
# ---------------

g_normalizeSignalToBkg = False;
g_plotData = True;

############################################################

def DrawStackPlot(hd, hmc, hmcsig, legendOrientation="right", normalizeSignalToBkg=True, plotData=True):
    
    binLoVal = 65.5;
    binHiVal = 104.5
    
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

    xtit = str(hmc[0].GetXaxis().GetTitle());
    ytit = str(hmc[0].GetYaxis().GetTitle());

    print "n_data = ",hd.Integral(),";"
    print "a_n_data = ",hd.Integral( hd.GetXaxis().FindBin( binLoVal ), hd.GetXaxis().FindBin( binHiVal ) )

    totbkgInt = 0;
    totbkgIntWindow = 0;
    ittbar = -1;
    isTop = -1;
    iVV = -1;
    hstack = ROOT.THStack("hstack",hmc[0].GetTitle());
    for i in range(len(hmc)):
        hmc[i].SetFillColor( mcColors[i] );
        hmc[i].SetFillStyle( 1001 );
        hstack.Add( hmc[i] );
        totbkgInt += hmc[i].Integral()
        totbkgIntWindow += hmc[i].Integral( hmc[i].GetXaxis().FindBin( binLoVal ), hmc[i].GetXaxis().FindBin( binHiVal ) );
    
        if hmc[i].GetName().find("ttbar") > 0: ittbar = i;
        if hmc[i].GetName().find("sTop") > 0: isTop = i;
        if hmc[i].GetName().find("VV") > 0: iVV = i;        
        print "n_"+mcTag[i]+" = ", hmc[i].Integral(),";";
        print "a_n_"+mcTag[i]+" = ", hmc[i].Integral( hmc[i].GetXaxis().FindBin( binLoVal ), hmc[i].GetXaxis().FindBin( binHiVal ) );

    print "bkg integral: ",totbkgInt
#    print "bkg integral window: ",totbkgIntWindow

    sigSF = 1;
    if normalizeSignalToBkg and hmcsig.Integral() > 0: hmcsig.Scale( hstack.GetStack().Last().Integral()/hmcsig.Integral() );
    hmcsig.SetLineStyle( 2 );
    hmcsig.SetLineWidth( 3 );
    hmcsig.SetLineColor( 4 );
    if not normalizeSignalToBkg: 
        sigSF = 5; 
        hmcsig.Scale( sigSF );
    print "signal integral: ", hmcsig.Integral();
    print "signal integral window: ", hmcsig.Integral( hmcsig.GetXaxis().FindBin( binLoVal ), hmcsig.GetXaxis().FindBin( binHiVal ) );

    maxmax = max(hstack.GetStack().Last().GetMaximum(),hd.GetMaximum());
#    maxmax = max(hstack.GetStack().Last().GetMaximum(),hmcsig.GetMaximum(),hd.GetMaximum());
    if not plotData: maxmax = max( hstack.GetStack().Last().GetMaximum(), hmcsig.GetMaximum() );

    #####################################Error_Band###########################
    lastStack = hstack.GetStack().Last();
    xvalue = array('d')
    yvalue = array('d')
    xlefterror = array('d')
    xrighterror = array('d')
    ylowerror = array('d')
    yhigherror = array('d')
    
    xlist = []
    ylist = []
    xleftlist = []
    xrightlist = []
    ylowlist = []
    yhighlist = []
    
    for ibin in range(1,lastStack.GetNbinsX() + 1):
        xlist.append(lastStack.GetBinCenter(ibin))
        ylist.append(lastStack.GetBinContent(ibin))
        xleftlist.append(0.5 * lastStack.GetBinWidth(ibin))
        xrightlist.append(0.5 * lastStack.GetBinWidth(ibin))
        ttbarerror = hmc[ittbar].GetBinContent(ibin) * 0.07
        sToperror = hmc[isTop].GetBinContent(ibin) * 0.30
        VVerror = hmc[iVV].GetBinContent(ibin) * 0.30        
        #wjeterror = wjet1D.GetBinContent(ibin) * 0.3
        lumierror = lastStack.GetBinContent(ibin) * 0.044
        statisticerror = lastStack.GetBinError(ibin)
        allerror = math.sqrt(math.pow(ttbarerror,2) + math.pow(sToperror,2) + math.pow(VVerror,2) + math.pow(lumierror,2) + math.pow(statisticerror,2))
        #allerror = statisticerror             
        ylowlist.append(allerror)
        yhighlist.append(allerror)
    
    xvalue.fromlist(xlist)
    yvalue.fromlist(ylist)
    xlefterror.fromlist(xleftlist)
    xrighterror.fromlist(xrightlist)
    ylowerror.fromlist(ylowlist)
    yhigherror.fromlist(yhighlist)
    
    mc1Derror = ROOT.TGraphAsymmErrors(lastStack.GetNbinsX(),xvalue,yvalue,xlefterror,xrighterror,ylowerror,yhigherror)
    mc1Derror.SetName("MC Uncerntainty")
    mc1Derror.SetFillColor(920+3)
    mc1Derror.SetFillStyle(3008)
    #####################################Error_Band###########################

    leg = ROOT.TLegend(xmin,0.6,xmax,0.85)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)                
    for i in range(len(hmc)):
        leg.AddEntry( hmc[i], mcTag[i], 'f' );
    leg.AddEntry( hmcsig, "MH600x"+str(sigSF), 'l' );
    if plotData: leg.AddEntry( hd, "data", 'p' );   
    leg.AddEntry( mc1Derror,"stat #oplus sys",'f');

    c = ROOT.TCanvas("c","c",800,800);
    hstack.SetMaximum( maxmax * 1.3 );
    hstack.Draw("hist");
    hmcsig.Draw("histsames");
    if plotData: hd.Draw("pesames");
    mc1Derror.Draw("2 same")

    hstack.GetXaxis().SetTitle( xtit );
    hstack.GetYaxis().SetTitle( ytit );
    leg.Draw();
#    ROOT.gPad.SetLogy();
    c.Update();
    name = "controls/"+options.channel+"_"+hd.GetName();
    c.SaveAs( name + ".eps");
    c.SaveAs( name + ".png");

# ==================

def DrawStackPlot_sig(hd, hmcsig, legendOrientation="right", normalizeSignalToBkg=True, fit=False):
    
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
    
    xtit = str(hd.GetXaxis().GetTitle());
    ytit = str(hd.GetYaxis().GetTitle());
    
    sigSF = 1;
    if normalizeSignalToBkg and hd.Integral() > 0: hmcsig.Scale( hd.Integral()/hmcsig.Integral() );
    hmcsig.SetLineStyle( 2 );
    hmcsig.SetLineWidth( 3 );
    hmcsig.SetLineColor( 4 );
    
    maxmax = max(hmcsig.GetMaximum(),hd.GetMaximum());
    
    leg = ROOT.TLegend(xmin,0.6,xmax,0.85)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)                
    leg.AddEntry( hd, sigTag, 'l' );                 
    leg.AddEntry( hmcsig, sigTag2, 'l' );

    fd = ROOT.TF1("fd","gaus",70,90);
    fmcsig = ROOT.TF1("fmcsig","gaus",70,90);
    fd.SetLineColor(ROOT.kBlack);
    fmcsig.SetLineColor(ROOT.kBlue);

    c = ROOT.TCanvas("c","c",800,800);
    hd.SetMaximum(maxmax*1.2);
    hd.Draw("hist");
    hmcsig.Draw("histsames");
#    hd.Fit(fd,"R","sames");
#    hmcsig.Fit(fmcsig,"R","sames");
    leg.Draw();
    #    ROOT.gPad.SetLogy();
    c.Update();
    name = "controls/sigCompare_"+options.channel+"_"+hd.GetName();
    c.SaveAs( name + ".eps");
    c.SaveAs( name + ".png");


###############################################################
###############################################################
###############################################################

if __name__ == '__main__':

    
    dF = ROOT.TFile(inputDir+dataFileName);
    mF = [];
    for i in range(len(mcFileNames)): mF.append( ROOT.TFile( inputDir+mcFileNames[i] ) );
    sF = ROOT.TFile(inputDir+signalFileName);
    sF2 = ROOT.TFile(inputDir+signalFileName2);

    
    # Roofit stuff for fitting later
    workspace = ROOT.RooWorkspace("workspace","workspace");
    rrv_mj = ROOT.RooRealVar("rrv_mj","rrv_mj",40,130,"GeV/c^{2}");
    rrv_weight = ROOT.RooRealVar("rrv_weight","rrv_weight",0.,100.);    
    rds_data = ROOT.RooDataSet("rds_data","rds_data",ROOT.RooArgSet(rrv_mj,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_mc_VV = ROOT.RooDataSet("rds_mc_VV","rds_mc_VV",ROOT.RooArgSet(rrv_mj,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_mc_sTop = ROOT.RooDataSet("rds_mc_sTop","rds_mc_sTop",ROOT.RooArgSet(rrv_mj,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_mc_ttbar = ROOT.RooDataSet("rds_mc_ttbar","rds_mc_ttbar",ROOT.RooArgSet(rrv_mj,rrv_weight),RooFit.WeightVar(rrv_weight));
    rds_mc_WJets = ROOT.RooDataSet("rds_mc_WJets","rds_mc_WJets",ROOT.RooArgSet(rrv_mj,rrv_weight),RooFit.WeightVar(rrv_weight));    
    rds_mc_all = ROOT.RooDataSet("rds_mc_all","rds_mc_all",ROOT.RooArgSet(rrv_mj,rrv_weight),RooFit.WeightVar(rrv_weight));    

    ###############################################################
    ###############################################################
    if not options.efficiency and not options.sigCompare:
                        
        # 2D array of histograms, observables and MC samples
        dataHistArray = [];
        mcHistArray = [];
        sigHistArray = [];
        
        # ---------------
        # data loop
        # ---------------
        dT = dF.Get("otree");
        print "dat file: ", dataFileName;
        for i in range(len(obs)):
            dataHistArray.append( ROOT.TH1F("h_"+obs[i]+"_"+dataTag,";%s;count"%(obs_xname[i]),obs_bin[i],obs_lo[i],obs_hi[i]) );
        
        for i in range(dT.GetEntries()):
            dT.GetEntry(i);
            # set cuts
            isTTBar_v2 = 0;
#            oppo1same1m = dT.ttb_nak5_same_csvm > 0 and dT.ttb_nak5_oppoveto_csvl > 0;
#            oppo2same0m = dT.ttb_nak5_same_csvl > 0 and dT.ttb_nak5_oppoveto_csvm > 0;                                        
            oppo1same1m = dT.ttb_nak5_same_csvm > 0 or dT.ttb_nak5_oppoveto_csvm > 0;
            oppo2same0m = dT.ttb_nak5_same_csvl == 0 and dT.ttb_nak5_oppo_csvl > 1;
            oppo2same0m = False;
            if oppo1same1m or oppo2same0m:
                isTTBar_v2 = 1

            if getattr(dT, "issignal") >= cutIssignal and getattr(dT, "isttbar") >= cutIsttbar and dT.v_pt > cutVpt and dT.l_pt > cutLpt and dT.pfMET > cutMET and dT.ttb_ca8_mass_pr > ttb_massCutLo and dT.ttb_ca8_mass_pr < ttb_massCutHi and dT.ttb_ca8_tau2tau1 < ttb_t2t1Cut  and dT.ttb_ca8_mu < ttb_muCut and dT.jet_tau2tau1 < cut_t2t1  and dT.jet_massdrop_pr < cut_mu and isTTBar_v2 >= cutIsttbar_v2 and dT.nbjetsSSVHE <= cutNbjetsSSVHEveto and dT.jet_mass_pr > cut_mJlo and dT.jet_mass_pr < cut_mJhi and dT.mass_lvj > cut_mlvJlo and dT.mass_lvj < cut_mlvJhi and dT.ttb_ca8_ungroomed_pt > ttb_ptCut and dT.ttb_ca8_ungroomed_pt< ttb_ptCutMax and dT.deltaphi_METca8jet > cut_deltaphi_METca8jet:

#            if getattr(dT, "issignal") >= cutIssignal and getattr(dT, "isttbar") >= cutIsttbar and isTTBar_v2 >= cutIsttbar_v2 and dT.l_pt > cutLpt and dT.pfMET > cutMET:

                for j in range(len(obs)):
                    val = 0
                    if obs[j] == "ttb_ht": val = getattr( dT, obs[j] ) - dT.v_pt - dT.ttb_ca8_ungroomed_pt
                    else: val = getattr( dT, obs[j] );
                    dataHistArray[j].Fill( val, getattr(dT, "totalEventWeight") );

                rrv_mj.setVal( dT.ttb_ca8_mass_pr );
                rds_data.add( ROOT.RooArgSet( rrv_mj ), 1. );

        # ---------------
        # mc loop
        # ---------------
        for a in range(len(mcFileNames)):
            print "mc file: ", mcFileNames[a];
            
            mT = mF[a].Get("otree");
            curMCHists = [];
            for i in range(len(obs)):
                curMCHists.append( ROOT.TH1F("h_"+obs[i]+"_"+mcTag[a],";%s;count"%(obs_xname[i]),obs_bin[i],obs_lo[i],obs_hi[i]) );

            for i in range(mT.GetEntries()):
                mT.GetEntry(i);
                # set cuts
                isTTBar_v2 = 0;
#                oppo1same1m = mT.ttb_nak5_same_csvm > 0 and mT.ttb_nak5_oppoveto_csvl > 0;
#                oppo2same0m = mT.ttb_nak5_same_csvl > 0 and mT.ttb_nak5_oppoveto_csvm > 0;  
                oppo1same1m = mT.ttb_nak5_same_csvm > 0 or mT.ttb_nak5_oppoveto_csvm > 0;                
                oppo2same0m = mT.ttb_nak5_same_csvl == 0 and mT.ttb_nak5_oppo_csvl > 1;
                oppo2same0m = False;                
                if oppo1same1m or oppo2same0m:
                    isTTBar_v2 = 1
    #                print isTTBar_v2            
                if getattr(mT, "issignal") >= cutIssignal and getattr(mT, "isttbar") >= cutIsttbar and mT.v_pt > cutVpt and mT.l_pt > cutLpt and mT.pfMET > cutMET and mT.ttb_ca8_mass_pr > ttb_massCutLo and mT.ttb_ca8_mass_pr < ttb_massCutHi and mT.ttb_ca8_tau2tau1 < ttb_t2t1Cut and mT.ttb_ca8_mu < ttb_muCut and mT.jet_tau2tau1 < cut_t2t1  and mT.jet_massdrop_pr < cut_mu and isTTBar_v2 >= cutIsttbar_v2 and mT.nbjetsSSVHE <= cutNbjetsSSVHEveto and mT.jet_mass_pr > cut_mJlo and mT.jet_mass_pr < cut_mJhi and mT.mass_lvj > cut_mlvJlo and mT.mass_lvj < cut_mlvJhi and mT.ttb_ca8_ungroomed_pt > ttb_ptCut and mT.ttb_ca8_ungroomed_pt< ttb_ptCutMax and mT.deltaphi_METca8jet > cut_deltaphi_METca8jet:
                    for j in range(len(obs)):
                        val = 0
                        if obs[j] == "ttb_ht": val = getattr( mT, obs[j] ) - mT.v_pt - mT.ttb_ca8_ungroomed_pt
                        else: val = getattr( mT, obs[j] );
                        curMCHists[j].Fill( val, getattr(mT, "totalEventWeight") );

                    curweight_mc = getattr(mT, "totalEventWeight");
                    rrv_mj.setVal( mT.ttb_ca8_mass_pr );
                    rds_mc_all.add( ROOT.RooArgSet( rrv_mj ), curweight_mc );
                    if mcTag[a] == "WJets": rds_mc_WJets.add( ROOT.RooArgSet( rrv_mj ), curweight_mc );
                    if mcTag[a] == "sTop": rds_mc_sTop.add( ROOT.RooArgSet( rrv_mj ), curweight_mc );
                    if mcTag[a] == "ttbar": rds_mc_ttbar.add( ROOT.RooArgSet( rrv_mj ), curweight_mc );
                    if mcTag[a] == "VV": rds_mc_VV.add( ROOT.RooArgSet( rrv_mj ), curweight_mc );                    

            mcHistArray.append( curMCHists );
                
        getattr( workspace, "import" )(rrv_mj);
        getattr( workspace, "import" )(rds_data);
        getattr( workspace, "import" )(rds_mc_VV);
        getattr( workspace, "import" )(rds_mc_sTop);
        getattr( workspace, "import" )(rds_mc_ttbar);
        getattr( workspace, "import" )(rds_mc_WJets);
        getattr( workspace, "import" )(rds_mc_all);
        workspace.writeToFile("workspace_rds_"+options.channel+".root");                

        # ---------------
        # sig loop
        # ---------------
        sT = sF.Get("otree");
        print "sig file: ", signalFileName;
        for i in range(len(obs)):
            sigHistArray.append( ROOT.TH1F("h_"+obs[i]+"_"+sigTag,";%s;count"%(obs_xname[i]),obs_bin[i],obs_lo[i],obs_hi[i]) );

        for i in range(sT.GetEntries()):
            sT.GetEntry(i);
            # set cuts
            isTTBar_v2 = 0;
#            oppo1same1m = sT.ttb_nak5_same_csvm > 0 and sT.ttb_nak5_oppoveto_csvl > 0;
#            oppo2same0m = sT.ttb_nak5_same_csvl > 0 and sT.ttb_nak5_oppoveto_csvm > 0;                            
            oppo1same1m = sT.ttb_nak5_same_csvm > 0 or sT.ttb_nak5_oppoveto_csvm > 0;
            oppo2same0m = sT.ttb_nak5_same_csvl == 0 and sT.ttb_nak5_oppo_csvl > 1;
            oppo2same0m = False;            
            if oppo1same1m or oppo2same0m:
                isTTBar_v2 = 1
            
            if getattr(sT, "issignal") >= cutIssignal and getattr(sT, "isttbar") >= cutIsttbar and sT.v_pt > cutVpt and sT.l_pt > cutLpt and sT.pfMET > cutMET and sT.ttb_ca8_mass_pr > ttb_massCutLo and sT.ttb_ca8_mass_pr < ttb_massCutHi and sT.ttb_ca8_tau2tau1 < ttb_t2t1Cut and sT.ttb_ca8_mu < ttb_muCut and sT.jet_tau2tau1 < cut_t2t1  and sT.jet_massdrop_pr < cut_mu and isTTBar_v2 >= cutIsttbar_v2 and sT.nbjetsSSVHE <= cutNbjetsSSVHEveto and sT.jet_mass_pr > cut_mJlo and sT.jet_mass_pr < cut_mJhi and sT.mass_lvj > cut_mlvJlo and sT.mass_lvj < cut_mlvJhi and sT.ttb_ca8_ungroomed_pt > ttb_ptCut and sT.ttb_ca8_ungroomed_pt< ttb_ptCutMax and sT.deltaphi_METca8jet > cut_deltaphi_METca8jet:
                for j in range(len(obs)):
                    val = 0
                    if obs[j] == "ttb_ht": val = getattr( sT, obs[j] ) - sT.v_pt - sT.ttb_ca8_ungroomed_pt
                    else: val = getattr( sT, obs[j] );                
                    sigHistArray[j].Fill( val, getattr(sT, "totalEventWeight") );
            
        # --------------------------------------
        # plotting
        # --------------------------------------

        print "plotting..."
        for i in range(len(obs)):
            vertSliceMC = [];
            for j in range(len(mcFileNames)):
                vertSliceMC.append(mcHistArray[j][i]);
    #        print dataHistArray[i]
    #        print vertSliceMC
    #        print sigHistArray[i]
            # pass all the histograms to a drawing module
            legorientation = "right";
            DrawStackPlot(dataHistArray[i], vertSliceMC, sigHistArray[i], legorientation, g_normalizeSignalToBkg, g_plotData);
        
    ###############################################################
    ###############################################################
    if options.efficiency:

        # loop over 2 files and apply cuts before and after mass + tau2/tau1 cut, create 4 histograms 
        nbin = 7;
        ptLo = 200;
        ptHi = 900;
        h_passed_herwig = ROOT.TH1F("h_passed_herwig","h_passed_herwig",nbin,ptLo,ptHi);
        h_total_herwig = ROOT.TH1F("h_total_herwig","h_total_herwig",nbin,ptLo,ptHi);
        h_passed_pythia = ROOT.TH1F("h_passed_pythia","h_passed_pythia",nbin,ptLo,ptHi);
        h_total_pythia = ROOT.TH1F("h_total_pythia","h_total_pythia",nbin,ptLo,ptHi);
    
        infile1 = ROOT.TFile(inputDir+"ofile_WJets_Herwig.root");
        infile2 = ROOT.TFile(inputDir+"ofile_WJets_Pythia.root");
        intree1 = infile1.Get("otree");
        intree2 = infile2.Get("otree");        

        for i in range(intree1.GetEntries()):
            intree1.GetEntry(i);
            
            if getattr(intree1, "issignal") >= cutIssignal and intree1.v_pt > cutVpt and intree1.l_pt > cutLpt and intree1.pfMET > cutMET and intree1.nbjetsSSVHE <= cutNbjetsSSVHEveto and intree1.jet_mass_pr > cut_mJlo and intree1.jet_mass_pr < cut_mJhi and intree1.mass_lvj > cut_mlvJlo and intree1.mass_lvj < cut_mlvJhi:

#                h_total_herwig.Fill( getattr( intree1, "ungroomed_jet_pt" ), getattr(intree1, "totalEventWeight") );
                h_total_herwig.Fill( getattr( intree1, "ungroomed_jet_pt" ) );
                
                if intree1.jet_mass_pr > 70 and intree1.jet_mass_pr < 100 and intree1.jet_tau2tau1 < cut_t2t1  and intree1.jet_massdrop_pr < cut_mu:
                
#                    h_passed_herwig.Fill( getattr( intree1, "ungroomed_jet_pt" ), getattr(intree1, "totalEventWeight") );
                     h_passed_herwig.Fill( getattr( intree1, "ungroomed_jet_pt" ) );       

        for i in range(intree2.GetEntries()):
            intree2.GetEntry(i);

            if getattr(intree2, "issignal") >= cutIssignal and intree2.v_pt > cutVpt and intree2.l_pt > cutLpt and intree2.pfMET > cutMET and intree2.nbjetsSSVHE <= cutNbjetsSSVHEveto and intree2.jet_mass_pr > cut_mJlo and intree2.jet_mass_pr < cut_mJhi and intree2.mass_lvj > cut_mlvJlo and intree2.mass_lvj < cut_mlvJhi:
                
#                h_total_pythia.Fill( getattr( intree2, "ungroomed_jet_pt" ), getattr(intree2, "totalEventWeight") );
                h_total_pythia.Fill( getattr( intree2, "ungroomed_jet_pt" ) );
                
                if intree2.jet_mass_pr > 70 and intree2.jet_mass_pr < 100 and intree2.jet_tau2tau1 < cut_t2t1  and intree2.jet_massdrop_pr < cut_mu:
                    
#                    h_passed_pythia.Fill( getattr( intree2, "ungroomed_jet_pt" ), getattr(intree2, "totalEventWeight") );
                    h_passed_pythia.Fill( getattr( intree2, "ungroomed_jet_pt" ) );

        pass_herwig = ROOT.TEfficiency.CheckConsistency(h_passed_herwig,h_total_herwig)                        
        pass_pythia = ROOT.TEfficiency.CheckConsistency(h_passed_pythia,h_total_pythia)   
                      
        if pass_herwig and pass_pythia:
            eff_herwig = ROOT.TEfficiency(h_passed_herwig,h_total_herwig)
            eff_pythia = ROOT.TEfficiency(h_passed_pythia,h_total_pythia)

            leg = ROOT.TLegend(0.25,0.7,0.55,0.85);
            leg.SetFillStyle(0);
            leg.SetBorderSize(0);
            leg.AddEntry(eff_herwig,"W+jets, herwig", "pe");
            leg.AddEntry(eff_pythia,"W+jets, pythia", "pe");
            
            canEff = ROOT.TCanvas("canEff","canEff",800,800);
            canEff.SetGrid();            
            h0 = canEff.DrawFrame(ptLo,0.,ptHi,0.3);
            h0.SetTitle("; jet pT (GeV); fake rate");
            #ROOT.gPad.SetLogy();
            eff_herwig.SetMarkerColor(ROOT.kRed);
            eff_herwig.SetLineColor(ROOT.kRed);            
            eff_herwig.Draw("pesames");
            eff_pythia.Draw("pesames");
            leg.Draw();
            
            canEff.SaveAs("testeff.eps");

    ###############################################################
    ###############################################################
    
    if options.sigCompare:
        
        sigHistArray = [];
        sigHistArray2 = [];

        # ---------------
        # sig loop
        # ---------------
        sT = sF.Get("otree");
        sT2 = sF2.Get("otree");        
        print "sig file: ", signalFileName;
        for i in range(len(obs)):
            sigHistArray.append( ROOT.TH1F("h_"+obs[i]+"_"+sigTag,";%s;count"%(obs_xname[i]),obs_bin[i],obs_lo[i],obs_hi[i]) );
            sigHistArray2.append( ROOT.TH1F("h_"+obs[i]+"_"+sigTag2,";%s;count"%(obs_xname[i]),obs_bin[i],obs_lo[i],obs_hi[i]) );
        
        for i in range(sT.GetEntries()):
            sT.GetEntry(i);
            # set cuts
            isTTBar_v2 = 0;
            oppo1same1m = sT.ttb_nak5_same_csvl > 0 or sT.ttb_nak5_oppo_csvl > 0;
            oppo2same0m = sT.ttb_nak5_same_csvl == 0 and sT.ttb_nak5_oppo_csvl > 1;
            oppo2same0m = False;        
            if oppo1same1m or oppo2same0m:
                isTTBar_v2 = 1
            
            if getattr(sT, "issignal") >= cutIssignal and getattr(sT, "isttbar") >= cutIsttbar and sT.v_pt > cutVpt and sT.l_pt > cutLpt and sT.pfMET > cutMET and sT.ttb_ca8_mass_pr > ttb_massCutLo and sT.ttb_ca8_mass_pr < ttb_massCutHi and sT.ttb_ca8_tau2tau1 < ttb_t2t1Cut and sT.ttb_ca8_mu < ttb_muCut and sT.jet_tau2tau1 < 100.53  and sT.jet_massdrop_pr < 100.35 and isTTBar_v2 >= cutIsttbar_v2 and sT.nbjetsSSVHE <= cutNbjetsSSVHEveto and sT.jet_mass_pr > cut_mJlo and sT.jet_mass_pr < cut_mJhi and sT.mass_lvj > cut_mlvJlo and sT.mass_lvj < cut_mlvJhi and sT.ttb_ca8_ungroomed_pt > ttb_ptCut and sT.ttb_ca8_ungroomed_pt< ttb_ptCutMax and sT.ungroomed_jet_pt > 200 and sT.ungroomed_jet_pt < 300:
                for j in range(len(obs)):
                    val = 0
                    if obs[j] == "ttb_ht": val = getattr( sT, obs[j] ) - sT.v_pt - sT.ttb_ca8_ungroomed_pt
                    else: val = getattr( sT, obs[j] );                
                    sigHistArray[j].Fill( val, getattr(sT, "totalEventWeight") );
        print "sig file 2: ", signalFileName2;
        ctrpass = 0;
        ctrpassFromOutside = 0;
        mHGen = 1000;
        twosigma = 6.47E+02*2;
#        mHGen = 900;
#        twosigma = 4.49E+02*2;
#        mHGen = 800;
#        twosigma = 3.04E+02*2;
#        mHGen = 700;
#        twosigma = 1.99E+02*2;
#        mHGen = 600;
#        twosigma = 1.23E+02*2;

        for i in range(sT2.GetEntries()):
            sT2.GetEntry(i);
            # set cuts
            isTTBar_v2 = 0;
            oppo1same1m = sT2.ttb_nak5_same_csvl > 0 or sT2.ttb_nak5_oppo_csvl > 0;
            oppo2same0m = sT2.ttb_nak5_same_csvl == 0 and sT2.ttb_nak5_oppo_csvl > 1;
            oppo2same0m = False;        
            if oppo1same1m or oppo2same0m:
                isTTBar_v2 = 1
            
            if getattr(sT2, "issignal") >= cutIssignal and getattr(sT2, "isttbar") >= cutIsttbar and sT2.v_pt > cutVpt and sT2.l_pt > cutLpt and sT2.pfMET > cutMET and sT2.ttb_ca8_mass_pr > ttb_massCutLo and sT2.ttb_ca8_mass_pr < ttb_massCutHi and sT2.ttb_ca8_tau2tau1 < ttb_t2t1Cut and sT2.ttb_ca8_mu < ttb_muCut and sT2.jet_tau2tau1 < 100.53  and sT2.jet_massdrop_pr < 100.35 and isTTBar_v2 >= cutIsttbar_v2 and sT2.nbjetsSSVHE <= cutNbjetsSSVHEveto and sT2.jet_mass_pr > cut_mJlo and sT2.jet_mass_pr < cut_mJhi and sT2.mass_lvj > cut_mlvJlo and sT2.mass_lvj < cut_mlvJhi and sT2.ttb_ca8_ungroomed_pt > ttb_ptCut and sT2.ttb_ca8_ungroomed_pt< ttb_ptCutMax and sT2.ungroomed_jet_pt > 450 and sT2.ungroomed_jet_pt < 550:
                for j in range(len(obs)):
                    val = 0
                    if obs[j] == "ttb_ht": val = getattr( sT2, obs[j] ) - sT2.v_pt - sT2.ttb_ca8_ungroomed_pt
                    else: val = getattr( sT2, obs[j] );                
                    sigHistArray2[j].Fill( val, getattr(sT2, "totalEventWeight") );

            ctrpass += 1;
            if getattr(sT2, "genHMass") > (mHGen + twosigma) or getattr(sT2, "genHMass") < (mHGen - twosigma): ctrpassFromOutside += 1;
        
        print "ctrpass = ", ctrpassFromOutside,"/",ctrpass,"=",float(ctrpassFromOutside)/float(ctrpass);

        print "plotting..."
        for i in range(len(obs)):
            fit = False;
            if obs == "jet_mass_pr": fit = True;
            DrawStackPlot_sig(sigHistArray[i], sigHistArray2[i], "right", g_normalizeSignalToBkg, fit);

