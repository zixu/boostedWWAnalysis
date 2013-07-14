#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

from array import array

import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex
from ROOT import *
import subprocess
from subprocess import Popen
from optparse import OptionParser

#from condorUtils import submitBatchJob

ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.16);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadTickX(1);
ROOT.gStyle.SetPadTickY(1); 
############################################################
############################################
#            Job steering                  #
############################################
parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')

parser.add_option('--makeCards', action='store_true', dest='makeCards', default=False,
                  help='no X11 windows')
parser.add_option('--computeLimits', action='store_true', dest='computeLimits', default=False,
                  help='no X11 windows')
parser.add_option('--plotLimits', action='store_true', dest='plotLimits', default=False,
                  help='no X11 windows')

# submit jobs to condor
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False,
                  help='no X11 windows')


parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--massPoint',action="store",type="int",dest="massPoint",default=-1)
parser.add_option('--cPrime',action="store",type="int",dest="cPrime",default=-1)
parser.add_option('--brNew',action="store",type="int",dest="brNew",default=-1)
parser.add_option('--odir',action="store",type="string",dest="odir",default=".")
parser.add_option('--sigChannel',action="store",type="string",dest="sigChannel",default='')


(options, args) = parser.parse_args()
############################################################

############################################################
def getAsymLimits(file):
    
    
    f = ROOT.TFile(file);
    t = f.Get("limit");
    entries = t.GetEntries();
    
    lims = [0,0,0,0,0,0];
    
    for i in range(entries):
        
        t.GetEntry(i);
        t_quantileExpected = t.quantileExpected;
        t_limit = t.limit;
        
        #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected;
        
        if t_quantileExpected == -1.: lims[0] = t_limit;
        elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit;
        elif t_quantileExpected >= 0.15 and t_quantileExpected <= 0.17: lims[2] = t_limit;            
        elif t_quantileExpected == 0.5: lims[3] = t_limit;            
        elif t_quantileExpected >= 0.83 and t_quantileExpected <= 0.85: lims[4] = t_limit;
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit;
        else: print "Unknown quantile!"
    
    return lims;

############################################################
############################################################
############################################################

def makeSMLimits(SIGCH):

    mass  = [ 600, 700, 800, 900,1000]    
    cprime = 10;
    brnew = 00;
    
    xbins = array('d', [])
    xbins_env = array('d', [])
    ybins_exp = array('d', [])
    ybins_obs = array('d', [])            
    ybins_1s = array('d', [])                        
    ybins_2s = array('d', [])    
    
    for i in range(len(mass)):
        curFile = "higgsCombinehwwlvj_ggH%03d_em%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[i],SIGCH,cprime,brnew,mass[i]);
        #print "curFile: ",curFile
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );                                
        ybins_exp.append( curAsymLimits[3] );                
        ybins_obs.append( curAsymLimits[0] );                                
        ybins_2s.append( curAsymLimits[1] );                                
        ybins_1s.append( curAsymLimits[2] ); 

    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombinehwwlvj_ggH%03d_em%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[i],SIGCH,cprime,brnew,mass[i]);
        curAsymLimits = getAsymLimits(curFile);   
        xbins_env.append( mass[i] );                                
        ybins_2s.append( curAsymLimits[5] );                                
        ybins_1s.append( curAsymLimits[4] ); 
    

    curGraph_exp = ROOT.TGraph(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraph(nPoints,xbins,ybins_obs);
    curGraph_1s = ROOT.TGraph(nPoints*2,xbins_env,ybins_1s);
    curGraph_2s = ROOT.TGraph(nPoints*2,xbins_env,ybins_2s);
    
    curGraph_exp.SetLineStyle(2);
    curGraph_exp.SetLineWidth(2);
    curGraph_obs.SetLineWidth(2);                    
    curGraph_exp.SetMarkerSize(2);
    curGraph_obs.SetMarkerSize(2);                    
    curGraph_1s.SetFillColor(ROOT.kGreen);
    curGraph_2s.SetFillColor(ROOT.kYellow);

    # -------
    banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
    banner.SetNDC(); banner.SetTextSize(0.028);
    oneLine = ROOT.TF1("oneLine","1",599,1001);
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    can_SM = ROOT.TCanvas("can_SM","can_SM",1000,800);
    hrl_SM = can_SM.DrawFrame(599,0.0,1001,10.0);
    hrl_SM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_SM.GetYaxis().SetTitleOffset(1.4);
    hrl_SM.GetXaxis().SetTitle("Higgs boson mass (GeV/c^{2})");
    #can_SM.SetGrid();
    
    curGraph_2s.Draw("F");
    curGraph_1s.Draw("F");
    curGraph_obs.Draw("PL");
    curGraph_exp.Draw("PL");
    
    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.AddEntry(curGraph_obs,"Observed","L")
    leg2.AddEntry(curGraph_exp,"Expected","L")
    leg2.AddEntry(curGraph_1s,"Expected, #pm 1#sigma","F")
    leg2.AddEntry(curGraph_2s,"Expected, #pm 2#sigma","F")
    leg2.AddEntry(oneLine,"SM Expected","L")

    leg2.Draw();
    banner.Draw();
    oneLine.Draw("LESAMES");
    
    #ROOT.gPad.SetLogy();
    can_SM.SaveAs("limitFigs/SMLim%s.eps"%(SIGCH));                      
    can_SM.SaveAs("limitFigs/SMLim%s.png"%(SIGCH));                      
    can_SM.SaveAs("limitFigs/SMLim%s.pdf"%(SIGCH));                      

############################################################
############################################################
############################################################

def makeBSMLimits_vsMass( SIGCH, cprimes ):

    print "module ===> makeBSMLimits_vsMass";
    
    mass  = [ 600, 700, 800, 900,1000]    
    curcolors = [1,2,4,6]
    brnew = 00;
    
    massCS  = [];
    if SIGCH == "":
        massCS.append( (0.5230 + 0.09688) );
        massCS.append( (0.2288 + 0.06330) );
        massCS.append( (0.1095 + 0.04365) );
        massCS.append( (0.05684 + 0.03164) );
        massCS.append( (0.03163 + 0.02399) );
    elif SIGCH == "_ggH":
        massCS.append( (0.5230) );
        massCS.append( (0.2288) );
        massCS.append( (0.1095) );
        massCS.append( (0.05684) );
        massCS.append( (0.03163) );
    elif SIGCH == "_vbfH":
        massCS.append( (0.09688) );
        massCS.append( (0.06330) );
        massCS.append( (0.04365) );
        massCS.append( (0.03164) );
        massCS.append( (0.02399) );
    else:
        print "problem!"
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01]    
    
    gridMax = -999;
    gridMaxSig = -999;

    tGraphs_exp = [];
    tGraphs_obs = [];    
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];    
    tGraphs_csXbr_th = [];    
    
    for j in range(len(cprimes)):
    
        xbins = array('d', [])
        ybins_exp = array('d', [])
        ybins_obs = array('d', [])            
        ybins_csXbr_exp = array('d', [])
        ybins_csXbr_obs = array('d', [])            
        ybins_csXbr_th = array('d', [])            
        
        for i in range(len(mass)):
            curFile = "higgsCombinehwwlvj_ggH%03d_em%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[i],SIGCH,cprimes[j],brnew,mass[i]);
            #print "curFile: ",curFile
            curAsymLimits = getAsymLimits(curFile);
            xbins.append( mass[i] );
            ybins_exp.append( curAsymLimits[3] );                
            ybins_obs.append( curAsymLimits[0] );    
            ybins_csXbr_exp.append( curAsymLimits[3]*massCS[i]*cprimes[j]*0.1*(1-brnew*0.1)*massBRWW[i] );                
            ybins_csXbr_obs.append( curAsymLimits[0]*massCS[i]*cprimes[j]*0.1*(1-brnew*0.1)*massBRWW[i] );    
            ybins_csXbr_th.append( 1.*massCS[i]*cprimes[j]*0.1*(1-brnew*0.1)*massBRWW[i] );    
        
            #print "curAsymLimits[0]: ",curAsymLimits[0]
            
            if gridMax < curAsymLimits[3]: gridMax = curAsymLimits[3];
            cscur = ( curAsymLimits[3]*massCS[i]*cprimes[j]*0.1*(1-brnew*0.1)*massBRWW[i] );
            if gridMaxSig < cscur: gridMaxSig = cscur;
        
        curGraph_exp = ROOT.TGraph(nPoints,xbins,ybins_exp);
        curGraph_obs = ROOT.TGraph(nPoints,xbins,ybins_obs);        
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(2);
        curGraph_obs.SetLineWidth(2);                    
        curGraph_exp.SetMarkerSize(2);
        curGraph_obs.SetMarkerSize(2);

        curGraph_csXbr_exp = ROOT.TGraph(nPoints,xbins,ybins_csXbr_exp);
        curGraph_csXbr_obs = ROOT.TGraph(nPoints,xbins,ybins_csXbr_obs);        
        curGraph_csXbr_th = ROOT.TGraph(nPoints,xbins,ybins_csXbr_th);
        curGraph_csXbr_exp.SetLineStyle(2);
        curGraph_csXbr_exp.SetLineWidth(2);
        curGraph_csXbr_obs.SetLineWidth(2);                    
        curGraph_csXbr_exp.SetMarkerSize(2);
        curGraph_csXbr_obs.SetMarkerSize(2);
        curGraph_csXbr_th.SetLineWidth(2);        

        tGraphs_exp.append( curGraph_exp );
        tGraphs_obs.append( curGraph_obs );
        tGraphs_csXbr_exp.append( curGraph_csXbr_exp );
        tGraphs_csXbr_obs.append( curGraph_csXbr_obs );
        tGraphs_csXbr_th.append( curGraph_csXbr_th );

    # -------
    banner = TLatex(0.47,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
    banner.SetNDC(); banner.SetTextSize(0.028);

    banner2 = TLatex(0.17,0.91,("BR_{new} = 0"));
    banner2.SetNDC(); banner2.SetTextSize(0.028);
    
    can_BSM = ROOT.TCanvas("can_BSM","can_BSM",1000,800);
    hrl_BSM = can_BSM.DrawFrame(599,0.0,1001,gridMax*1.5);
    hrl_BSM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_BSM.GetYaxis().SetTitleOffset(1.4);
    hrl_BSM.GetXaxis().SetTitle("Higgs boson mass (GeV/c^{2})");
    #can_BSM.SetGrid();
    can_BSM.SetRightMargin(0.1);
    
    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetNColumns(2);

    for k in range(len(cprimes)):
        #print cprime[k]
        tGraphs_exp[k].SetLineStyle(2);
        tGraphs_exp[k].SetLineColor(curcolors[k]);
        tGraphs_obs[k].SetLineColor(curcolors[k]);                
        tGraphs_exp[k].SetLineWidth(2);
        tGraphs_obs[k].SetLineWidth(2);                
        tGraphs_exp[k].Draw("PL");
        tGraphs_obs[k].Draw("PL");   

        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_exp[k],tmplabel,"L")        
        tmplabel = "obs., C'^{ 2} = %1.1f"%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_obs[k],tmplabel,"L");        
    
    leg2.Draw();
    banner.Draw();
    banner2.Draw();

    #ROOT.gPad.SetLogy();
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu.eps"%(SIGCH));                      
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu.png"%(SIGCH));                      
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu.pdf"%(SIGCH));          

    ##----

    can_BSMsig = ROOT.TCanvas("can_BSMsig","can_BSMsig",1000,800);
    hrl_BSMsig = can_BSMsig.DrawFrame(599,0.0,1001,gridMaxSig*1.8);
    hrl_BSMsig.GetYaxis().SetTitle("#sigma_{95%} #times BR_{WW} (pb)");
    hrl_BSMsig.GetYaxis().SetTitleOffset(1.4);
    hrl_BSMsig.GetXaxis().SetTitle("Higgs boson mass (GeV/c^{2})");
    #can_BSMsig.SetGrid();
    can_BSMsig.SetRightMargin(0.1);

    leg2 = ROOT.TLegend(0.2,0.65,0.85,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetNColumns(3);

    for k in range(len(cprimes)):
        tGraphs_csXbr_exp[k].SetLineStyle(2);
        tGraphs_csXbr_exp[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_obs[k].SetLineColor(curcolors[k]);                
        tGraphs_csXbr_th[k].SetLineStyle(3);
        tGraphs_csXbr_th[k].SetLineColor(curcolors[k]);                        
        tGraphs_csXbr_exp[k].SetLineWidth(2);
        tGraphs_csXbr_obs[k].SetLineWidth(2);                
        tGraphs_csXbr_exp[k].Draw("PL");
        tGraphs_csXbr_obs[k].Draw("PL");   
        tGraphs_csXbr_th[k].Draw("PL");   
        
        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_exp[k],tmplabel,"L")        
        tmplabel = "obs., C'^{ 2} = %1.1f    "%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_obs[k],tmplabel,"L");        
        tmplabel = "th., C'^{ 2} = %1.1f"%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_th[k],tmplabel,"L");        

    leg2.Draw();
    banner.Draw();
    banner2.Draw();

    #ROOT.gPad.SetLogy();
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma.eps"%(SIGCH));                      
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma.png"%(SIGCH));                      
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma.pdf"%(SIGCH));  

############################################################
############################################################
############################################################

def makeBSMLimits_vsBRnew( SIGCH, cprimes, mass ):
    
    print "module ===> makeBSMLimits_vsBRnew";
    
    curcolors = [1,2,4,6]
    brnews = [00,01,02,03,04,05];
    massindex = {600:0,700:1,800:2,900:3,1000:4}
    massCS  = [];
    if SIGCH == "":
        massCS.append( (0.5230 + 0.09688) );
        massCS.append( (0.2288 + 0.06330) );
        massCS.append( (0.1095 + 0.04365) );
        massCS.append( (0.05684 + 0.03164) );
        massCS.append( (0.03163 + 0.02399) );
    elif SIGCH == "_ggH":
        massCS.append( (0.5230) );
        massCS.append( (0.2288) );
        massCS.append( (0.1095) );
        massCS.append( (0.05684) );
        massCS.append( (0.03163) );
    elif SIGCH == "_vbfH":
        massCS.append( (0.09688) );
        massCS.append( (0.06330) );
        massCS.append( (0.04365) );
        massCS.append( (0.03164) );
        massCS.append( (0.02399) );
    else:
        print "problem!"
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01]    
    
    gridMax = -999;
    gridMaxSig = -999;
    
    tGraphs_exp = [];
    tGraphs_obs = [];    
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];    
    tGraphs_csXbr_th = [];    
    
    for j in range(len(cprimes)):
        
        xbins = array('d', [])
        ybins_exp = array('d', [])
        ybins_obs = array('d', [])            
        ybins_csXbr_exp = array('d', [])
        ybins_csXbr_obs = array('d', [])            
        ybins_csXbr_th = array('d', [])            
        
        for i in range(len(brnews)):
            curFile = "higgsCombinehwwlvj_ggH%03d_em%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass,SIGCH,cprimes[j],brnews[i],mass);
            #print "curFile: ",curFile
            curAsymLimits = getAsymLimits(curFile);
            xbins.append( brnews[i] );
            ybins_exp.append( curAsymLimits[3] );                
            ybins_obs.append( curAsymLimits[0] );    
            ybins_csXbr_exp.append( curAsymLimits[3]*massCS[massindex[mass]]*cprimes[j]*0.1*(1-brnews[i]*0.1)*massBRWW[massindex[mass]] );                
            ybins_csXbr_obs.append( curAsymLimits[0]*massCS[massindex[mass]]*cprimes[j]*0.1*(1-brnews[i]*0.1)*massBRWW[massindex[mass]] );    
            ybins_csXbr_th.append( 1.*massCS[massindex[mass]]*cprimes[j]*0.1*(1-brnews[i]*0.1)*massBRWW[massindex[mass]] );    
            
            #print "curAsymLimits[0]: ",curAsymLimits[0]
            
            if gridMax < curAsymLimits[3]: gridMax = curAsymLimits[3];
            cscur = ( curAsymLimits[3]*massCS[massindex[mass]]*cprimes[j]*0.1*(1-brnews[i]*0.1)*massBRWW[massindex[mass]] );
            if gridMaxSig < cscur: gridMaxSig = cscur;
        
        curGraph_exp = ROOT.TGraph(nPoints,xbins,ybins_exp);
        curGraph_obs = ROOT.TGraph(nPoints,xbins,ybins_obs);        
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(2);
        curGraph_obs.SetLineWidth(2);                    
        curGraph_exp.SetMarkerSize(2);
        curGraph_obs.SetMarkerSize(2);
        
        curGraph_csXbr_exp = ROOT.TGraph(nPoints,xbins,ybins_csXbr_exp);
        curGraph_csXbr_obs = ROOT.TGraph(nPoints,xbins,ybins_csXbr_obs);        
        curGraph_csXbr_th = ROOT.TGraph(nPoints,xbins,ybins_csXbr_th);
        curGraph_csXbr_exp.SetLineStyle(2);
        curGraph_csXbr_exp.SetLineWidth(2);
        curGraph_csXbr_obs.SetLineWidth(2);                    
        curGraph_csXbr_exp.SetMarkerSize(2);
        curGraph_csXbr_obs.SetMarkerSize(2);
        curGraph_csXbr_th.SetLineWidth(2);        
        
        tGraphs_exp.append( curGraph_exp );
        tGraphs_obs.append( curGraph_obs );
        tGraphs_csXbr_exp.append( curGraph_csXbr_exp );
        tGraphs_csXbr_obs.append( curGraph_csXbr_obs );
        tGraphs_csXbr_th.append( curGraph_csXbr_th );
    
    # -------
    banner = TLatex(0.47,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
    banner.SetNDC(); banner.SetTextSize(0.028);

    banner2 = TLatex(0.17,0.91,("Higgs mass, %i GeV/c^{2}"%(mass)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);
    
    can_BSM = ROOT.TCanvas("can_BSM","can_BSM",1000,800);
    hrl_BSM = can_BSM.DrawFrame(-0.01,0.0,0.51,gridMax*1.5);
    hrl_BSM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_BSM.GetYaxis().SetTitleOffset(1.4);
    hrl_BSM.GetXaxis().SetTitle("BR_{new}");
    #can_BSM.SetGrid();
    can_BSM.SetRightMargin(0.1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetNColumns(2);
    
    for k in range(len(cprimes)):
        tGraphs_exp[k].SetLineStyle(2);
        tGraphs_exp[k].SetLineColor(curcolors[k]);
        tGraphs_obs[k].SetLineColor(curcolors[k]);                
        tGraphs_exp[k].SetLineWidth(2);
        tGraphs_obs[k].SetLineWidth(2);                
        tGraphs_exp[k].Draw("PL");
        tGraphs_obs[k].Draw("PL");   
        
        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_exp[k],tmplabel,"L")        
        tmplabel = "obs., C'^{ 2} = %1.1f"%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_obs[k],tmplabel,"L");        
    
    leg2.Draw();
    banner.Draw();
    banner2.Draw();

    #ROOT.gPad.SetLogy();
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.eps"%(SIGCH,mass));                      
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.png"%(SIGCH,mass));                      
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.pdf"%(SIGCH,mass));          
    
    ##----
    
    can_BSMsig = ROOT.TCanvas("can_BSMsig","can_BSMsig",1000,800);
    hrl_BSMsig = can_BSMsig.DrawFrame(-0.01,0.0,0.51,gridMaxSig*1.8);
    hrl_BSMsig.GetYaxis().SetTitle("#sigma_{95%} #times BR_{WW} (pb)");
    hrl_BSMsig.GetYaxis().SetTitleOffset(1.4);
    hrl_BSMsig.GetXaxis().SetTitle("BR_{new}");
    #can_BSMsig.SetGrid();
    can_BSMsig.SetRightMargin(0.1);
    
    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.SetNColumns(2);
    
    for k in range(len(cprimes)):
        tGraphs_csXbr_exp[k].SetLineStyle(2);
        tGraphs_csXbr_exp[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_obs[k].SetLineColor(curcolors[k]);                
        tGraphs_csXbr_th[k].SetLineStyle(3);
        tGraphs_csXbr_th[k].SetLineColor(curcolors[k]);                        
        tGraphs_csXbr_exp[k].SetLineWidth(2);
        tGraphs_csXbr_obs[k].SetLineWidth(2);                
        tGraphs_csXbr_exp[k].Draw("PL");
        tGraphs_csXbr_obs[k].Draw("PL");   
#        tGraphs_csXbr_th[k].Draw("PL");   
        
        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_exp[k],tmplabel,"L")        
        tmplabel = "obs., C'^{ 2} = %1.1f    "%( float((cprimes[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_obs[k],tmplabel,"L");        
#        tmplabel = "th., C'^{ 2} = %1.1f"%( float((cprimes[k])/10.) )
#        leg2.AddEntry(tGraphs_csXbr_th[k],tmplabel,"L");        
    
    leg2.Draw();
    banner.Draw();
    banner2.Draw();

    #ROOT.gPad.SetLogy();
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.eps"%(SIGCH,mass));                      
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.png"%(SIGCH,mass));                      
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.pdf"%(SIGCH,mass));  


############################################################
############################################################
############################################################

def makeBSMLimits_2D( SIGCH, mass ):
    
    print "module ===> makeBSMLimits_2D";
    
    cprimes = [1,2,3,4,5,6,7,8,9,10];
    brnews = [00,01,02,03,04,05];
    massindex = {600:0,700:1,800:2,900:3,1000:4}
    massCS  = [];
    if SIGCH == "":
        massCS.append( (0.5230 + 0.09688) );
        massCS.append( (0.2288 + 0.06330) );
        massCS.append( (0.1095 + 0.04365) );
        massCS.append( (0.05684 + 0.03164) );
        massCS.append( (0.03163 + 0.02399) );
    elif SIGCH == "_ggH":
        massCS.append( (0.5230) );
        massCS.append( (0.2288) );
        massCS.append( (0.1095) );
        massCS.append( (0.05684) );
        massCS.append( (0.03163) );
    elif SIGCH == "_vbfH":
        massCS.append( (0.09688) );
        massCS.append( (0.06330) );
        massCS.append( (0.04365) );
        massCS.append( (0.03164) );
        massCS.append( (0.02399) );
    else:
        print "problem!"
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01]    
            
    h2d_exp = TH2D("h2d_exp",";BR_{new};C'^{ 2};#mu = #sigma_{95%} / #sigma_{SM}",6,-0.05,0.55,10,0.05,1.05);
    h2d_obs = TH2D("h2d_obs",";BR_{new};C'^{ 2};#mu = #sigma_{95%} / #sigma_{SM}",6,-0.05,0.55,10,0.05,1.05);
    h2d_csXbr_exp = TH2D("h2d_csXbr_exp",";BR_{new};C'^{ 2};#sigma_{95%} #times BR_{WW} (pb)",6,-0.05,0.55,10,0.05,1.05);
    h2d_csXbr_obs = TH2D("h2d_csXbr_obs",";BR_{new};C'^{ 2};#sigma_{95%} #times BR_{WW} (pb)",6,-0.05,0.55,10,0.05,1.05);    
    
    for j in range(len(cprimes)):
        for i in range(len(brnews)):
            curFile = "higgsCombinehwwlvj_ggH%03d_em%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass,SIGCH,cprimes[j],brnews[i],mass);
            curAsymLimits = getAsymLimits(curFile);
            h2d_exp.SetBinContent( i+1, j+1, curAsymLimits[3] );
            h2d_obs.SetBinContent( i+1, j+1, curAsymLimits[0] );
            h2d_csXbr_exp.SetBinContent( i+1, j+1, curAsymLimits[3]*massCS[massindex[mass]]*cprimes[j]*0.1*(1-brnews[i]*0.1)*massBRWW[massindex[mass]] );
            h2d_csXbr_obs.SetBinContent( i+1, j+1, curAsymLimits[0]*massCS[massindex[mass]]*cprimes[j]*0.1*(1-brnews[i]*0.1)*massBRWW[massindex[mass]] );

    h2d_exp.GetZaxis().SetTitleOffset(1.4);
    h2d_obs.GetZaxis().SetTitleOffset(1.4);
    h2d_csXbr_exp.GetZaxis().SetTitleOffset(1.4);
    h2d_csXbr_obs.GetZaxis().SetTitleOffset(1.4);

    h2d_exp.GetZaxis().RotateTitle(1);
    h2d_obs.GetZaxis().RotateTitle(1);
    h2d_csXbr_exp.GetZaxis().RotateTitle(1);
    h2d_csXbr_obs.GetZaxis().RotateTitle(1);

    # -------
    banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
    banner.SetNDC(); banner.SetTextSize(0.028);

    banner2 = TLatex(0.17,0.91,("Higgs mass, %i GeV"%(mass)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    can1_BSM2D = ROOT.TCanvas("can1_BSM2D","can1_BSM2D",1000,800);
    h2d_exp.Draw("colz");
    banner.Draw();
    banner2.Draw();
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.eps"%(SIGCH,mass));                      
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.png"%(SIGCH,mass));                      
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.pdf"%(SIGCH,mass));          

    can2_BSM2D = ROOT.TCanvas("can2_BSM2D","can2_BSM2D",1000,800);
    h2d_obs.Draw("colz");
    banner.Draw();
    banner2.Draw();
    can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.eps"%(SIGCH,mass));                      
    can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.png"%(SIGCH,mass));                      
    can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.pdf"%(SIGCH,mass));          
    
    can3_BSM2D = ROOT.TCanvas("can3_BSM2D","can3_BSM2D",1000,800);
    h2d_csXbr_exp.Draw("colz");
    banner.Draw();
    banner2.Draw();
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_%i.eps"%(SIGCH,mass));                      
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_%i.png"%(SIGCH,mass));                      
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_%i.pdf"%(SIGCH,mass));          

    
    can4_BSM2D = ROOT.TCanvas("can4_BSM2D","can4_BSM2D",1000,800);
    h2d_csXbr_obs.Draw("colz");
    banner.Draw();
    banner2.Draw();
    can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_%i.eps"%(SIGCH,mass));                      
    can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_%i.png"%(SIGCH,mass));                      
    can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_%i.pdf"%(SIGCH,mass));          

    ##----



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
################################ M A I N #################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

if __name__ == '__main__':
    
    
    ###############
    
    CHAN = options.channel;
    DIR = "cards_"+CHAN;
    SIGCH = "";
    if options.sigChannel.find("H") >= 0: SIGCH = "_"+options.sigChannel;
        
    
    mass  = [ 600, 700, 800, 900,1000]
    ccmlo = [ 550, 600, 700, 750, 800]  
    ccmhi = [ 700, 850, 950,1100,1200]  
    mjlo  = [  40,  40,  40,  40,  40]  
    mjhi  = [ 130, 130, 130, 130, 130]  
    mlo   = [ 400, 400, 600, 600, 600]      
    mhi   = [1000,1000,1400,1400,1400]          
    shape    = ["ErfPowExp_v1","ErfPowExp_v1","Exp","Exp","Exp"]
    shapeAlt = [  "ErfPow2_v1",  "ErfPow2_v1","Pow","Pow","Pow"]
        
    BRnew = [0,1,2,3,4,5];
    cprime = [1,2,3,4,5,6,7,8,9,10];    
    massCS  = [];
    massCS.append( (0.5230 + 0.09688) );
    massCS.append( (0.2288 + 0.06330) );
    massCS.append( (0.1095 + 0.04365) );
    massCS.append( (0.05684 + 0.03164) );
    massCS.append( (0.03163 + 0.02399) );
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01]

    moreCombineOpts = "";
    
    ###############

    if options.computeLimits or options.plotLimits: os.chdir("cards_em");

    # put in functionality to test just one mass point or just one cprime
    nMasses = len(mass);
    mLo = 0;
    mHi = nMasses;
    nCprimes = len(cprime);
    cpLo = 0;
    cpHi = nCprimes;
    nBRnews = len(BRnew);
    brLo = 0;
    brHi = nBRnews;
    
    if options.massPoint > 0:   
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;
    if options.cPrime > 0:   
        curIndex = cprime.index(options.cPrime);
        cpLo = curIndex;
        cpHi = curIndex+1;
        nCprimes = 1;
    if options.brNew >= 0:   
        curIndex = BRnew.index(options.brNew);
        brLo = curIndex;
        brHi = curIndex+1;
        nBRnews = 1;  


# =====================================
    
    nGraphs = nCprimes*2 + 2;
    
    tGraphs = [];
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];
    tGraphs_csXbr_vsBRnew_exp = [];
    tGraphs_csXbr_vsBRnew_obs = [];    
    tGraphs_csXbr_obs = [];

    nPoints = len(mass);
    
    if options.plotLimits:
        
        makeSMLimits( SIGCH );
        cprimes = [3,6,8,10];    
        makeBSMLimits_vsMass( SIGCH, cprimes );
        makeBSMLimits_vsBRnew( SIGCH, cprimes, 600 );
        makeBSMLimits_vsBRnew( SIGCH, cprimes, 700 );
        makeBSMLimits_vsBRnew( SIGCH, cprimes, 800 );
        makeBSMLimits_vsBRnew( SIGCH, cprimes, 900 );
        makeBSMLimits_vsBRnew( SIGCH, cprimes, 1000 );

        makeBSMLimits_2D( SIGCH, 600 );
        
