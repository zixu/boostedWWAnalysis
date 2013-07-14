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
import subprocess
from subprocess import Popen
from optparse import OptionParser

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

parser.add_option('--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('--massPoint',action="store",type="int",dest="massPoint",default=-1)
parser.add_option('--cPrime',action="store",type="int",dest="cPrime",default=-1)


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

if __name__ == '__main__':


    ###############
    
    CHAN = options.channel;
    DIR = "cards_"+CHAN;

    #mass  = [ 600,1000,1200,1500,1600]
    #ccmlo = [ 500, 900,1100,1400,1500]  
    #ccmhi = [ 700,1100,1300,1600,1700]  
    #mjlo  = [  40,  40,  40,  40,  40]  
    #mjhi  = [ 130, 130, 130, 130, 130]  
    #mlo   = [ 400, 800, 800, 800, 800]      
    #mhi   = [1400,2500,2500,2500,2500]          
    #shape    = ["ExpPowExp_v1","Exp","Exp","Exp","Exp"]
    #shapeAlt = ["ExpPow2_v1"  ,"Pow","Pow","Pow","Pow"]

    mass  = [1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2400,2500]
    ccmlo = [ 900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1800,1900,2000,2200,2300]  
    ccmhi = [1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2300,2400,2600,2700]  
    mjlo  = [  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40]  
    mjhi  = [ 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130]  
    mlo   = [ 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]      
    mhi   = [2800,2800,2800,2800,2800,2800,2800,2800,2800,2800,2800,2800,2800,2800,2800]          
    #shape    = [   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN",   "ExpN"]
    shape    = [   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp",   "Exp"]
    shapeAlt = ["ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail"]
    #shape    = ["Exp","Exp","Exp","Exp","Exp"]
    #shapeAlt = ["Pow","Pow","Pow","Pow","Pow"]

    BRnew = 00;
    #cprime = [10,07,05,04,03,02,01];
    cprime = [10];
    
    moreCombineOpts = "";
    
    ###############
    
    if options.makeCards:
        if not os.path.isdir(DIR): os.system("mkdir "+DIR);
        else: 
            print "Directory "+DIR+" already exists...";
            #sys.exit(0);

    if options.computeLimits or options.plotLimits: os.chdir("cards_em");


    # put in functionality to test just one mass point or just one cprime
    nMasses = len(mass);
    mLo = 0;
    mHi = nMasses;
    nCprimes = len(cprime);
    cpLo = 0;
    cpHi = nCprimes;
    if options.massPoint > 0:   
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;
    if options.cPrime > 0:   
        curIndex = cprime.index(options.cPrime);
        cpLo = curIndex;
        cpHi = curIndex+1;

    # =====================================

    if options.makeCards:
        if not os.path.isdir("log"): os.system("mkdir log" );
        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                
                print "--------------------------------------------------";                
                print "--------------------------------------------------";                
                print "R U N N I N G   F I T S" 
                print "mass = ",mass[i],", cprime = ",cprime[j];
                print "--------------------------------------------------";                
                print "--------------------------------------------------";  
                
                time.sleep(1);
                
                #command = "nohup python EXO_doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %02d --BRnew 00  > log/log_%s_ggH%03d_%02d_%02d_%02d_%02d_%02d_%02d_%s_%s_cprime_%02d_BRnew_00 &"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j]);
                command = "nohup python exo_doFit_class.py %s BulkG_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %02d --BRnew 00  &"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j]);
                #command = "python doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %02d --BRnew 00  >> log/log_%s_ggH%03d_%02d_%02d_%02d_%02d_%02d_%02d_%s_%s_cprime_%02d_BRnew_00 "%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j]);
                print command #raw_input("ENTER");
                os.system(command);
                
                #mvcmmd1 = "mv plots_"+CHAN+"* "+DIR+"/plots_"+CHAN+"_"+str(mass[i])+"_"+str(cprime[j])+"_00";
                ##print mvcmmd1;
                #os.system(mvcmmd1);
            
            #mvcmmd0 = "mv hwwlvj_*_"+CHAN+"* other_*_"+CHAN+"* "+DIR+"/.";
            #print mvcmmd0;
            #os.system(mvcmmd0);


    # =====================================

    if options.computeLimits:

        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):

                print "--------------------------------------------------";
                print "--------------------------------------------------";                
                print "creating card: hwwlvj_ggH%03d_em_%02d_00_unbin.txt"%(mass[i],cprime[j]);
                combineCmmd = "combineCards.py hwwlvj_ggH%03d_el_%02d_00_unbin.txt hwwlvj_ggH%03d_mu_%02d_00_unbin.txt > hwwlvj_ggH%03d_em_%02d_00_unbin.txt"%(mass[i],cprime[j],mass[i],cprime[j],mass[i],cprime[j]);
                os.system(combineCmmd);

                
                print "running card: hwwlvj_ggH%03d_em_%02d_00_unbin.txt"%(mass[i],cprime[j]);
                runCmmd = "combine -M Asymptotic -n hwwlvj_ggH%03d_em_%02d_00_unbin -m %03d -d hwwlvj_ggH%03d_em_%02d_00_unbin.txt %s"%(mass[i],cprime[j],mass[i],mass[i],cprime[j],moreCombineOpts);
                os.system(runCmmd);

    # =====================================
    
    nGraphs = nCprimes*2 + 2;

    tGraphs = [];
    nPoints = len(mass);
    if options.plotLimits:

#        for i in range(mLo,mHi):
#            for j in range(cpLo,cpHi):
        for j in range(cpLo,cpHi):
            
            xbins = array('d', [])
            ybins_exp = array('d', [])
            ybins_obs = array('d', [])            
            xbins_env = array('d', [])
            ybins_1s = array('d', [])                        
            ybins_2s = array('d', [])                        
            
            for i in range(mLo,mHi):

                curFile = "higgsCombinehwwlvj_ggH%03d_em_%02d_00_unbin.Asymptotic.mH%03d.root"%(mass[i],cprime[j],mass[i]);
                curAsymLimits = getAsymLimits(curFile);
                
                xbins.append( mass[i] );
                ybins_exp.append( curAsymLimits[3] );                
                ybins_obs.append( curAsymLimits[0] );                                
                xbins_env.append( mass[i] );                                
                ybins_2s.append( curAsymLimits[1] );                                
                ybins_1s.append( curAsymLimits[2] );                                                

                print "mass: ", mass[i], ", cprime: ", cprime[j], " -- obs: ", curAsymLimits[0], ", exp: ", curAsymLimits[3];

            for i in range(mHi-1,mLo-1,-1):
                curFile = "higgsCombinehwwlvj_ggH%03d_em_%02d_00_unbin.Asymptotic.mH%03d.root"%(mass[i],cprime[j],mass[i]);
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

            tGraphs.append(curGraph_exp);
            tGraphs.append(curGraph_obs);
            tGraphs.append(curGraph_1s);
            tGraphs.append(curGraph_2s);

#        MakePlots(tGraphs);
        leg = ROOT.TLegend(0.15,0.7,0.85,0.9);
        leg.SetFillStyle(0);
        leg.SetBorderSize(0);
        leg.SetNColumns(2);
        for k in range(nCprimes):
            if k == 0: leg.AddEntry(tGraphs[0],"expected, C' = 1","L")
            tmplabel = "observed, C' = %1.1f"%( float(cprime[k]/10.) )
            leg.AddEntry(tGraphs[k*4+1],tmplabel,"L");        
                
                
        drawColors = [];
        drawColors.append(ROOT.kRed+2);
        drawColors.append(ROOT.kGreen+2);
        drawColors.append(ROOT.kBlue+2);
        drawColors.append(ROOT.kMagenta+2);
        drawColors.append(ROOT.kCyan+2);
        drawColors.append(ROOT.kGray+2);                          
        drawColors.append(ROOT.kOrange+2);                                          
                
        oneLine = ROOT.TF1("oneLine","1",599,1001);
        oneLine.SetLineColor(ROOT.kRed);
        oneLine.SetLineWidth(3);

        can = ROOT.TCanvas("can","can",800,800);
        hrl = can.DrawFrame(599,0.5,1001,100.0);
        hrl.GetYaxis().SetTitle("#mu' = C' #times #mu");
        hrl.GetXaxis().SetTitle("mass (GeV)");
        can.SetGrid();
        tGraphs[3].Draw("f");
        tGraphs[2].Draw("f");
        tGraphs[1].Draw("PL");
        tGraphs[0].Draw("PL");
            
        oneLine.Draw("LSAMES");

        for k in range(1,nCprimes):
#            tGraphs[k*4].SetLineColor(k+1);
#            tGraphs[k*4+1].SetLineColor(k+1);
#            tGraphs[k*4].SetMarkerColor(k+1);
#            tGraphs[k*4+1].SetMarkerColor(k+1);
            tGraphs[k*4].SetLineColor(drawColors[k]);
            tGraphs[k*4+1].SetLineColor(drawColors[k]);
            tGraphs[k*4].SetMarkerColor(drawColors[k]);
            tGraphs[k*4+1].SetMarkerColor(drawColors[k]);
            tGraphs[k*4].Draw("PL");
            tGraphs[k*4+1].Draw("PL");   
                         
        leg.Draw();

        ROOT.gPad.SetLogy();
        can.SaveAs("test.eps");



