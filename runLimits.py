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

#from condorUtils import submitBatchJob

ROOT.gStyle.SetPadRightMargin(0.16);

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
parser.add_option('--sigChannel',action="store",type="string",dest="sigChannel",default="")


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

def submitBatchJob( command, fn ):
    
    currentDir = os.getcwd();
    
    # create a dummy bash/csh
    outScript=open(fn+".sh","w");
    
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'date');
    outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
    outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');    
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd -');    
    outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
    outScript.write("\n"+'echo ${PATH}');
    outScript.write("\n"+'ls'); 
    outScript.write("\n"+command);
    outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *');    
    outScript.close();
    
    # link a condor script to your shell script
    condorScript=open("condor_"+fn,"w");
    condorScript.write('universe = vanilla')
    condorScript.write("\n"+"Executable = "+fn+".sh")
    condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
    condorScript.write("\n"+'Should_Transfer_Files = YES')
    condorScript.write("\n"+'Transfer_Input_Files = doFit_class.py')    
    condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log    = out_$(Cluster).log')
    condorScript.write("\n"+'Notification    = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();
    
    # submit the condor job 
    
    os.system("condor_submit "+"condor_"+fn)

# ----------------------------------------
def submitBatchJobCombine( command, fn, mass, cprime, BRnew ):
    
    
    SIGCH = "";
    if options.sigChannel.find("H") >= 0: SIGCH = "_"+options.sigChannel;
    
    currentDir = os.getcwd();
    
    # create a dummy bash/csh
    outScript=open(fn+".sh","w");
    
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'date');
    outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd -');
    outScript.write("\n"+'ls');    
    outScript.write("\n"+command);
    outScript.close();
    
    file1 = "hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass,SIGCH,cprime,BRnew);
    file2 = "hwwlvj_ggH%03d_mu_%02d_%02d_workspace.root"%(mass,cprime,BRnew);
    file3 = "hwwlvj_ggH%03d_el_%02d_%02d_workspace.root"%(mass,cprime,BRnew);    
    # link a condor script to your shell script
    condorScript=open("subCondor_"+fn,"w");
    condorScript.write('universe = vanilla')
    condorScript.write("\n"+"Executable = "+fn+".sh")
    condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
    condorScript.write("\n"+'Should_Transfer_Files = YES')
    condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2+', '+file3)    
    condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log    = out_$(Cluster).log')
    condorScript.write("\n"+'Notification    = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();
    
    # submit the condor job 
    
    os.system("condor_submit "+"subCondor_"+fn)

############################################################

if __name__ == '__main__':


    ###############
    
    CHAN = options.channel;
    DIR = "cards_"+CHAN;
    SIGCH = "";
    if options.sigChannel.find("H") >= 0: SIGCH = "_"+options.sigChannel;

    #mass  = [ 600, 700, 800, 900,1000]
    #ccmlo = [ 550, 600, 600, 750, 800]  
    #ccmhi = [ 700, 850, 900,1100,1200]  
    #mjlo  = [  40,  40,  40,  40,  40]  
    #mjhi  = [ 130, 130, 130, 130, 130]  
    #mlo   = [ 400, 400, 400, 600, 600]      
    #mhi   = [1000,1000,1000,1400,1400]          
    #shape    = ["ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","Exp","Exp"]
    #shapeAlt = [  "ErfPow2_v1",  "ErfPow2_v1",  "ErfPow2_v1","Pow","Pow"]
     
    mass  = [ 600, 700, 800, 900,1000]
    ccmlo = [ 550, 600, 700, 750, 800]  
    ccmhi = [ 700, 850, 950,1100,1200]  
    mjlo  = [  40,  40,  40,  40,  40]  
    mjhi  = [ 130, 130, 130, 130, 130]  
    mlo   = [ 400, 400, 600, 600, 600]      
    mhi   = [1000,1000,1400,1400,1400]          
    shape    = ["ErfPowExp_v1","ErfPowExp_v1","Exp","Exp","Exp"]
    shapeAlt = [  "ErfPow2_v1",  "ErfPow2_v1","Pow","Pow","Pow"]

    #mass  = [ 600, 800]
    #ccmlo = [ 550, 700]  
    #ccmhi = [ 700, 950]  
    #mjlo  = [  40,  40]  
    #mjhi  = [ 130, 130]  
    #mlo   = [ 400, 600]      
    #mhi   = [1000,1400]          
    #shape    = ["ErfPowExp_v1","Exp"]
    #shapeAlt = [  "ErfPow2_v1","Pow"]
    
    BRnew = [0,1,2,3,4,5];
    cprime = [1,2,3,4,5,6,7,8,9,10];    
    
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

    if options.makeCards:
        if not os.path.isdir("log"): os.system("mkdir log" );
        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";                
                    print "--------------------------------------------------";                
                    print "R U N N I N G   F I T S" 
                    print "mass = ",mass[i],", cprime = ",cprime[j],", brnew = ",BRnew[k],", channel: ",options.channel;
                    print "--------------------------------------------------";                
                    print "--------------------------------------------------";  
                    
                    time.sleep(0.3);
                    
                    #command = "python doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew %01d --odir %s > %s/log/log_%s_ggH%03d_%02d_%02d_%02d_%02d_%02d_%02d_%s_%s_cprime_%02d_BRnew_%02d "%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k], options.odir, options.odir, CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k]);
                    command = "python doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew %01d --inPath %s"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k], os.getcwd());

                    print command #raw_input("ENTER");
#                    
                    unbinnedCard = options.odir+"/cards_%s/hwwlvj_ggH%03d_%s_%02d_%02d_unbin.txt"%(options.channel,mass[i],options.channel,cprime[j],BRnew[k]);
                    fileExists = os.path.isfile(unbinnedCard)
                    print "fileExists: ",fileExists,", cards: ", unbinnedCard
                    if options.batchMode and not fileExists:
                        fn = "fitScript_%s_%03d_%02d_%02d"%(options.channel,mass[i],cprime[j],BRnew[k]);
                        submitBatchJob( command, fn );
                    if not options.batchMode: 
                        os.system(command);
#                        print "commented out!"

    # =====================================

    if options.computeLimits:

        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";
                    print "--------------------------------------------------";                
                    print "creating card: hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                    combineCmmd = "combineCards.py hwwlvj_ggH%03d_el%s_%02d_%02d_unbin.txt hwwlvj_ggH%03d_mu%s_%02d_%02d_unbin.txt > hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k]);
                    os.system(combineCmmd);

                    
                    print "running card: hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                    runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n hwwlvj_ggH%03d_em%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt %s -v 0"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],mass[i],SIGCH,cprime[j],BRnew[k],moreCombineOpts);                    
                    
                    if options.batchMode:
                        fn = "combineScript_%03d%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                        cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                        submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                    else: 
                        os.system(runCmmd);

    # =====================================
    
    nGraphs = nCprimes*2 + 2;

    tGraphs = [];
    nPoints = len(mass);

    if options.plotLimits:

##        for i in range(mLo,mHi):
##            for j in range(cpLo,cpHi):
        for j in range(cpHi-1,cpLo-1,-1):
            
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

#                if nCprimes > 1:
#                    limitsHist_exp.SetBinContent(i+1,j+1,curAsymLimits[3]);
#                    limitsHist_obs.SetBinContent(i+1,j+1,curAsymLimits[0]);                    
            
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
#
##        MakePlots(tGraphs);
        # ----------------------------
        # Drawing time
    
        banner = TLatex(0.25,0.92,("CMS Preliminary, 19.2-19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
        banner.SetNDC(); banner.SetTextSize(0.028);
        
        leg = ROOT.TLegend(0.15,0.7,0.85,0.9);
        leg.SetFillStyle(0);
        leg.SetBorderSize(0);
        leg.SetNColumns(2);
        for k in range(nCprimes):
            print "k: ",k
            if k == 0: leg.AddEntry(tGraphs[0],"expected, C' = 1","L")
            tmplabel = "observed, C' = %1.1f"%( float((11.-cprime[k])/10.) )
            leg.AddEntry(tGraphs[k*4+1],tmplabel,"L");        
                
                
        drawColors = [];
        drawColors.append(ROOT.kRed+2);
        drawColors.append(ROOT.kGreen+2);
        drawColors.append(ROOT.kBlue+2);
        drawColors.append(ROOT.kMagenta+2);
        drawColors.append(ROOT.kCyan+2);
        drawColors.append(ROOT.kGray+2);                          
        drawColors.append(ROOT.kOrange+2);                                          
        drawColors.append(ROOT.kOrange+2);                                          
        drawColors.append(ROOT.kOrange+2);                                          
        drawColors.append(ROOT.kOrange+2);                                          
                
        oneLine = ROOT.TF1("oneLine","1",599,1001);
        oneLine.SetLineColor(ROOT.kRed);
        oneLine.SetLineWidth(3);

        can0 = ROOT.TCanvas("can0","can0",800,800);
        hrl0 = can0.DrawFrame(599,0.5,1001,100.0);
        hrl0.GetYaxis().SetTitle("#mu' = C' #times #mu");
        hrl0.GetXaxis().SetTitle("mass (GeV)");
        can0.SetGrid();
        tGraphs[3].Draw("f");
        tGraphs[2].Draw("f");
        tGraphs[1].Draw("PL");
        tGraphs[0].Draw("PL");
        oneLine.Draw("LSAMES");
        leg.Draw();
        ROOT.gPad.SetLogy();
        can0.SaveAs("limitFigs/test.eps");
    
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

        for k in range(nCprimes):
#            if k == cpHi-1: 
#                tGraphs[3].Draw("f");
#                tGraphs[2].Draw("f");            
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
        can.SaveAs("limitFigs/test_cPrime.eps");

        # ----------
        
        banner = TLatex(0.25,0.92,("CMS Preliminary, 19.2-19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
        banner.SetNDC(); banner.SetTextSize(0.028);
        
        cprimeHist = [0.05,0.15,0.25,0.35,0.45,0.6,0.8,1.2];
        massHist = [550,650,750,850,950,1050];
        limitsHist_exp = ROOT.TH2F("limitsHist_exp","",nPoints,550,1050,nCprimes,0.05,1.05);
        limitsHist_obs = ROOT.TH2F("limitsHist_obs","",nPoints,550,1050,nCprimes,0.05,1.05);

        limitsHist_exp_600_2D = ROOT.TH2F("limitsHist_exp_600_2D","",nCprimes,0.05,1.05,nBRnews,-0.05,0.55);
        limitsHist_obs_600_2D = ROOT.TH2F("limitsHist_obs_600_2D","",nCprimes,0.05,1.05,nBRnews,-0.05,0.55); 
        
        for j in range(cpLo,cpHi):
            for i in range(mLo,mHi):
                curFile = "higgsCombinehwwlvj_ggH%03d_em_%02d_00_unbin.Asymptotic.mH%03d.root"%(mass[i],cprime[j],mass[i]);
                curAsymLimits = getAsymLimits(curFile);
                limitsHist_exp.SetBinContent(i+1,j+1,curAsymLimits[3]);
                limitsHist_obs.SetBinContent(i+1,j+1,curAsymLimits[0]);          

        for j in range(cpLo,cpHi):
            for i in range(brLo,brHi):
                curFile = "higgsCombinehwwlvj_ggH%03d_em_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(600,cprime[j],BRnew[i],600);
                curAsymLimits = getAsymLimits(curFile);
                limitsHist_exp_600_2D.SetBinContent(j+1,i+1,curAsymLimits[3]);
                limitsHist_obs_600_2D.SetBinContent(j+1,i+1,curAsymLimits[0]);          


        limitsHist_exp.SetTitle("; mass (GeV); C'^{2}; Expected Upper Limit, #mu'");
        limitsHist_obs.SetTitle("; mass (GeV); C'^{2}; Observed Upper Limit, #mu'");
        limitsHist_exp.GetZaxis().RotateTitle(1);
        limitsHist_obs.GetZaxis().RotateTitle(1);
    
        limitsHist_exp_600_2D.SetTitle("; C'^{2}; BR_{New}; Expected Upper Limit, #mu'");
        limitsHist_obs_600_2D.SetTitle("; C'^{2}; BR_{New}; Observed Upper Limit, #mu'");
        limitsHist_exp_600_2D.GetZaxis().RotateTitle(1);
        limitsHist_obs_600_2D.GetZaxis().RotateTitle(1);

        can2d_exp = ROOT.TCanvas("can2d_exp","can2d_exp",1000,800);
        limitsHist_exp.SetStats(0);
        limitsHist_exp.Draw("colz");
        banner.Draw();
        can2d_exp.SaveAs("limitFigs/test2D_exp.eps");

        can2d_obs = ROOT.TCanvas("can2d_obs","can2d_obs",1000,800);
        limitsHist_obs.SetStats(0);
        limitsHist_obs.Draw("colz");
        banner.Draw();
        can2d_obs.SaveAs("limitFigs/test2D_obs.eps");

        can2dbr_exp = ROOT.TCanvas("can2dbr_exp","can2dbr_exp",1000,800);
        limitsHist_exp_600_2D.SetStats(0);
        limitsHist_exp_600_2D.Draw("colz");
        banner.Draw();
        can2dbr_exp.SaveAs("limitFigs/test2Dbr_exp.eps");

        can2dbr_obs = ROOT.TCanvas("can2dbr_obs","can2dbr_obs",1000,800);
        limitsHist_obs_600_2D.SetStats(0);
        limitsHist_obs_600_2D.Draw("colz");
        banner.Draw();
        can2dbr_obs.SaveAs("limitFigs/test2Dbr_obs.eps");



