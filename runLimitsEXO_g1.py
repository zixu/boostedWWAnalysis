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
parser.add_option('--odir',action="store",type="string",dest="odir",default=".")
parser.add_option('--category',action="store",type="string",dest="category",default="HP")
parser.add_option('--closuretest', action='store',type="int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False,
                                    help='no X11 windows')


(options, args) = parser.parse_args()
############################################################

#globals 

mass  = [1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500] 

ccmlo = [ 900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400]  
ccmhi = [1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000]  

mjlo  = [  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40]  
mjhi  = [ 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130]  

mlo   = [ 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800,800]      
mhi   = [3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000]          

shapeAlt = ["ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail"]
shape    = ["ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN"]
#shape    = ["Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp","Exp"]
#shapeAlt = ["Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow","Pow"]
BRnew = 00;
cprime = [10];

#cross-sections for kappa = 0.5
xsDict2 = {600:0.32298, 
          700:0.11827,
          800:0.04931,
          900:0.022506,
          1000:0.011035,
          1100:0.0056883,
          1200:0.0030626,
          1300:0.0017003,
          1400:0.00097456,
          1500:0.00056979,
          1600:0.00034149,
          1700:0.00020677,
          1800:0.000127,
          1900:7.9677e-05,
          2000:5.0345e-05,    
          2100:3.198e-05,    
          2200:2.0502e-05,    
          2300:1.324e-05,    
          2400:8.6099e-06,    
          2500:5.6338e-06}
#cross-sections for kappa = 0.2
xsDict  = {600:0.052087, 
          700:0.019006,
          800:0.0079064,
          900:0.0036364,
          1000:0.0017742,
          1100:0.00091785,
          1200:0.00049262,
          1300:0.00027418,
          1400:0.00015697,
          1500:9.2073e-05,
          1600:5.4715e-05,
          1700:3.3199e-05,
          1800:2.0367e-05,
          1900:1.2723e-05,
          2000:8.0046e-06,    
          2100:5.0566e-06,    
          2200:3.2608e-06,    
          2300:2.0938e-06,    
          2400:1.3566e-06,    
          2500:8.8518e-07}

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

    # create a dummy bash/sh
    outScript=open(fn+".sh","w");

    outScript.write('#!/bin/bash');
    outScript.write("\n"+'date');
    outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
    outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');
    #outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'cd /uscms_data/d3/zixu/BoostJet/CMSSW_6_1_1/src ');
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
    condorScript.write("\n"+'Transfer_Input_Files = g1_exo_doFit_class.py')
    condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log    = out_$(Cluster).log')
    condorScript.write("\n"+'Notification    = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();

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
    outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');
    #outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'cd /uscms_data/d3/zixu/BoostJet/CMSSW_6_1_1/src ');
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd -');
    outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
    outScript.write("\n"+'echo ${PATH}');
    outScript.write("\n"+'ls');
    outScript.write("\n"+command);
    outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *');
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
                                               

def getPValueFromCard( file ):

    f = ROOT.TFile(file);
    t = f.Get("limit");
    entries = t.GetEntries();
    
    lims = 1;
    
    for i in range(1):
        
        t.GetEntry(i);
        lims = t.limit
    
    return lims;

def doULPlot( suffix ):

    
    xbins = array('d', [])
    xbins_env = array('d', [])
    ybins_exp = array('d', [])
    ybins_obs = array('d', [])            
    ybins_1s = array('d', [])                        
    ybins_2s = array('d', [])    
    ybins_xs_02 = array('d', [])    
    ybins_xs_05 = array('d', [])    
    
    for i in range(len(mass)):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        sf = xsDict[mass[i]];
        sf2 = xsDict2[mass[i]];
        #print "curFile: ",curFile
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );                                
        ybins_exp.append( curAsymLimits[3]*sf );                
        ybins_obs.append( curAsymLimits[0]*sf );                                
        ybins_2s.append( curAsymLimits[1]*sf );                                
        ybins_1s.append( curAsymLimits[2]*sf ); 
        ybins_xs_02.append(sf); 
        ybins_xs_05.append(sf2); 
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        sf = xsDict[mass[i]];
        curAsymLimits = getAsymLimits(curFile);   
        xbins_env.append( mass[i] );                                
        ybins_2s.append( curAsymLimits[5]*sf );                                
        ybins_1s.append( curAsymLimits[4]*sf ); 
    
    
    curGraph_exp = ROOT.TGraph(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraph(nPoints,xbins,ybins_obs);
    curGraph_xs_02 = ROOT.TGraph(nPoints,xbins,ybins_xs_02);
    curGraph_xs_05 = ROOT.TGraph(nPoints,xbins,ybins_xs_05);
    curGraph_1s = ROOT.TGraph(nPoints*2,xbins_env,ybins_1s);
    curGraph_2s = ROOT.TGraph(nPoints*2,xbins_env,ybins_2s);
    
    curGraph_obs.SetMarkerStyle(20);
    curGraph_exp.SetLineStyle(2);
    curGraph_xs_02.SetLineStyle(2);
    curGraph_xs_05.SetLineStyle(1);
    curGraph_exp.SetLineWidth(2);
    curGraph_obs.SetLineWidth(2);                    
    curGraph_xs_02.SetLineWidth(3);                    
    curGraph_xs_02.SetLineWidth(3);                    
    curGraph_exp.SetMarkerSize(2);
    curGraph_xs_02.SetMarkerSize(2);
    curGraph_xs_05.SetMarkerSize(2);
    curGraph_obs.SetMarkerSize(2);                    
    curGraph_xs_02.SetLineColor(ROOT.kRed+1);
    curGraph_xs_05.SetLineColor(ROOT.kRed+1);
    curGraph_1s.SetFillColor(ROOT.kGreen);
    curGraph_2s.SetFillColor(ROOT.kYellow);

    # -------
    if suffix =="_el_HP" :        
        banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, ele HP"));
        banner.SetNDC(); banner.SetTextSize(0.028);
    if suffix =="_el_LP" :        
        banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, ele LP"));
        banner.SetNDC(); banner.SetTextSize(0.028);
    if suffix =="_mu_HP" :        
        banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, #mu HP"));
        banner.SetNDC(); banner.SetTextSize(0.028);
    if suffix =="_mu_LP" :        
        banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, #mu LP"));
        banner.SetNDC(); banner.SetTextSize(0.028);
    if suffix =="_em" :        
        banner = TLatex(0.44,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
        banner.SetNDC(); banner.SetTextSize(0.028);
        
    oneLine = ROOT.TF1("oneLine","1",599,1001);
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    can_SM = ROOT.TCanvas("can_SM","can_SM",1000,800);
    hrl_SM = can_SM.DrawFrame(999,1e-3,2501,5.0);
    hrl_SM.GetYaxis().SetTitle("#sigma_{95%} #times BR_{WW} (pb)");
    hrl_SM.GetYaxis().SetTitleOffset(1.4);
    hrl_SM.GetXaxis().SetTitle("mass (GeV/c^{2})");
    hrl_SM.GetYaxis().SetRangeUser(0.0001,5);
    can_SM.SetGrid();

    curGraph_2s.Draw("F");
    curGraph_1s.Draw("F");
    curGraph_obs.Draw("PL");
    curGraph_exp.Draw("PL");
    curGraph_xs_02.Draw("PL");
    curGraph_xs_05.Draw("PL");

    leg2 = ROOT.TLegend(0.55,0.65,0.85,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(0);
    leg2.AddEntry(curGraph_obs,"Observed","L")
    leg2.AddEntry(curGraph_exp,"Expected","L")
    leg2.AddEntry(curGraph_1s,"Expected, #pm 1#sigma","F")
    leg2.AddEntry(curGraph_2s,"Expected, #pm 2#sigma","F")
#    leg2.AddEntry(oneLine,"SM Expected","L")
    leg2.AddEntry(curGraph_xs_02,"#sigma_{Bulk}#times BR(G#rightarrow WW), #tilde{k}=0.2","L")
    leg2.AddEntry(curGraph_xs_05,"#sigma_{Bulk}#times BR(G#rightarrow WW), #tilde{k}=0.5","L")

    leg2.Draw();
    banner.Draw();
#    oneLine.Draw("LESAMES");

    ROOT.gPad.SetLogy();
    can_SM.SaveAs("~/limitFigs/Lim%s.png"%(suffix));                      
    can_SM.SaveAs("~/limitFigs/Lim%s.pdf"%(suffix));       

############################################################

if __name__ == '__main__':


    ###############
    
    CHAN = options.channel;
    DIR = "cards_exo_"+CHAN;


    
    moreCombineOpts = "";
    
    ###############
    
    if options.makeCards:
        if not os.path.isdir(DIR): os.system("mkdir "+DIR);
        else: 
            print "Directory "+DIR+" already exists...";
            #sys.exit(0);

    if options.computeLimits or options.plotLimits: os.chdir("cards_em_EXO_allCat_v2_ExpN");


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
                
                time.sleep(0.3);
                
                #command = "nohup python EXO_doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %02d --BRnew 00  > log/log_%s_ggH%03d_%02d_%02d_%02d_%02d_%02d_%02d_%s_%s_cprime_%02d_BRnew_00 &"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j]);
                #command = "nohup python g1_exo_doFit_class.py %s BulkG_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %02d --BRnew 00  &"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j]);
                command = "python g1_exo_doFit_class.py %s BulkG_WW_lvjj_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], os.getcwd(), options.category,options.closuretest);

                print command #raw_input("ENTER");
                unbinnedCard = options.odir+"/cards_%s/hwwlvj_ggH%03d_%s_%02d_%02d_unbin.txt"%(options.channel,mass[i],options.channel,cprime[j],BRnew);

                fileExists = os.path.isfile(unbinnedCard)
                print "fileExists: ",fileExists,", cards: ", unbinnedCard
                if options.batchMode and not fileExists:
                 fn = "fitScript_%s_%03d_%02d_%02d_%s"%(options.channel,mass[i],cprime[j],BRnew,options.category);
                 submitBatchJob( command, fn );
                if not options.batchMode:
                 os.system(command);


    # =====================================

    if options.computeLimits:

        for i in range(mLo,mHi):
            print "------------------"+str(mass[i])+"-------------------------";                
#            print "creating card: wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_counting.txt"%(mass[i],CHAN);
            combineCmmd = "combineCards.py wwlvj_BulkG_WW_lvjj_c0p2_M%03d_el_10_00_HP_unbin.txt wwlvj_BulkG_WW_lvjj_c0p2_M%03d_mu_10_00_HP_unbin.txt wwlvj_BulkG_WW_lvjj_c0p2_M%03d_el_10_00_LP_unbin.txt wwlvj_BulkG_WW_lvjj_c0p2_M%03d_mu_10_00_LP_unbin.txt > wwlvj_BulkG_WW_lvjj_c0p2_M%03d_em_10_00_unbin.txt"%(mass[i],mass[i],mass[i],mass[i],mass[i]);
            print combineCmmd
            os.system(combineCmmd);
            
                
            # mu HP
            runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt"%(mass[i],mass[i],"mu",mass[i],"mu");
            os.system(runCmmd);
            # mu LP
            runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt"%(mass[i],mass[i],"mu",mass[i],"mu");
            os.system(runCmmd);
            # el HP            
            runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt"%(mass[i],mass[i],"el",mass[i],"el");
            os.system(runCmmd);
            # el LP            
            runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt"%(mass[i],mass[i],"el",mass[i],"el");
            os.system(runCmmd);
            # combined            
            runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt"%(mass[i],mass[i],"em",mass[i],"em");
            os.system(runCmmd);

            #############
            runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping --rMin 1 --rMax 100000 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 0"%(mass[i],mass[i],"el",mass[i],"el");
            os.system(runCmmd2);
            runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping --rMin 1 --rMax 100000 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 0"%(mass[i],mass[i],"mu",mass[i],"mu");
            os.system(runCmmd2);
            runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping --rMin 1 --rMax 100000 -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 0"%(mass[i],mass[i],"el",mass[i],"el");
            os.system(runCmmd2);
            runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping --rMin 1 --rMax 100000 -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 0"%(mass[i],mass[i],"mu",mass[i],"mu");
            os.system(runCmmd2);
            runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping --rMin 1 --rMax 100000 -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt -v 0"%(mass[i],mass[i],"em",mass[i],"em");
            os.system(runCmmd2);

    # =====================================
    
    nGraphs = nCprimes*2 + 2;

    tGraphs = [];
    nPoints = len(mass);
    if options.plotLimits:

        
        ####################################################
        
        xbins = array('d', [])
        ybins_el_hp = array('d', [])
        ybins_mu_hp = array('d', [])
        ybins_el_lp = array('d', [])
        ybins_mu_lp = array('d', [])
        ybins_em = array('d', [])

        for i in range(mLo,mHi):
            print "------------------"+str(mass[i])+"-------------------------";                
            print "getting card: higgsCombine_pval_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
            print "getting card: higgsCombine_pval_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
            print "getting card: higgsCombine_pval_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
            print "getting card: higgsCombine_pval_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
            print "getting card: higgsCombine_pval_%03d_%s.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i])

            orootname_el_hp = "higgsCombine_pval_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
            orootname_mu_hp = "higgsCombine_pval_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
            orootname_el_lp = "higgsCombine_pval_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
            orootname_mu_lp = "higgsCombine_pval_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
            orootname_em = "higgsCombine_pval_%03d_%s.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

            xbins.append( mass[i] );
            ybins_el_hp.append( getPValueFromCard(orootname_el_hp) );
            ybins_mu_hp.append( getPValueFromCard(orootname_mu_hp) );
            ybins_el_lp.append( getPValueFromCard(orootname_el_lp) );
            ybins_mu_lp.append( getPValueFromCard(orootname_mu_lp) );
            ybins_em.append( getPValueFromCard(orootname_em) );

        gr_el_hp = ROOT.TGraph(nPoints,xbins,ybins_el_hp);
        gr_mu_hp = ROOT.TGraph(nPoints,xbins,ybins_mu_hp);
        gr_el_hp.SetLineColor( 1 ); gr_el_hp.SetMarkerColor( 1 ); gr_el_hp.SetMarkerStyle( 20 );  gr_el_hp.SetLineWidth( 3 );gr_el_hp.SetMarkerSize( 1.6 );
        gr_mu_hp.SetLineColor( 1 ); gr_mu_hp.SetMarkerColor( 1 ); gr_mu_hp.SetMarkerStyle( 20 );  gr_mu_hp.SetLineWidth( 3 );gr_mu_hp.SetMarkerSize( 1.6 );
        gr_el_lp = ROOT.TGraph(nPoints,xbins,ybins_el_lp);
        gr_mu_lp = ROOT.TGraph(nPoints,xbins,ybins_mu_lp);
        gr_el_lp.SetLineColor( 1 ); gr_el_lp.SetMarkerColor( 1 ); gr_el_lp.SetMarkerStyle( 20 ); gr_el_lp.SetLineWidth( 3 );gr_el_lp.SetMarkerSize( 1.6 );
        gr_mu_lp.SetLineColor( 1 ); gr_mu_lp.SetMarkerColor( 1 ); gr_mu_lp.SetMarkerStyle( 20 ); gr_mu_lp.SetLineWidth( 3 );gr_mu_lp.SetMarkerSize( 1.6 );
        gr_em = ROOT.TGraph(nPoints,xbins,ybins_em);
        gr_em.SetLineColor( 1 ); gr_em.SetMarkerColor( 1 ); gr_em.SetMarkerStyle( 20 );gr_em.SetLineWidth( 3 ); gr_em.SetMarkerSize( 1.4 );
    
        oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",1000,2500);
        oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
        twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",1000,2500);
        twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
        threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",1000,2500);
        threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
        fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",1000,2500);
        fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);           
    
        banner = TLatex(0.43,0.91,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV"));
        banner.SetNDC(); banner.SetTextSize(0.028);
    
        ban1s = TLatex(2400,1.58655253931457074e-01,("1 #sigma"));
        ban1s.SetTextSize(0.028); ban1s.SetTextColor(2)    
        ban2s = TLatex(2400,2.27501319481792155e-02,("2 #sigma"));
        ban2s.SetTextSize(0.028); ban2s.SetTextColor(2)    
        ban3s = TLatex(2400,1.34989803163009588e-03,("3 #sigma"));
        ban3s.SetTextSize(0.028); ban3s.SetTextColor(2);    
        ban4s = TLatex(2400,3.16712418331199785e-05,("4 #sigma"));
        ban4s.SetTextSize(0.028); ban4s.SetTextColor(2)    
    
        leg2 = ROOT.TLegend(0.2,0.2,0.5,0.35);
        leg2.SetFillStyle(0);
        leg2.SetBorderSize(1);
        leg2.AddEntry( gr_el_hp, "obs signif, e (HP)", "pl" );

        leg3 = ROOT.TLegend(0.2,0.2,0.5,0.35);
        leg3.SetFillStyle(0);
        leg3.SetBorderSize(1);
        leg3.AddEntry( gr_mu_hp, "obs signif, #mu (HP)", "pl" );

        leg4 = ROOT.TLegend(0.2,0.2,0.5,0.35);
        leg4.SetFillStyle(0);
        leg4.SetBorderSize(1);
        leg4.AddEntry( gr_el_lp, "obs signif, e (LP)", "pl" );

        leg5 = ROOT.TLegend(0.2,0.2,0.5,0.35);
        leg5.SetFillStyle(0);
        leg5.SetBorderSize(1);
        leg5.AddEntry( gr_mu_lp, "obs signif, #mu (LP)", "pl" );

        leg6 = ROOT.TLegend(0.2,0.2,0.5,0.35);
        leg6.SetFillStyle(0);
        leg6.SetBorderSize(1);
        leg6.AddEntry( gr_em, "obs signif, all", "pl" );
    
        can = ROOT.TCanvas("can","can",800,800);
        hrl = can.DrawFrame(1000,1e-5,2500,0.6);
        hrl.GetYaxis().SetTitle("p-value");
        hrl.GetXaxis().SetTitle("mass (GeV)");
        can.SetGrid();
        ROOT.gPad.SetLogy();
        gr_el_hp.Draw("PL");
        oneSLine.Draw("same");
        twoSLine.Draw("same");
        threeSLine.Draw("same");
        fourSLine.Draw("same");
        banner.Draw();
        leg2.Draw();
        ban1s.Draw();
        ban2s.Draw();
        ban3s.Draw();
        ban4s.Draw();    
        can.SaveAs("~/limitFigs/pvals_el_HP.pdf","pdf");
        can.SaveAs("~/limitFigs/pvals_el_HP.png","png");

        can2 = ROOT.TCanvas("can2","can2",800,800);
        hrl2 = can2.DrawFrame(1000,1e-5,2500,0.6);
        hrl2.GetYaxis().SetTitle("p-value");
        hrl2.GetXaxis().SetTitle("mass (GeV)");
        can2.SetGrid();
        ROOT.gPad.SetLogy();
        gr_mu_hp.Draw("PL");
        oneSLine.Draw("same");
        twoSLine.Draw("same");
        threeSLine.Draw("same");
        fourSLine.Draw("same");
        banner.Draw();
        leg3.Draw();
        ban1s.Draw();
        ban2s.Draw();
        ban3s.Draw();
        ban4s.Draw();    
        can2.SaveAs("~/limitFigs/pvals_mu_HP.pdf","pdf");
        can2.SaveAs("~/limitFigs/pvals_mu_HP.png","png");

        can3 = ROOT.TCanvas("can3","can3",800,800);
        hrl3 = can3.DrawFrame(1000,1e-5,2500,0.6);
        hrl3.GetYaxis().SetTitle("p-value");
        hrl3.GetXaxis().SetTitle("mass (GeV)");
        can3.SetGrid();
        ROOT.gPad.SetLogy();
        gr_el_lp.Draw("PL");
        oneSLine.Draw("same");
        twoSLine.Draw("same");
        threeSLine.Draw("same");
        fourSLine.Draw("same");
        banner.Draw();
        leg4.Draw();
        ban1s.Draw();
        ban2s.Draw();
        ban3s.Draw();
        ban4s.Draw();    
        can3.SaveAs("~/limitFigs/pvals_el_LP.pdf","pdf");
        can3.SaveAs("~/limitFigs/pvals_el_LP.png","png");

        can4 = ROOT.TCanvas("can4","can4",800,800);
        hrl4 = can4.DrawFrame(1000,1e-5,2500,0.6);
        hrl4.GetYaxis().SetTitle("p-value");
        hrl4.GetXaxis().SetTitle("mass (GeV)");
        can4.SetGrid();
        ROOT.gPad.SetLogy();
        gr_mu_lp.Draw("PL");
        oneSLine.Draw("same");
        twoSLine.Draw("same");
        threeSLine.Draw("same");
        fourSLine.Draw("same");
        banner.Draw();
        leg5.Draw();
        ban1s.Draw();
        ban2s.Draw();
        ban3s.Draw();
        ban4s.Draw();    
        can4.SaveAs("~/limitFigs/pvals_mu_LP.pdf","pdf");
        can4.SaveAs("~/limitFigs/pvals_mu_LP.png","png");

        can5 = ROOT.TCanvas("can5","can5",800,800);
        hrl5 = can5.DrawFrame(1000,1e-5,2500,0.6);
        hrl5.GetYaxis().SetTitle("p-value");
        hrl5.GetXaxis().SetTitle("mass (GeV)");
        can5.SetGrid();
        ROOT.gPad.SetLogy();
        gr_em.Draw("PL");
        oneSLine.Draw("same");
        twoSLine.Draw("same");
        threeSLine.Draw("same");
        fourSLine.Draw("same");
        banner.Draw();
        leg6.Draw();
        ban1s.Draw();
        ban2s.Draw();
        ban3s.Draw();
        ban4s.Draw();    
        can5.SaveAs("~/limitFigs/pvals_em.pdf","pdf");
        can5.SaveAs("~/limitFigs/pvals_em.png","png");



        ####################################################

        # Upper limit plots
        doULPlot("_el_HP");
        doULPlot("_mu_HP");
        doULPlot("_el_LP");
        doULPlot("_mu_LP");
        doULPlot("_em");

        ####################################################





