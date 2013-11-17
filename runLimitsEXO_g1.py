#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

from array import array

import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TH1D
import subprocess
from subprocess import Popen
from optparse import OptionParser

############################################
#             Job steering                 #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

### to make the full analysis fit + datacards
parser.add_option('--makeCards', action='store_true', dest='makeCards', default=False, help='make datacards plus whole analysis')

### to compute limits
parser.add_option('--computeLimits', action='store_true', dest='computeLimits', default=False, help='compute limits')

### to plot limits
parser.add_option('--plotLimits', action='store_true', dest='plotLimits', default=False, help='plot limits')

### to do just signal lineshape fits
parser.add_option('--fitSignal', action='store_true', dest='fitSignal', default=False, help='do signal lineshape fits')

### other options 
parser.add_option('--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('--massPoint',action="store",type="int",dest="massPoint",default=-1)
parser.add_option('--cPrime',action="store",type="int",dest="cPrime",default=-1)
parser.add_option('--odir',action="store",type="string",dest="odir",default=".")
parser.add_option('--category',action="store",type="string",dest="category",default="HP")
parser.add_option('--closuretest', action='store',type="int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False, help='no X11 windows')
parser.add_option('--limitMode', action='store',type="int", dest='limitMode', default=0, help='limit Mode; 0: Asymptotic ; 1: ProfileLikelihood ; 2: FullCLs ; 3: MaxLikelihoodFit')
parser.add_option('--isReReco', action='store',type="int", dest='isReReco', default=0, help='limit Mode; 0: Old signal samples ; 1: New signal Samples')
parser.add_option('--noSys', action='store',type="int", dest='noSys', help='run limit without systematic')
parser.add_option('--plotPvalue', action='store',type="int", default=0, dest='plotPvalue', help='plot p value')
parser.add_option('--signalWidth', action='store',type="int", default=0, dest='signalWidth', help='analysis on non-narrow signals')


(options, args) = parser.parse_args()

#########################################
### Global Variables for running jobs ###
#########################################

### mass point for signal to be fitted
#mass      = [600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500]
mass       = [800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500]
mass_width = [1000,1500,2100]

### mass window for couting analysis
ccmlo = [700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400]
ccmhi = [900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000]

ccmlo_width = [850,1300,1600]
ccmhi_width = [1150,1700,2800]

### jet mass range
mjlo = [ 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]
mjhi = [ 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130]

mjlo_width = [40,40,40]
mjhi_width = [130,130,130]

### mlvj range min and max used when run with option --makeCards
mlo = [ 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700, 700,700]
mhi = [ 3000, 3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000]

mlo_width = [ 700, 700, 700]
mhi_width = [ 3000, 3000,3000]

### mlvj range min and max used when run with option --fitSignal
mlo_sig = [ 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1800, 1800, 1900, 2000]
mhi_sig = [ 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2600, 2700, 2900, 3000]

mlo_sig_width = [ 500, 800, 1000]
mhi_sig_width = [ 1500, 2200, 3200]

### shape to be used for bkg when --makeCards
shapeAlt = ["ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN","ExpN"]
shape    = ["ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail","ExpTail"]

### shape to be used for bkg when --fitSignal

shape_sig_width  = ["BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB","BWDoubleCB"]

#shape_sig_narrow = ["CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1","CB_v1"]

shape_sig_narrow = ["DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1","DoubleCB_v1"]


### signal mass fraction for non narrow samples
mass_fraction = [0.05,0.15,0.3]

BRnew  = [00];
cprime = [10];

##################################################
#  cross-sections for bulk graviton kappa = 0.5  #
##################################################

xsDict2 = {600:0.32298,700:0.11827,800:0.04931,900:0.022506,
          1000:0.011035,1100:0.0056883,1200:0.0030626,1300:0.0017003,
          1400:0.00097456,1500:0.00056979,1600:0.00034149,1700:0.00020677,
          1800:0.000127,1900:7.9677e-05,2000:5.0345e-05,2100:3.198e-05,2200:2.0502e-05,
          2300:1.324e-05,2400:8.6099e-06,2500:5.6338e-06}

##################################################
#  cross-sections for bulk graviton kappa = 0.2  #
##################################################

xsDict = {600:0.052087,700:0.019006,800:0.0079064,
          900:0.0036364,1000:0.0017742,1100:0.00091785,
          1200:0.00049262,1300:0.00027418,1400:0.00015697,
          1500:9.2073e-05,1600:5.4715e-05,1700:3.3199e-05,
          1800:2.0367e-05,1900:1.2723e-05,2000:8.0046e-06,
          2100:5.0566e-06,2200:3.2608e-06,2300:2.0938e-06,
          2400:1.3566e-06,2500:8.8518e-07}


###########################################################
### submit fit jobs on condor batch -> makeCards option ###
###########################################################

def submitBatchJob( command, fn ):

    currentDir = os.getcwd();

    # create a dummy bash/sh -> give to condor the instruction
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
    outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *');  ## pruce a tar file as output
    outScript.close();

    #### link a condor script to your shell script
    condorScript=open("condor_"+fn,"w");
    condorScript.write('universe = vanilla')
    condorScript.write("\n"+"Executable = "+fn+".sh")
    condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
    condorScript.write("\n"+'Should_Transfer_Files = YES')
    condorScript.write("\n"+'Transfer_Input_Files = g1_exo_doFit_class.py')
    condorScript.write("\n"+'WhenToTransferOutput = ON_EXIT_OR_EVICT')
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log = out_$(Cluster).log')
    condorScript.write("\n"+'Notification = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();

    os.system("condor_submit "+"condor_"+fn)


##################################################################
### submit combine jobs on condor batch -> computeLimit option ###
##################################################################


def submitBatchJobCombine( command, fn, channel, mass, cprime, BRnew, purity, isReReco ):

    currentDir = os.getcwd();
    #### create a dummy bash/csh
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
    if not purity=="combo" :
     outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *'+str(mass)+'*'+str(channel)+'*'+str(purity)+'*');
    else:
     outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *'+str(mass)+'*'+str(purity)+'*');
        
    outScript.close();

    if not purity=="combo" :
        if isReReco == 0:
         file1 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_%s_unbin.txt"%(mass,channel,purity);
         file2 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_%s_workspace.root"%(mass,channel,purity);
        else: 
         file1 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_%s_unbin.txt"%(mass,channel,purity);
         file2 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_%s_workspace.root"%(mass,channel,purity);
    else:
       if isReReco == 0: 
          file1 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt"%(mass,purity);
          file2 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_mu_10_00_HP_workspace.root"%(mass);
          file3 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_el_10_00_HP_workspace.root"%(mass);
          file4 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_mu_10_00_LP_workspace.root"%(mass);
          file5 = "wwlvj_BulkG_WW_lvjj_c0p2_M%03d_el_10_00_LP_workspace.root"%(mass);
       else:   
          file1 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt"%(mass,purity);
          file2 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_mu_10_00_HP_workspace.root"%(mass);
          file3 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_el_10_00_HP_workspace.root"%(mass);
          file4 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_mu_10_00_LP_workspace.root"%(mass);
          file5 = "wwlvj_BulkG_WW_inclusive_c0p2_M%03d_el_10_00_LP_workspace.root"%(mass);
        
    # link a condor script to your shell script
    condorScript=open("subCondor_"+fn,"w");
    condorScript.write('universe = vanilla')
    condorScript.write("\n"+"Executable = "+fn+".sh")
    condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
    condorScript.write("\n"+'Should_Transfer_Files = YES')
    if not purity=="combo" : condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2)
    else : condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2+', '+file3+', '+file4+', '+file5)
    condorScript.write("\n"+'WhenToTransferOutput = ON_EXIT_OR_EVICT')
                                                                    
    condorScript.write("\n"+'Output = out_$(Cluster).stdout')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Error = out_$(Cluster).stderr')
    condorScript.write("\n"+'Log = out_$(Cluster).log')
    condorScript.write("\n"+'Notification = Error')
    condorScript.write("\n"+'Queue 1')
    condorScript.close();


    os.system("condor_submit "+"subCondor_"+fn)



##################################################
### Get Limit Value from Combine -M Asymptotic ###
##################################################

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
        elif t_quantileExpected >= 0.15  and t_quantileExpected <= 0.17:  lims[2] = t_limit;
        elif t_quantileExpected == 0.5: lims[3] = t_limit;
        elif t_quantileExpected >= 0.83  and t_quantileExpected <= 0.85:  lims[4] = t_limit;
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit;
        else: print "Unknown quantile!"
    
    return lims;


####################################################
### Get PValue from combine -M ProfileLikelihood ###
####################################################

def getPValueFromCard( file ):

    f = ROOT.TFile(file);
    t = f.Get("limit");
    entries = t.GetEntries();
    
    lims = 1;
    
    for i in range(1):
        
        t.GetEntry(i);
        lims = t.limit
    
    return lims;


##########################################
### Make Limit Plot --> Brazilian Plot ###
##########################################

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
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3]*sf2 );
        ybins_obs.append( curAsymLimits[0]*sf2 );
        ybins_2s.append( curAsymLimits[1]*sf2 );
        ybins_1s.append( curAsymLimits[2]*sf2 );
        ybins_xs_02.append(sf);
        ybins_xs_05.append(sf2);
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        sf = xsDict[mass[i]];
        sf2 = xsDict2[mass[i]];
        curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5]*sf2 );
        ybins_1s.append( curAsymLimits[4]*sf2 );
    
    
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

    if suffix =="_el_HP" :
        banner = TLatex(0.31,0.91,("CMS Preliminary, 19.5 fb^{-1} at #sqrt{s}=8TeV, W#rightarrow e#nu HP"));
        banner.SetNDC(); banner.SetTextSize(0.032);
    if suffix =="_el_LP" :
        banner = TLatex(0.31,0.91,("CMS Preliminary, 19.5 fb^{-1} at #sqrt{s}=8TeV, W#rightarrow e#nu LP"));
        banner.SetNDC(); banner.SetTextSize(0.032);
    if suffix =="_mu_HP" :
        banner = TLatex(0.31,0.91,("CMS Preliminary, 19.5 fb^{-1} at #sqrt{s}=8TeV, W#rightarrow #mu#nu HP"));
        banner.SetNDC(); banner.SetTextSize(0.032);
    if suffix =="_mu_LP" :
        banner = TLatex(0.31,0.91,("CMS Preliminary, 19.5 fb^{-1} at #sqrt{s}=8TeV, W#rightarrow #mu#nu LP"));
        banner.SetNDC(); banner.SetTextSize(0.032); 
    if suffix =="_combo" :
        banner = TLatex(0.31,0.91,("CMS Preliminary, 19.5 fb^{-1} at #sqrt{s}=8TeV, e+#mu combined"));
        banner.SetNDC(); banner.SetTextSize(0.032);
    if suffix =="_em_HP" :
        banner = TLatex(0.31,0.91,("CMS Preliminary, 19.5 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
        banner.SetNDC(); banner.SetTextSize(0.032);
        
    oneLine = ROOT.TF1("oneLine","1",799,1001);
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    can_SM = ROOT.TCanvas("can_SM","can_SM",1000,800);
    hrl_SM = can_SM.DrawFrame(799,1e-5,2501,5.0);
    hrl_SM.GetYaxis().SetTitle("#sigma_{95%} #times BR_{WW} (pb)");
    hrl_SM.GetYaxis().SetTitleOffset(1.2);
    hrl_SM.GetXaxis().SetTitle("m_{WW} (GeV)");
    hrl_SM.GetXaxis().SetTitleSize(0.036);
    hrl_SM.GetXaxis().SetTitleSize(0.036);
    hrl_SM.GetXaxis().SetLabelSize(0.036);
    hrl_SM.GetYaxis().SetRangeUser(0.0001,10);
    hrl_SM.GetYaxis().SetNdivisions(504);
    can_SM.SetGrid();

    curGraph_2s.Draw("F");
    curGraph_1s.Draw("F");
    curGraph_obs.Draw("PL");
    curGraph_exp.Draw("PL");
    curGraph_xs_02.Draw("PL");
    curGraph_xs_05.Draw("PL");

       
    leg2 = ROOT.TLegend(0.52,0.60,0.88,0.85);
    leg2.SetFillColor(0);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(3001);
    leg2.AddEntry(curGraph_obs,"Observed","LP")
    leg2.AddEntry(curGraph_exp,"Expected","L")
    leg2.AddEntry(curGraph_1s,"Expected, #pm 1#sigma","F")
    leg2.AddEntry(curGraph_2s,"Expected, #pm 2#sigma","F")
    leg2.AddEntry(curGraph_xs_02,"#sigma_{Bulk}#times BR(G#rightarrow WW), #tilde{k}=0.2","L")
    leg2.AddEntry(curGraph_xs_05,"#sigma_{Bulk}#times BR(G#rightarrow WW), #tilde{k}=0.5","L")

    ROOT.gPad.SetLogy();
    can_SM.Update();
    can_SM.RedrawAxis();
    can_SM.RedrawAxis("g");
    can_SM.Update();

    leg2.Draw();
    banner.Draw();

    can_SM.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/Lim%s.png"%(suffix));
    can_SM.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/Lim%s.pdf"%(suffix));
    can_SM.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/Lim%s.root"%(suffix));
    can_SM.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/Lim%s.C"%(suffix));


#################
### Main Code ###    
#################
    
if __name__ == '__main__':

    
    CHAN = options.channel;
    
    moreCombineOpts = "";

    ### Set the working directory
        
    if (options.computeLimits or options.plotLimits) and options.limitMode == 2 : os.chdir("cards_em_EXO_allCat_v2_ExpTail_g1_rereco_c0p5_modelIndependent");
    elif (options.computeLimits or options.plotLimits) : os.chdir("cards_em_EXO_allCat_v2_ExpTail_g1_rereco_c0p5_modelIndependent");

    
    ### put in functionality to test just one mass point or just one cprime

    nMasses = len(mass);
    mLo = 0;
    mHi = nMasses;

    nMasses_width = len(mass_width);
    mLo_width = 0;
    mHi_width = nMasses_width;

    nCprimes = len(cprime);
    cpLo = 0;
    cpHi = nCprimes;

    if options.massPoint > 0 and not options.signalWidth:
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;

    elif options.massPoint > 0 and options.signalWidth:
        curIndex = mass_width.index(options.massPoint);
        mLo_width = curIndex;
        mHi_width = curIndex+1;
        

    if options.cPrime > 0:
        curIndex = cprime.index(options.cPrime);
        cpLo = curIndex;
        cpHi = curIndex+1;

    ### Make cards option analysis
    if options.makeCards:
     if options.signalWidth == 0:
        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                
                print "##################################################";
                print "##################################################";
                print "############# R U N N I N G F I T S ##############";
                print "mass = ",mass[i],", cprime = ",cprime[j];
                print "###################################################";
                print "###################################################";
                
                time.sleep(0.3);
                if options.isReReco == 0 :
                 command = "python g1_exo_doFit_class.py %s BulkG_WW_lvjj_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], os.getcwd(), options.category,options.closuretest);
                else : 
                 command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], os.getcwd(), options.category,options.closuretest);


                print command ;

                if options.batchMode :
                 fn = "fitScript_%s_%03d_%02d_%02d_%s"%(options.channel,mass[i],cprime[j],BRnew[0],options.category);
                 submitBatchJob(command,fn);

                if not options.batchMode:
                 os.system(command);

     elif options.signalWidth == 1:
      for k in range(len(mass_fraction)):
        for i in range(mLo_width,mHi_width):
            for j in range(cpLo,cpHi):
                
                print "###################################################";
                print "###################################################";
                print "####### R U N N I N G F I T S   W I D T H #########";
                print "mass = ",mass_width[i],", cprime = ",cprime[j]," mass_fraction = ",mass_fraction[k];
                print "##################################################";
                print "##################################################";
                
                time.sleep(0.3);
                if options.isReReco == 1 :
                 if int(mass_width[i]*mass_fraction[k]) < 100:
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%02d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d"%(CHAN, mass_width[i], mass_width[i]*mass_fraction[k],ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_width[i], mhi_width[i], shape[i], shapeAlt[i], cprime[j], os.getcwd(), options.category,options.closuretest);
                 elif int(mass_width[i]*mass_fraction[k]) >= 100 and int(mass_width[i]*mass_fraction[k]) < 1000:
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W03%d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d"%(CHAN, mass_width[i], mass_width[i]*mass_fraction[k],ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_width[i], mhi_width[i], shape[i], shapeAlt[i], cprime[j], os.getcwd(), options.category,options.closuretest);
                 else:
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%04d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d"%(CHAN, mass_width[i], mass_width[i]*mass_fraction[k],ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_width[i], mhi_width[i], shape[i], shapeAlt[i], cprime[j], os.getcwd(), options.category,options.closuretest);

                print command ;

                if options.batchMode :
                 fn = "fitScript_%s_W%s_%03d_%02d_%02d_%s"%(options.channel,mass_width[i],int(mass_width[i]*mass_fraction[k]),cprime[j],BRnew[0],options.category);
                 submitBatchJob(command,fn);
                if not options.batchMode:
                 os.system(command);



    ### Make signal fit option
    if options.fitSignal:
     if options.signalWidth == 0:
        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                
                print "##############################################";
                print "##############################################";
                print "############# FIT SIGNAL SHAPES ##############";
                print "mass = ",mass[i],", cprime = ",cprime[j];
                print "##############################################";
                print "##############################################";
                
                time.sleep(0.3);
                if options.isReReco == 0 :
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_lvjj_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo_sig[i], mhi_sig[i], shape_sig_narrow[i], shape_sig_width[i], cprime[j], os.getcwd(), options.category,options.closuretest);
                elif options.isReReco == 1 : 
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_c0p2_M%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo_sig[i], mhi_sig[i], shape_sig_narrow[i], shape_sig_width[i], cprime[j], os.getcwd(), options.category,options.closuretest);
          
                print command ;

                if options.batchMode :
                 fn = "fitScript_signal_%s_%03d_%02d_%02d_%s"%(options.channel,mass[i],cprime[j],BRnew[0],options.category);
                 submitBatchJob(command,fn);

                if not options.batchMode:
                 os.system(command);


     elif options.signalWidth == 1:
      for k in range(len(mass_fraction)):
        for i in range(mLo_width,mHi_width):
            for j in range(cpLo,cpHi):
                
                print "###################################################";
                print "###################################################";
                print "####### R U N N I N G F I T S   W I D T H #########";
                print "mass = ",mass_width[i],", cprime = ",cprime[j]," mass_fraction = ",mass_fraction[k];
                print "##################################################";
                print "##################################################";
                
                time.sleep(0.3);
                if options.isReReco == 1 :
                 if int(mass_width[i]*mass_fraction[k]) < 100:
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%02d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass_width[i], mass_width[i]*mass_fraction[k],ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_sig_width[i], mhi_sig_width[i], shape_sig_narrow[i], shape_sig_width[i], cprime[j], os.getcwd(), options.category,options.closuretest);
                 elif int(mass_width[i]*mass_fraction[k]) >= 100 and int(mass_width[i]*mass_fraction[k]) <= 1000:
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass_width[i], mass_width[i]*mass_fraction[k],ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_sig_width[i], mhi_sig_width[i], shape_sig_narrow[i], shape_sig_width[i], cprime[j], os.getcwd(), options.category,options.closuretest);
                 else :
                  command = "python g1_exo_doFit_class.py %s BulkG_WW_inclusive_M%03d_W%04d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew 00 --inPath %s --category %s --closuretest %01d --fitSignal 1"%(CHAN, mass_width[i], mass_width[i]*mass_fraction[k],ccmlo_width[i], ccmhi_width[i], mjlo_width[i], mjhi_width[i], mlo_sig_width[i], mhi_sig_width[i], shape_sig_narrow[i], shape_sig_width[i], cprime[j], os.getcwd(), options.category,options.closuretest);

                print command ;

                if options.batchMode :
                 fn = "fitScript_%s_W%s_%03d_%02d_%02d_%s"%(options.channel,mass_width[i],int(mass_width[i]*mass_fraction[k]),cprime[j],BRnew[0],options.category);
                 submitBatchJob(command,fn);
                if not options.batchMode:
                 os.system(command);


                  
    ### Compute Limits
    if options.computeLimits:

        for i in range(mLo,mHi):
            print "##################"+str(mass[i])+"#####################";

            time.sleep(0.3);

            ### combine different categories HP and LP
            if options.channel != "em" :
             if options.isReReco == 0 :
              combineCmmd = "combineCards.py wwlvj_BulkG_WW_lvjj_c0p2_M%03d_el_10_00_HP_unbin.txt wwlvj_BulkG_WW_lvjj_c0p2_M%03d_mu_10_00_HP_unbin.txt wwlvj_BulkG_WW_lvjj_c0p2_M%03d_el_10_00_LP_unbin.txt wwlvj_BulkG_WW_lvjj_c0p2_M%03d_mu_10_00_LP_unbin.txt > wwlvj_BulkG_WW_lvjj_c0p2_M%03d_combo_10_00_unbin.txt"%(mass[i],mass[i],mass[i],mass[i],mass[i]);
             elif options.isReReco == 1 and options.channel != "em" : 
              combineCmmd = "combineCards.py wwlvj_BulkG_WW_inclusive_c0p2_M%03d_el_10_00_HP_unbin.txt wwlvj_BulkG_WW_inclusive_c0p2_M%03d_mu_10_00_HP_unbin.txt wwlvj_BulkG_WW_inclusive_c0p2_M%03d_el_10_00_LP_unbin.txt wwlvj_BulkG_WW_inclusive_c0p2_M%03d_mu_10_00_LP_unbin.txt > wwlvj_BulkG_WW_inclusive_c0p2_M%03d_combo_10_00_unbin.txt"%(mass[i],mass[i],mass[i],mass[i],mass[i]);

             print combineCmmd
             os.system(combineCmmd);
            
            ### Asymptotic Limit + profileLikelihood to have an hint of the r value
            if options.limitMode==0:

             ## ele HP
             if options.channel != "em":   
              if options.isReReco == 0 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");
              elif options.isReReco == 0 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");
              elif options.isReReco == 1 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");
              elif options.isReReco == 1 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");

              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd2,fn,"el",mass[i],cprime[0],BRnew[0],"HP",options.isReReco);
              else:  
               os.system(runCmmd2);

              time.sleep(0.1);

              ## mu HP
              if options.isReReco == 0 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              elif options.isReReco == 0 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              elif options.isReReco == 1 and not  options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              elif options.isReReco == 1 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");

              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd2,fn,"mu",mass[i],cprime[0],BRnew[0],"HP",options.isReReco);
              else:  
               os.system(runCmmd2);

              time.sleep(0.1);

              ## ele LP
              if options.isReReco == 0 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");
              elif options.isReReco == 0 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");
              elif options.isReReco == 1 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");
              elif options.isReReco == 1 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"el",mass[i],"el");

              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine(runCmmd2,fn,"el",mass[i],cprime[0],BRnew[0],"LP",options.isReReco);
              else:  
               os.system(runCmmd2);

              time.sleep(0.1);

              ## mu LP
              if options.isReReco == 0 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              elif options.isReReco == 0 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");              
              elif options.isReReco == 1 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");
              elif options.isReReco == 1 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt -v 2"%(mass[i],mass[i],"mu",mass[i],"mu");

              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine(runCmmd2,fn,"mu",mass[i],cprime[0],BRnew[0],"LP",options.isReReco);
              else:  
               os.system(runCmmd2);

              time.sleep(0.1);

              ## combined
              if options.isReReco == 0 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt -v 2"%(mass[i],mass[i],"combo",mass[i],"combo");
              elif options.isReReco == 0 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt -v 2"%(mass[i],mass[i],"combo",mass[i],"combo");
              elif options.isReReco == 1 and not  options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt -v 2"%(mass[i],mass[i],"combo",mass[i],"combo");
              elif options.isReReco == 1 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt -v 2"%(mass[i],mass[i],"combo",mass[i],"combo");
              
              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%02d_%02d"%("combo",mass[i],"combo",cprime[0],BRnew[0]);
               submitBatchJobCombine(runCmmd2,fn,"",mass[i],cprime[0],BRnew[0],"combo",options.isReReco);
              else:  
               os.system(runCmmd2);

             elif options.channel == "em" :
              ## combined
              if options.isReReco == 0 and not options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"em",mass[i],"em");
              elif options.isReReco == 0 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"em",mass[i],"em");
              elif options.isReReco == 1 and not  options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"em",mass[i],"em");
              elif options.isReReco == 1 and options.noSys:
               runCmmd2 = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -H ProfileLikelihood -S 0 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -v 2"%(mass[i],mass[i],"em",mass[i],"em");
              
              if options.batchMode:
               fn = "combineScript_Asymptotic_%s_%03d%s_%02d_%02d_%s"%("em",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd2,fn,"em",mass[i],cprime[0],BRnew[0],"HP",options.isReReco);
              else:  
               os.system(runCmmd2);

              time.sleep(0.1);

             ### pvalue evaluation
            if options.limitMode == 1 :
             if options.channel != "em":
                
              ## mu HP
              print "##################### mu HP #####################";   
              if options.isReReco == 0:   
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -t -1 --expectSignal=1 --toysFreq"%(mass[i],mass[i],"mu",mass[i],"mu");
              else: 
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -t -1 --expectSignal=1 --toysFreq"%(mass[i],mass[i],"mu",mass[i],"mu");
              
              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine( runCmmd, fn, "mu",mass[i], cprime[0], BRnew[0],"HP",options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

              ## mu LP
              print "##################### mu LP #####################";   
              runCmmd ="";
             
              if options.isReReco == 0:
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"mu",mass[i],"mu");
              else: 
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"mu",mass[i],"mu");

              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine( runCmmd, fn, "mu",mass[i], cprime[0], BRnew[0],"LP", options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

              ## el HP
              print "##################### el HP #####################";   
              runCmmd ="";
              if options.isReReco == 0:
               runCmmd = "combine -M ProfileLikelihood  --signif --pvalue -m %03d -n _pval_obs_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(mass[i],mass[i],"el",mass[i],"el");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"el",mass[i],"el");
              else: 
               runCmmd = "combine -M ProfileLikelihood  --signif --pvalue -m %03d -n _pval_obs_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(mass[i],mass[i],"el",mass[i],"el");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"el",mass[i],"el");
              
              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd,fn,"el",mass[i],cprime[0],BRnew[0],"HP",options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

              ## el LP
              print "##################### el LP #####################";   
              runCmmd ="";
              if options.isReReco == 0 :
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(mass[i],mass[i],"el",mass[i],"el");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"el",mass[i],"el");
              else: 
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(mass[i],mass[i],"el",mass[i],"el");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"el",mass[i],"el");
             
              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine(runCmmd,fn,"el",mass[i],cprime[0],BRnew[0],"LP",options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

              ## combined
              runCmmd ="";
              print "##################### combo #####################";   
              if options.isReReco == 0 :
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt\n"%(mass[i],mass[i],"combo",mass[i],"combo");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"combo",mass[i],"combo");
              else: 
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt\n"%(mass[i],mass[i],"combo",mass[i],"combo");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"combo",mass[i],"combo");
              
              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%02d_%02d"%("combo",mass[i],"",cprime[0],BRnew[0]);
               submitBatchJobCombine(runCmmd,fn,"combo",mass[i],cprime[0], BRnew[0],"combo",options.isReReco);
              else:  
               os.system(runCmmd);

             elif options.channel == "em":

              ## semplified for el+mu
              runCmmd ="";
              if options.isReReco == 0 :
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(mass[i],mass[i],"em",mass[i],"em");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"em",mass[i],"em");
              else: 
               runCmmd = "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_obs_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(mass[i],mass[i],"em",mass[i],"em");
               runCmmd += "combine -M ProfileLikelihood --signif --pvalue -m %03d -n _pval_exp_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt --expectSignal=1 --toysFreq -t -1"%(mass[i],mass[i],"em",mass[i],"em");
              
              if options.batchMode:
               fn = "combineScript_ProfileLikelihood_%s_%03d%s_%02d_%02d_%s"%("em",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd,fn,"em",mass[i],cprime[0], BRnew[0],"HP",options.isReReco);
              else:  
               os.system(runCmmd);


            ### Full CLs -> naive setup
            if options.limitMode ==2:

             if options.channel != "em":

              ## el HP
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2  -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4  -v 0 -s -1\n"%(mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4  -v 0 -s -1 --expectedFromGrid 0.5\n"%(mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4  -v 0 -s -1 --expectedFromGrid 0.16\n"%(mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4  -v 0 -s -1 --expectedFromGrid 0.84\n"%(mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4  -v 0 -s -1 --expectedFromGrid 0.025\n"%(mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4  -v 0 -s -1 --expectedFromGrid 0.975\n"%(mass[i],mass[i],"el",mass[i],"el");
                                                                     
              if options.batchMode:
                fn = "combineScript_HybridNew_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"HP");
                submitBatchJobCombine( runCmmd3,fn, "el",mass[i], cprime[0], BRnew[0],"HP" );
              else:  
                os.system(runCmmd3);

              time.sleep(0.1);

              ## mu HP
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");

                                                                     
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine( runCmmd3,fn, "mu",mass[i], cprime[0], BRnew[0],"HP" );
              else:  
               os.system(runCmmd3);

              time.sleep(0.1);

              ## el LP              
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"el",mass[i],"el");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"el",mass[i],"el");


              time.sleep(0.1);
              

              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine( runCmmd3,fn, "el",mass[i], cprime[0], BRnew[0],"LP" );
              else:  
               os.system(runCmmd3);

              ## mu LP              
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_LP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"mu",mass[i],"mu");

                                                                     
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine( runCmmd3,fn, "mu",mass[i], cprime[0], BRnew[0],"LP" );
              else:  
               os.system(runCmmd3);

              time.sleep(0.1);

              ## combined
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"combo",mass[i],"combo");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"combo",mass[i],"combo");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"combo",mass[i],"combo");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"combo",mass[i],"combo");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"combo",mass[i],"combo");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"combo",mass[i],"combo");

                                                                     
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%02d_%02d"%("combo",mass[i],"",cprime[0],BRnew[0]);
               submitBatchJobCombine( runCmmd3,fn, "",mass[i], cprime[0], BRnew[0],"combo" );
              else:  
               os.system(runCmmd3);


             elif options.channel == "em":

              ## combined
              runCmmd3="";
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 \n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.5\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.16\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.84\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.025\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");
              runCmmd3 = runCmmd3 + "combine -M HybridNew --frequentist --minimizerAlgo Minuit2 --rMin %s --rMax %s -m %03d -n _lim_%03d_%s_HP -d wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt --testStat LHC -H ProfileLikelihood --fork 4 -v 0 -s -1 --expectedFromGrid 0.975\n"%(tmp_rMin,tmp_rMax, mass[i],mass[i],"em",mass[i],"em");

                                                                     
              if options.batchMode:
               fn = "combineScript_HybridNew_%s_%03d%s_%02d_%02d"%("em",mass[i],"",cprime[0],BRnew[0]);
               submitBatchJobCombine(runCmmd3,fn, "",mass[i], cprime[0], BRnew[0],"em");
              else:  
               os.system(runCmmd3);

            ### maximum likelihood  fit plus toys for S+B -> get bias estimation
            if options.limitMode == 3 :

             tmp_rMin=0.01; tmp_rMax=1000;
             if mass[i]>=2000: tmp_rMin=1000; tmp_rMax=50000;
             if mass[i]>=1500 and mass[i]<2000: tmp_rMin=100; tmp_rMax=10000;
             if mass[i]>=700 and mass[i]<1500:  tmp_rMin=10; tmp_rMax=5000;
             if mass[i]<700: tmp_rMin=0.001; tmp_rMax=1000;

             if options.channel != "em":
                
              ## mu HP
              if options.isReReco == 0:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2  --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -t 500 --expectSignal=600"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
              else :
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2  --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -t 500 --expectSignal=600"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
                 
              print runCmmd
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd,fn,"mu",mass[i],cprime[0],BRnew[0],"HP",options.isReReco);
              else:  
                os.system(runCmmd);

              time.sleep(0.1);

              ## mu LP
              if options.isReReco == 0:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
              else:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"mu",mass[i],"mu");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%02d_%02d_%s"%("mu",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine(runCmmd,fn,"mu",mass[i],cprime[0],BRnew[0],"LP",options.isReReco);
              else:  
               os.system(runCmmd);

              time.sleep(0.1);

              ## el HP
              if options.isReReco == 0:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_HP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_HP_unbin.txt -t 500 --expectSignal=600 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
              else:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_HP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_HP_unbin.txt -t 500 --expectSignal=600 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"HP");
               submitBatchJobCombine(runCmmd,fn,"el",mass[i],cprime[0],BRnew[0],"HP",options.isReReco);
              else:  
                os.system(runCmmd);

              time.sleep(0.1);

              ## el LP
              if options.isReReco == 0:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_LP wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
              else:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s_LP wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_LP_unbin.txt -t 500 --expectSignal=2000 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"el",mass[i],"el");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%02d_%02d_%s"%("el",mass[i],"",cprime[0],BRnew[0],"LP");
               submitBatchJobCombine(runCmmd,fn,"el",mass[i],cprime[0],BRnew[0],"LP",options.isReReco);
              else :  
               os.system(runCmmd);

              time.sleep(0.1);

              ## combined
              if options.isReReco == 0:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"combo",mass[i],"combo");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"combo",mass[i],"combo");
              else:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"combo",mass[i],"combo");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"combo",mass[i],"combo");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%02d_%02d"%("combo",mass[i],"",cprime[0],BRnew[0]);
               submitBatchJobCombine(runCmmd,fn,"combo",mass[i],cprime[0],BRnew[0],"combo",options.isReReco);
              else:  
               os.system(runCmmd);


             elif options.channel == "em":

              ## combined
              if options.isReReco == 0:
               #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"em",mass[i],"em");
               runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_BulkG_WW_lvjj_c0p2_M%03d_%s_10_00_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"em",mass[i],"em");
             else:
              #runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_%03d_%s wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt\n"%(tmp_rMin,tmp_rMax,mass[i],mass[i],"em",mass[i],"em");
              runCmmd = "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %s  --rMax %s -m %03d -n _fit_toy_%03d_%s wwlvj_BulkG_WW_inclusive_c0p2_M%03d_%s_10_00_unbin.txt -t 500 --expectSignal=350 -s -1 "%(tmp_rMin,tmp_rMax,mass[i],mass[i],"em",mass[i],"em");
                 
              if options.batchMode:
               fn = "combineScript_MaxLikelihoodFit_%s_%03d%s_%02d_%02d"%("em",mass[i],"",cprime[0],BRnew[0]);
               submitBatchJobCombine(runCmmd,fn,"em",mass[i],cprime[0],BRnew[0],"em",options.isReReco);
              else:  
               os.system(runCmmd);


                 
    ### make the plots    
    if options.plotLimits:

        nGraphs = nCprimes*2 + 2;

        tGraphs = [];
        nPoints = len(mass);
        
        xbins = array('d', [])
        ybins_el_hp = array('d', [])
        ybins_mu_hp = array('d', [])
        ybins_el_lp = array('d', [])
        ybins_mu_lp = array('d', [])
        ybins_em_hp = array('d', [])
        ybins_combo = array('d', [])

        xbins_exp       = array('d', [])
        ybins_el_hp_exp = array('d', [])
        ybins_mu_hp_exp = array('d', [])
        ybins_el_lp_exp = array('d', [])
        ybins_mu_lp_exp = array('d', [])
        ybins_em_hp_exp = array('d', [])
        ybins_combo_exp = array('d', [])

        if options.channel != "em":

         for i in range(mLo,mHi):
             print "############################"+str(mass[i])+"#############################";            
             print "getting card: higgsCombine_pval_obs_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             print "getting card: higgsCombine_pval_obs_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             print "getting card: higgsCombine_pval_obs_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             print "getting card: higgsCombine_pval_obs_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             print "getting card: higgsCombine_pval_obs_%03d_%s.ProfileLikelihood.mH%03d.root"%(mass[i],"combo",mass[i]);

             orootname_el_hp = "higgsCombine_pval_obs_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             orootname_mu_hp = "higgsCombine_pval_obs_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             orootname_el_lp = "higgsCombine_pval_obs_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             orootname_mu_lp = "higgsCombine_pval_obs_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             orootname_combo = "higgsCombine_pval_obs_%03d_%s.ProfileLikelihood.mH%03d.root"%(mass[i],"combo",mass[i]);

             print "getting card: higgsCombine_pval_exp_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             print "getting card: higgsCombine_pval_exp_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             print "getting card: higgsCombine_pval_exp_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             print "getting card: higgsCombine_pval_exp_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             print "getting card: higgsCombine_pval_exp_%03d_%s.ProfileLikelihood.mH%03d.root"%(mass[i],"combo",mass[i]);

             orootname_el_hp_exp = "higgsCombine_pval_exp_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             orootname_mu_hp_exp = "higgsCombine_pval_exp_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             orootname_el_lp_exp = "higgsCombine_pval_exp_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"el",mass[i]);
             orootname_mu_lp_exp = "higgsCombine_pval_exp_%03d_%s_LP.ProfileLikelihood.mH%03d.root"%(mass[i],"mu",mass[i]);
             orootname_combo_exp = "higgsCombine_pval_exp_%03d_%s.ProfileLikelihood.mH%03d.root"%(mass[i],"combo",mass[i]); 

             xbins.append( mass[i] );
             xbins_exp.append( mass[i] );

             if options.plotPvalue :

              ybins_el_hp.append( getPValueFromCard(orootname_el_hp) );
              ybins_mu_hp.append( getPValueFromCard(orootname_mu_hp) );
              ybins_el_lp.append( getPValueFromCard(orootname_el_lp) );
              ybins_mu_lp.append( getPValueFromCard(orootname_mu_lp) );
              ybins_combo.append( getPValueFromCard(orootname_combo) );

              ybins_el_hp_exp.append( getPValueFromCard(orootname_el_hp_exp) );
              ybins_mu_hp_exp.append( getPValueFromCard(orootname_mu_hp_exp) );
              ybins_el_lp_exp.append( getPValueFromCard(orootname_el_lp_exp) );
              ybins_mu_lp_exp.append( getPValueFromCard(orootname_mu_lp_exp) );
              ybins_combo_exp.append( getPValueFromCard(orootname_combo_exp) ); 

         if options.plotPvalue : 

          gr_el_hp = ROOT.TGraph(nPoints,xbins,ybins_el_hp);
          gr_mu_hp = ROOT.TGraph(nPoints,xbins,ybins_mu_hp);
          gr_el_hp.SetLineColor( 1 ); gr_el_hp.SetMarkerColor( 1 ); gr_el_hp.SetMarkerStyle( 20 ); gr_el_hp.SetLineWidth( 3 );gr_el_hp.SetMarkerSize( 1.6 );
          gr_mu_hp.SetLineColor( 1 ); gr_mu_hp.SetMarkerColor( 1 ); gr_mu_hp.SetMarkerStyle( 20 ); gr_mu_hp.SetLineWidth( 3 );gr_mu_hp.SetMarkerSize( 1.6 );
          gr_el_lp = ROOT.TGraph(nPoints,xbins,ybins_el_lp);
          gr_mu_lp = ROOT.TGraph(nPoints,xbins,ybins_mu_lp);
          gr_el_lp.SetLineColor( 1 ); gr_el_lp.SetMarkerColor( 1 ); gr_el_lp.SetMarkerStyle( 20 ); gr_el_lp.SetLineWidth( 3 );gr_el_lp.SetMarkerSize( 1.6 );
          gr_mu_lp.SetLineColor( 1 ); gr_mu_lp.SetMarkerColor( 1 ); gr_mu_lp.SetMarkerStyle( 20 ); gr_mu_lp.SetLineWidth( 3 );gr_mu_lp.SetMarkerSize( 1.6 );
          gr_combo = ROOT.TGraph(nPoints,xbins,ybins_combo);
          gr_combo.SetLineColor( 1 ); gr_combo.SetMarkerColor( 1 ); gr_combo.SetMarkerStyle( 20 );gr_combo.SetLineWidth( 3 ); gr_combo.SetMarkerSize( 1.4 );

          gr_el_hp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_el_hp_exp);
          gr_mu_hp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_mu_hp_exp);
          gr_el_hp_exp.SetLineColor( 2 ); gr_el_hp_exp.SetMarkerColor( 2 ); gr_el_hp_exp.SetMarkerStyle( 20 ); gr_el_hp_exp.SetLineWidth( 3 );gr_el_hp_exp.SetMarkerSize( 1.6 );
          gr_el_hp_exp.SetLineStyle(8);
          gr_mu_hp_exp.SetLineColor( 2 ); gr_mu_hp_exp.SetMarkerColor( 2 ); gr_mu_hp_exp.SetMarkerStyle( 20 ); gr_mu_hp_exp.SetLineWidth( 3 );gr_mu_hp_exp.SetMarkerSize( 1.6 );
          gr_mu_hp_exp.SetLineStyle(8);
          gr_el_lp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_el_lp_exp);
          gr_mu_lp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_mu_lp_exp);
          gr_el_lp_exp.SetLineColor( 2 ); gr_el_lp_exp.SetMarkerColor( 2 ); gr_el_lp_exp.SetMarkerStyle( 20 ); gr_el_lp_exp.SetLineWidth( 3 );gr_el_lp_exp.SetMarkerSize( 1.6 );
          gr_el_lp_exp.SetLineStyle(8);
          gr_mu_lp_exp.SetLineColor( 2 ); gr_mu_lp_exp.SetMarkerColor( 2 ); gr_mu_lp_exp.SetMarkerStyle( 20 ); gr_mu_lp_exp.SetLineWidth( 3 );gr_mu_lp_exp.SetMarkerSize( 1.6 );
          gr_mu_lp_exp.SetLineStyle(8);
          gr_combo_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_combo_exp);
          gr_combo_exp.SetLineColor( 2 ); gr_combo_exp.SetMarkerColor( 2 ); gr_combo_exp.SetMarkerStyle( 20 );gr_combo_exp.SetLineWidth( 3 ); gr_combo_exp.SetMarkerSize( 1.6 ); 
          gr_combo_exp.SetLineStyle(8);

          oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",600,2500);
          oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
          twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",600,2500);
          twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
          threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",600,2500);
          threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
          fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",600,2500);
          fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);
    
          banner = TLatex(0.43,0.91,("CMS Preliminary, 19.7 fb^{-1} at #sqrt{s}=8TeV"));
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
          leg2.AddEntry( gr_el_hp_exp, "exp signif, #tilde{k}=0.5", "pl" );

          leg3 = ROOT.TLegend(0.2,0.2,0.5,0.35);
          leg3.SetFillStyle(0);
          leg3.SetBorderSize(1);
          leg3.AddEntry( gr_mu_hp, "obs signif, #mu (HP)", "pl" );
          leg3.AddEntry( gr_mu_hp_exp, "exp signif, #tilde{k}=0.5", "pl" );

          leg4 = ROOT.TLegend(0.2,0.2,0.5,0.35);
          leg4.SetFillStyle(0);
          leg4.SetBorderSize(1);
          leg4.AddEntry( gr_el_lp, "obs signif, e (LP)", "pl" );
          leg4.AddEntry( gr_el_lp_exp, "exp signif, #tilde{k}=0.5", "pl" );

          leg5 = ROOT.TLegend(0.2,0.2,0.5,0.35);
          leg5.SetFillStyle(0);
          leg5.SetBorderSize(1);
          leg5.AddEntry( gr_mu_lp, "obs signif, #mu (LP)", "pl" );
          leg5.AddEntry( gr_mu_lp_exp, "exp signif, #tilde{k}=0.5", "pl" );

          leg6 = ROOT.TLegend(0.2,0.2,0.5,0.35);
          leg6.SetFillStyle(0);
          leg6.SetBorderSize(1);
          leg6.AddEntry( gr_combo, "obs signif, e+#mu combined", "pl" );
          leg6.AddEntry( gr_combo_exp, "exp signif, #tilde{k}=0.5", "pl" );
    
          can = ROOT.TCanvas("can","can",800,800);
          hrl = can.DrawFrame(799,1e-3,2500,0.6);
          hrl.GetYaxis().SetTitle("p-value");
          hrl.GetXaxis().SetTitle("mass (GeV)");
          can.SetGrid();
          ROOT.gPad.SetLogy();
          gr_el_hp.Draw("PL");
          gr_el_hp_exp.Draw("PLsame");
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
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_HP.pdf","pdf");
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_HP.png","png");
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_HP.root","root");
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_HP.C","C");

          can2 = ROOT.TCanvas("can2","can2",800,800);
          hrl2 = can2.DrawFrame(799,1e-3,2500,0.6);
          hrl2.GetYaxis().SetTitle("p-value");
          hrl2.GetXaxis().SetTitle("mass (GeV)");
          can2.SetGrid();
          ROOT.gPad.SetLogy();
          gr_mu_hp.Draw("PL");
          gr_mu_hp_exp.Draw("PLsame");
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
          can2.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_HP.pdf","pdf");
          can2.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_HP.png","png");
          can2.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_HP.root","root");
          can2.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_HP.C","C");

          can3 = ROOT.TCanvas("can3","can3",800,800);
          hrl3 = can3.DrawFrame(799,1e-3,2500,0.6);
          hrl3.GetYaxis().SetTitle("p-value");
          hrl3.GetXaxis().SetTitle("mass (GeV)");
          can3.SetGrid();
          ROOT.gPad.SetLogy();
          gr_el_lp.Draw("PL");
          gr_el_lp_exp.Draw("PLsame");
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
          can3.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_LP.pdf","pdf");
          can3.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_LP.png","png");
          can3.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_LP.root","root");
          can3.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_el_LP.C","C");

          can4 = ROOT.TCanvas("can4","can4",800,800);
          hrl4 = can4.DrawFrame(799,1e-3,2500,0.6);
          hrl4.GetYaxis().SetTitle("p-value");
          hrl4.GetXaxis().SetTitle("mass (GeV)");
          can4.SetGrid();
          ROOT.gPad.SetLogy();
          gr_mu_lp.Draw("PL");
          gr_mu_lp_exp.Draw("PLsame");
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
          can4.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_LP.pdf","pdf");
          can4.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_LP.png","png");
          can4.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_LP.root","root");
          can4.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_mu_LP.C","C");

          can5 = ROOT.TCanvas("can5","can5",800,800);
          hrl5 = can5.DrawFrame(799,1e-3,2500,0.6);
          hrl5.GetYaxis().SetTitle("p-value");
          hrl5.GetXaxis().SetTitle("mass (GeV)");
          can5.SetGrid();
          ROOT.gPad.SetLogy();
          gr_combo.Draw("PL");
          gr_combo_exp.Draw("PLsame");
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
          can5.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_combo.pdf","pdf");
          can5.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_combo.png","png");
          can5.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_combo.root","root");
          can5.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_combo.C","C");


         doULPlot("_el_HP");
         doULPlot("_mu_HP");
         doULPlot("_el_LP");
         doULPlot("_mu_LP");
         doULPlot("_combo");

         ####################################################

        elif options.channel == "em":

         for i in range(mLo,mHi):
             print "############################"+str(mass[i])+"#############################";            

             print "getting card: higgsCombine_pval_obs_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             orootname_em_hp = "higgsCombine_pval_obs_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             print "getting card: higgsCombine_pval_exp_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             orootname_em_hp_exp = "higgsCombine_pval_exp_%03d_%s_HP.ProfileLikelihood.mH%03d.root"%(mass[i],"em",mass[i]);

             xbins.append( mass[i] );
             xbins_exp.append( mass[i] );
             if options.plotPvalue :
              ybins_em_hp.append( getPValueFromCard(orootname_em_hp) );
              ybins_em_hp_exp.append( getPValueFromCard(orootname_em_hp_exp) );

         if options.plotPvalue :
          gr_em_hp = ROOT.TGraph(nPoints,xbins,ybins_em_hp);
          gr_em_hp.SetLineColor( 1 ); gr_em_hp.SetMarkerColor( 1 ); gr_em_hp.SetMarkerStyle( 20 ); gr_em_hp.SetLineWidth( 3 );gr_em_hp.SetMarkerSize( 1.6 );

          gr_em_hp_exp = ROOT.TGraph(nPoints,xbins_exp,ybins_em_hp_exp);
          gr_em_hp_exp.SetLineColor( 2 ); gr_em_hp_exp.SetMarkerColor( 2 ); gr_em_hp_exp.SetMarkerStyle( 20 ); gr_em_hp_exp.SetLineWidth( 3 );gr_em_hp_exp.SetMarkerSize( 1.6 );

          oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",600,2500);
          oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
          twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",600,2500);
          twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
          threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",600,2500);
          threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
          fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",600,2500);
          fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);
    
          banner = TLatex(0.43,0.91,("CMS Preliminary, 19.7 fb^{-1} at #sqrt{s}=8TeV"));
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
          leg2.AddEntry( gr_em_hp, "obs signif, e+#mu (HP)", "pl" );
          leg2.AddEntry( gr_em_hp_exp, "exp signif, #tilde{k}=0.5", "pl" );

    
          can = ROOT.TCanvas("can","can",800,800);
          hrl = can.DrawFrame(799,1e-3,2500,0.6);
          hrl.GetYaxis().SetTitle("p-value");
          hrl.GetXaxis().SetTitle("mass (GeV)");
          can.SetGrid();
          ROOT.gPad.SetLogy();
          gr_em_hp.Draw("PL");
          gr_em_hp_exp.Draw("PLsame");
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
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_em_HP.pdf","pdf");
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_em_HP.png","png");
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_em_HP.root","root");
          can.SaveAs("~/LimitResult/Limit_ExpTail_modelIndependent/pvals_em_HP.C","C");

         doULPlot("_em_HP");
