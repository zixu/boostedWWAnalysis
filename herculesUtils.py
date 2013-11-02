import os
import glob
import math
import time
import subprocess
from subprocess import Popen
from optparse import OptionParser

import ROOT
from ROOT import *
from ROOT import gROOT, gStyle, gSystem, TLatex

#### to execute the scripts 
#### python herculesUtils.py -b --batchMode --createTree --sampleToProcess all  --channel mu --queque shortcms
#### python herculesUtils.py -b --batchMode --createTree --sampleToProcess all  --channel el --queque shortcms

############################################
#            Job steering                  #
############################################
parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

parser.add_option('--createTree', action='store_true', dest='createTree', default=False,help='no X11 windows')

# submit jobs to condor
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False, help='no X11 windows')
parser.add_option('--sampleToProcess',action="store",type="string",dest="sampleToProcess",default=None)
parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--queque',action="store",type="string",dest="queque",default="longcms")


(options, args) = parser.parse_args()

######################################


def submitBatchJob( command, fn ):

    currentDir = os.getcwd();
    # create a dummy bash/csh
    outScript=open(fn+".sh","w");
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+command);
    outScript.close();

    # submit the condor job 

    Lunch = "lancia_%s"%options.channel;
    outLunch=open(Lunch+".sh","a");
    outLunch.write("qsub -V -d "+currentDir+" -q "+options.queque+" "+currentDir+"/"+fn+".sh \n") ;
    

if __name__ == '__main__':

    all = [ "data",
            "ggH600","ggH700","ggH800","ggH900","ggH1000",
            "vbfH600","vbfH700","vbfH800","vbfH900","vbfH1000",
#           "rsg1000_kMpl01_py","rsg1000_kMpl01_hw","rsg1500_kMpl01_py","rsg1500_kMpl01_hw","rsg2000_kmpl01_py",
            "BulkG_c0p2_M600","BulkG_c0p2_M700","BulkG_c0p2_M800","BulkG_c0p2_M900","BulkG_c0p2_M1000","BulkG_c0p2_M1100","BulkG_c0p2_M1200","BulkG_c0p2_M1300",
            "BulkG_c0p2_M1400","BulkG_c0p2_M1500","BulkG_c0p2_M1600","BulkG_c0p2_M1700","BulkG_c0p2_M1800","BulkG_c0p2_M1900","BulkG_c0p2_M2000","BulkG_c0p2_M2100",
            "BulkG_c0p2_M2200","BulkG_c0p2_M2300","BulkG_c0p2_M2400","BulkG_c0p2_M2500",
            "WJets_Pythia","WJets_Pythia_matchDw","WJets_Pythia_matchUp","WJets_Pythia_scaleUp","WJets_Pythia_scaleDw",
            "WJets_Herwig","WJets_Pythia180","WJets_Pythia100",
            "ZJets","WW","WZ","ZZ",
            "TTbar","TTbar_Powheg","TTbar_matchDn","TTbar_matchUp","TTbar_scaleDn","TTbar_scaleUp","TTbar_mcatnlo",
            "tch","tWch","sch","tch_bar","tWch_bar","sch_bar",
             "X_WW_lvjj_PS1000m","X_WW_lvjj_PS2000m","X_WW_lvjj_PS600m","X_WW_lvjj_PS600p",
             "X_WW_lvjj_SM1000m","X_WW_lvjj_SM2000m","X_WW_lvjj_SM600m","X_WW_lvjj_SM600p"
           ]


    sourcefiledirectory = "/gwteray/users/gerosa/VBF_Trees/VBF_Trees_v1/";
    outputfiledirectory = "/gwteray/users/gerosa/otrees/otrees_VBF_v1/";

    os.system("rm lancia_%s.sh"%options.channel);
    os.system("touch lancia_%s.sh"%options.channel);
        
    if options.sampleToProcess == "all" and options.batchMode and options.createTree:
    
        for i in range(len(all)):
            
            cmmd = "python runAnalysis.py -b --createTrees --channel "+options.channel+" --sampleToProcess "+ all[i] +" --sourcefiledirectory "+sourcefiledirectory+" --outputfiledirectory "+outputfiledirectory;
            print cmmd
            fn = "createTreeScript_%s_%s"%(options.channel,all[i]);
            submitBatchJob( cmmd, fn );

    
    elif options.createTree and options.batchMode and not options.sampleToProcess == None:

        cmmd = "python runAnalysis.py -b --createTrees --channel "+options.channel+" --sampleToProcess "+ options.sampleToProcess + " --sourcefiledirectory "+sourcefiledirectory+" --outputfiledirectory "+outputfiledirectory;
        print cmmd
        fn = "createTreeScript_%s_%s"%(options.channel,options.sampleToProcess);
        submitBatchJob( cmmd, fn );

    else:
        print "do nothing"
