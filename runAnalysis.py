#!/usr/bin/env python
import os
import glob
import math
import array
import subprocess
from subprocess import Popen

from ROOT import gROOT, gStyle, gSystem, TLatex

from sampleWrapperClass import *
from BoostedWSamples import * 
from trainingClass import *

############################################
#            Job steering                  #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False,help='no X11 windows')

parser.add_option('--createTrees', action='store_true', dest='createTrainingTrees', default=False,help='creates Training Trees')

parser.add_option('--makeControlPlots', action='store_true', dest='makeControlPlots', default=False, help='makeControlPlots')
parser.add_option('--makeTMVAPlots', action='store_true', dest='makeTMVAPlots', default=False, help='makeTMVAPlots')
parser.add_option('--makeSignalRegionControlPlots', action='store_true', dest='makeSignalRegionControlPlots', default=False,help='makeSignalRegionControlPlots')
parser.add_option('--makeTTBarControlPlots', action='store_true', dest='makeTTBarControlPlots', default=False, help='makeTTBarControlPlots')
parser.add_option('--makeFinalTree', action='store_true', dest='makeFinalTree', default=False, help='make Final Tree')
parser.add_option('--doTraining', action='store_true', dest='doTraining', default=False,help='does training')
parser.add_option('-m', '--trainingMethod',action="store",type="string",dest="trainingMethod",default="Likelihood")


parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--sampleToProcess',action="store",type="string",dest="sampleToProcess",default=None)

parser.add_option('-i', '--graphindex',action="store",type="int",dest="graphindex",default=0)

parser.add_option('--sourcefiledirectory',action="store",type="string",dest"sourcefiledirectory",defaul="/gwteray/users/gerosa/RD_Trees_Production/RD_Trees_ExoLeptonID/")
parser.add_option('--outputfiledirectory',action="store",type="string",dest"outputfiledirectory",defaul="/gwteray/users/gerosa/RD_Trees_Production/RD_Trees_ExoLeptonID/")

(option, args) = parser.parse_args()

############################################################

if __name__ == '__main__':

    
    print "Welcome to the boosted analysis..."

    sourcefiledirectory = options.sourcefiledirectory;    
    outputfiledirectory = options.outputfiledirectory;    

    # ---------------------------------------------------
    # check if directories exists
    if not os.path.isdir(outputfiledirectory+"trainingtrees_mu"): os.system("mkdir "+outputfiledirectory+"trainingtrees_mu");
    if not os.path.isdir(outputfiledirectory+"trainingtrees_el"): os.system("mkdir "+outputfiledirectory+"trainingtrees_el");
    if not os.path.isdir("classifier"): os.system("mkdir classifier");    

    # ---------------------------------------------------
    # define samples, this creates some trees in the "trainingtrees" directory

    isData = True;
    notData = False;
    CHANNEL = options.channel
    if CHANNEL == 'el':
        LUMI = 19.2;
    elif CHANNEL == 'mu':
        LUMI = 19.3;


    # --------------------------------------------------
    treename = ""
    
    if options.makeControlPlots or options.makeTTBarControlPlots: 
       treename = "WJet"
    if options.makeTMVAPlots:
       sourcefiledirectory = "/uscms_data/d3/weizou/VBFHiggsAnalysis/BoostedWAnalysis2012/boostedWWAnalysis/trainingtrees_%s/"%(CHANNEL)
       treename = "otree"
    lumifile = "MCScaleFactors.txt"

    # --------------------------------------------------- take the root RD files from what is passed to the main
    boostedWSamples = Samples(CHANNEL)
    boostedWSamples.SetFilePath(sourcefiledirectory)
    boostedWSamples.SetTreeName(treename)
    boostedWSamples.SetFileNames()
    boostedWSamples.SetLumi(LUMI)

    # ----------------------------------------------------
    
    datasample = sampleWrapperClass("data",boostedWSamples.GetFileNames()["data"],CHANNEL,LUMI,LUMI,boostedWSamples.GetTreeName(),isData,outputfiledirectory)
    
    ggH600Sample = sampleWrapperClass("ggH600",boostedWSamples.GetFileNames()["ggH600"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ggH600")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
    
    ggH700Sample = sampleWrapperClass("ggH700",boostedWSamples.GetFileNames()["ggH700"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ggH700")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
    
    ggH800Sample = sampleWrapperClass("ggH800",boostedWSamples.GetFileNames()["ggH800"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ggH800")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
    
    ggH900Sample = sampleWrapperClass("ggH900",boostedWSamples.GetFileNames()["ggH900"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ggH900")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
    
    ggH1000Sample = sampleWrapperClass("ggH1000",boostedWSamples.GetFileNames()["ggH1000"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ggH1000")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    vbfH600Sample = sampleWrapperClass("vbfH600",boostedWSamples.GetFileNames()["vbfH600"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"vbfH600")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    vbfH700Sample = sampleWrapperClass("vbfH700",boostedWSamples.GetFileNames()["vbfH700"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"vbfH700")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    vbfH800Sample = sampleWrapperClass("vbfH800",boostedWSamples.GetFileNames()["vbfH800"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"vbfH800")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    vbfH900Sample = sampleWrapperClass("vbfH900",boostedWSamples.GetFileNames()["vbfH900"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"vbfH900")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    vbfH1000Sample = sampleWrapperClass("vbfH1000",boostedWSamples.GetFileNames()["vbfH1000"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"vbfH1000")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    rsg1000Sample_kMpl01_py = sampleWrapperClass("rsg1000_kMpl01_py",boostedWSamples.GetFileNames()["rsg1000_kMpl01_py"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"rsg1000_kMpl01_py")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    rsg1000Sample_kMpl01_hw = sampleWrapperClass("rsg1000_kMpl01_hw",boostedWSamples.GetFileNames()["rsg1000_kMpl01_hw"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"rsg1000_kMpl01_hw")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    rsg1500Sample_kMpl01_py = sampleWrapperClass("rsg1500_kMpl01_py",boostedWSamples.GetFileNames()["rsg1500_kMpl01_py"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"rsg1500_kMpl01_py")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    rsg1500Sample_kMpl01_hw = sampleWrapperClass("rsg1500_kMpl01_hw",boostedWSamples.GetFileNames()["rsg1500_kMpl01_hw"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"rsg1500_kMpl01_hw")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    rsg2000Sample_kMpl01_py = sampleWrapperClass("rsg2000_kMpl01_py",boostedWSamples.GetFileNames()["rsg2000_kMpl01_py"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"rsg2000_kMpl01_py")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M600 =  sampleWrapperClass("BulkG_c0p2_M600",boostedWSamples.GetFileNames()["BulkG_c0p2_M600"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M600")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M700 =  sampleWrapperClass("BulkG_c0p2_M700",boostedWSamples.GetFileNames()["BulkG_c0p2_M700"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M700")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M800 =  sampleWrapperClass("BulkG_c0p2_M800",boostedWSamples.GetFileNames()["BulkG_c0p2_M800"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M800")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M900 =  sampleWrapperClass("BulkG_c0p2_M900",boostedWSamples.GetFileNames()["BulkG_c0p2_M900"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M900")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1000 =  sampleWrapperClass("BulkG_c0p2_M1000",boostedWSamples.GetFileNames()["BulkG_c0p2_M1000"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1000")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1100 =  sampleWrapperClass("BulkG_c0p2_M1100",boostedWSamples.GetFileNames()["BulkG_c0p2_M1100"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1100")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1200 =  sampleWrapperClass("BulkG_c0p2_M1200",boostedWSamples.GetFileNames()["BulkG_c0p2_M1200"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1200")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1300 =  sampleWrapperClass("BulkG_c0p2_M1300",boostedWSamples.GetFileNames()["BulkG_c0p2_M1300"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1300")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1400 =  sampleWrapperClass("BulkG_c0p2_M1400",boostedWSamples.GetFileNames()["BulkG_c0p2_M1400"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1400")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1500 =  sampleWrapperClass("BulkG_c0p2_M1500",boostedWSamples.GetFileNames()["BulkG_c0p2_M1500"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1500")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1600 =  sampleWrapperClass("BulkG_c0p2_M1600",boostedWSamples.GetFileNames()["BulkG_c0p2_M1600"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1600")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1700 =  sampleWrapperClass("BulkG_c0p2_M1700",boostedWSamples.GetFileNames()["BulkG_c0p2_M1700"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1700")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M1800 =  sampleWrapperClass("BulkG_c0p2_M1800",boostedWSamples.GetFileNames()["BulkG_c0p2_M1800"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1800")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
    
    BulkG_c0p2_M1900 =  sampleWrapperClass("BulkG_c0p2_M1900",boostedWSamples.GetFileNames()["BulkG_c0p2_M1900"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M1900")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M2000 =  sampleWrapperClass("BulkG_c0p2_M2000",boostedWSamples.GetFileNames()["BulkG_c0p2_M2000"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M2000")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M2100 =  sampleWrapperClass("BulkG_c0p2_M2100",boostedWSamples.GetFileNames()["BulkG_c0p2_M2100"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M2100")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
  
    BulkG_c0p2_M2200 =  sampleWrapperClass("BulkG_c0p2_M2200",boostedWSamples.GetFileNames()["BulkG_c0p2_M2200"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M2200")),    LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M2300 =  sampleWrapperClass("BulkG_c0p2_M2300",boostedWSamples.GetFileNames()["BulkG_c0p2_M2300"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M2300")),    LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M2400 =  sampleWrapperClass("BulkG_c0p2_M2400",boostedWSamples.GetFileNames()["BulkG_c0p2_M2400"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M2400")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    BulkG_c0p2_M2500 =  sampleWrapperClass("BulkG_c0p2_M2500",boostedWSamples.GetFileNames()["BulkG_c0p2_M2500"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"BulkG_c0p2_M2500")), LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)


    WJets_PythiaSample = sampleWrapperClass("WJets_Pythia",boostedWSamples.GetFileNames()["WJets_Pythia"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_Pythia_matchDwSample = sampleWrapperClass("WJets_Pythia_matchDw",boostedWSamples.GetFileNames()["WJets_Pythia_matchDw"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia_matchDw")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_Pythia_matchUpSample = sampleWrapperClass("WJets_Pythia_matchUp",boostedWSamples.GetFileNames()["WJets_Pythia_matchUp"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia_matchUp")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_Pythia_scaleDwSample = sampleWrapperClass("WJets_Pythia_scaleDw",boostedWSamples.GetFileNames()["WJets_Pythia_scaleDw"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia_scaleDw")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_Pythia_scaleUpSample = sampleWrapperClass("WJets_Pythia_scaleUp",boostedWSamples.GetFileNames()["WJets_Pythia_scaleUp"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia_scaleUp")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    W1Jets_Pythia = sampleWrapperClass("W1Jets_Pythia",boostedWSamples.GetFileNames()["W1Jets_Pythia"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"W1Jets_Pythia")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    W2Jets_Pythia = sampleWrapperClass("W2Jets_Pythia",boostedWSamples.GetFileNames()["W2Jets_Pythia"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"W2Jets_Pythia")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    W3Jets_Pythia = sampleWrapperClass("W3Jets_Pythia",boostedWSamples.GetFileNames()["W3Jets_Pythia"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"W3Jets_Pythia")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    W4Jets_Pythia = sampleWrapperClass("W4Jets_Pythia",boostedWSamples.GetFileNames()["W4Jets_Pythia"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"W4Jets_Pythia")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_Pythia100Sample = sampleWrapperClass("WJets_Pythia100",boostedWSamples.GetFileNames()["WJets_Pythia100"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia100")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_Pythia180Sample = sampleWrapperClass("WJets_Pythia180",boostedWSamples.GetFileNames()["WJets_Pythia180"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Pythia180")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    WJets_HerwigSample = sampleWrapperClass("WJets_Herwig",boostedWSamples.GetFileNames()["WJets_Herwig"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WJets_Herwig")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)
    
    ZJetsSample = sampleWrapperClass("ZJets",boostedWSamples.GetFileNames()["ZJets"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ZJets")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory)

    TTbarSample = sampleWrapperClass("TTbar",boostedWSamples.GetFileNames()["TTbar"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"TTbar")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);
    
    TTbarSample_matchDn = sampleWrapperClass("TTbar_matchDn",boostedWSamples.GetFileNames()["TTbar_matchDn"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"TTbar_matchDn")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    TTbarSample_matchUp = sampleWrapperClass("TTbar_matchUp",boostedWSamples.GetFileNames()["TTbar_matchUp"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"TTbar_matchUp")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    TTbarSample_Powheg = sampleWrapperClass("TTbar_Powheg",boostedWSamples.GetFileNames()["TTbar_Powheg"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"TTbar_Powheg")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    TTbarSample_scaleDn = sampleWrapperClass("TTbar_scaleDn",boostedWSamples.GetFileNames()["TTbar_scaleDn"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"TTbar_scaleDn")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    TTbarSample_scaleUp = sampleWrapperClass("TTbar_scaleUp",boostedWSamples.GetFileNames()["TTbar_scaleUp"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"TTbar_scaleUp")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    WWSample = sampleWrapperClass("WW",boostedWSamples.GetFileNames()["WW"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WW")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    WZSample = sampleWrapperClass("WZ",boostedWSamples.GetFileNames()["WZ"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"WZ")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    ZZSample = sampleWrapperClass("ZZ",boostedWSamples.GetFileNames()["ZZ"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"ZZ")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);
   
    tchSample = sampleWrapperClass("tch",boostedWSamples.GetFileNames()["tch"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"tch")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    tWchSample = sampleWrapperClass("tWch",boostedWSamples.GetFileNames()["tWch"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"tWch")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    schSample = sampleWrapperClass("sch",boostedWSamples.GetFileNames()["sch"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"sch")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    tch_barSample = sampleWrapperClass("tch_bar",boostedWSamples.GetFileNames()["tch_bar"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"tch_bar")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    tWch_barSample = sampleWrapperClass("tWch_bar",boostedWSamples.GetFileNames()["tWch_bar"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"tWch_bar")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);

    sch_barSample = sampleWrapperClass("sch_bar",boostedWSamples.GetFileNames()["sch_bar"],CHANNEL,1.0/(boostedWSamples.GetLumiScaleFactor(lumifile,"sch_bar")),LUMI,boostedWSamples.GetTreeName(),notData,outputfiledirectory);


    allSamples = [datasample,
                  ggH600Sample,ggH700Sample,ggH800Sample,ggH900Sample,ggH1000Sample,
                  vbfH600Sample,vbfH700Sample,vbfH800Sample,vbfH900Sample,vbfH1000Sample,
                  rsg1000Sample_kMpl01_py,rsg1000Sample_kMpl01_hw,rsg1500Sample_kMpl01_py,rsg1500Sample_kMpl01_hw,rsg2000Sample_kMpl01_py,
                  BulkG_c0p2_M600,BulkG_c0p2_M700,BulkG_c0p2_M800,BulkG_c0p2_M900,BulkG_c0p2_M1000,BulkG_c0p2_M1100,BulkG_c0p2_M1200,
                  BulkG_c0p2_M1300,BulkG_c0p2_M1400,BulkG_c0p2_M1500,BulkG_c0p2_M1600,BulkG_c0p2_M1700,BulkG_c0p2_M1800,BulkG_c0p2_M1900,BulkG_c0p2_M2000,
                  BulkG_c0p2_M2100,BulkG_c0p2_M2200,BulkG_c0p2_M2300,BulkG_c0p2_M2400,BulkG_c0p2_M2500,
                  WJets_PythiaSample,WJets_Pythia_matchUpSample,WJets_Pythia_matchDwSample,WJets_Pythia_scaleUpSample,WJets_Pythia_scaleDwSample,W1Jets_PythiaSample,
                  W2Jets_PythiaSample,W3Jets_PythiaSample,W4Jets_PythiaSample,WJets_Pythia100Sample
                  WJets_Pythia180Sample,WJets_HerwigSample,ZJetsSample,TTbarSample,
                  TTbarSample_matchDn,TTbarSample_matchUp,TTbarSample_scaleDn,TTbarSample_scaleUp,TTbarSample_Powheg,
                  WWSample,WZSample,ZZSample,
                  tchSample,tWchSample,schSample,tch_barSample,tWch_barSample,sch_barSample];


    if options.createTrainingTrees:
        
        # ---------------------------------------------------
        # create training tree

        if options.sampleToProcess == None:

            datasample.createTrainingTree();
            WJets_Pythia180Sample.createTrainingTree();
            WJets_Pythia100Sample.createTrainingTree();

            ZJetsSample.createTrainingTree();

            TTbarSample.createTrainingTree();
            TTbarSample_matchDn.createTrainingTree();
            TTbarSample_matchUp.createTrainingTree();
            TTbarSample_Powheg.createTrainingTree();
            TTbarSample_scaleDn.createTrainingTree();            
            TTbarSample_scaleUp.createTrainingTree();            

            WWSample.createTrainingTree();
            WZSample.createTrainingTree();
            ZZSample.createTrainingTree();

            tchSample.createTrainingTree();
            tWchSample.createTrainingTree();
            schSample.createTrainingTree();
            tch_barSample.createTrainingTree();
            tWch_barSample.createTrainingTree();
            sch_barSample.createTrainingTree();

        else:
            for i in range(len(allSamples)):
                if options.sampleToProcess == allSamples[i].getLabel():
                    allSamples[i].createTrainingTree();


    if options.makeControlPlots:
                
        # ---------------------------------------------------
        # make control plots based on the same cuts as the training trees

        if not os.path.isdir("controlPlots"): os.system("mkdir controlPlots");
        #myPlotter.makeControlPlots("controlPlots","nocuts");

        print "Please Check the Cuts used on the BoostedWControlPlots.py is reasonable"
        Cuts = "W_pt > 200 && GroomedJet_CA8_pt[0] > 200 && ggdboostedWevt == 1 && event_met_pfmet > 50 && GroomedJet_CA8_deltaphi_METca8jet > 2.0 && GroomedJet_CA8_deltaR_lca8jet > 1.57";
        if CHANNEL == "el": Cuts = Cuts + " && event_met_pfmet > 70 && W_electron_pt > 35";   
        if CHANNEL == "mu": Cuts = Cuts + " && W_muon_pt > 30"; 

        Cuts = Cuts + " && GroomedJet_numberbjets_csvm > 0";
        
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt_pr[0] > 180 && ggdboostedWevt == 1 && event_metMVA_met > 50"
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt_pr[0] > 180 && ggdboostedWevt == 1 && event_metMVA_met > 50 && GroomedJet_numberjets <= 1"
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt_pr[0] > 180 && ggdboostedWevt == 1 && event_metMVA_met > 50 && GroomedJet_numberbjets == 0"

        print "Cuts we apply: " + Cuts

        if options.noX:
           p = subprocess.Popen(["python","BoostedWControlPlots.py","-b","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts,"-n","%s"%CHANNEL])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )
        else:
           p = subprocess.Popen(["python","BoostedWControlPlots.py","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts,"-n","%s"%CHANNEL])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )
    
    if options.makeTTBarControlPlots:        

        if not os.path.isdir("controlPlots_ttbar"): os.system("mkdir controlPlots_ttbar");
        print "Please Check the Cuts used on the BoostedWControlPlots.py is reasonable"
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt[0] > 200 && ggdboostedWevt == 1 && event_metMVA_met > 50 && GroomedJet_CA8_tau2tau1[0]< 0.524   "
        Cuts = "W_pt > 180 && GroomedJet_CA8_pt[0] > 200 && ggdboostedWevt == 1 && event_metMVA_met > 50 && GroomedJet_CA8_tau2tau1[0]< 0.524  && GroomedJet_numberbjets >= 1 "
        addition  = "ttbar";
        print "Cuts we apply: " + Cuts
        print "We don't put the cuts on the signal and try to compare the W jet performance with TTbar and Data!!!!"
        print "Make Sure the numberbjets cut is the last cut in the cut sequence"
        if options.noX:
           p = subprocess.Popen(["python","BoostedWControlPlots.py","-b","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts,"-n","%s"%CHANNEL,"-a",addition])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )
        else:
           p = subprocess.Popen(["python","BoostedWControlPlots.py","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts,"-n","%s"%CHANNEL,"-a",addition])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )


    if options.makeSignalRegionControlPlots:        
        if not os.path.isdir("controlPlots_signalregion"): os.system("mkdir controlPlots_signalregion");
        print "Please Check the Cuts used on the BoostedWControlPlots.py is reasonable"
        Cuts = "W_pt > 180 && GroomedJet_CA8_pt[0] > 200 && ggdboostedWevt == 1 && event_metMVA_met > 50  && GroomedJet_CA8_mass_pr[0]>=70 && GroomedJet_CA8_mass_pr[0]<=100 && boostedW_lvj_m>400 && boostedW_lvj_m<1400 && GroomedJet_CA8_tau2tau1[0]< 0.525 && GroomedJet_numberjets <1"
        addition  = "signalregion";
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt[0] > 200 && ggdboostedWevt == 1 && event_metMVA_met > 50  && GroomedJet_CA8_mass_pr[0]>=70 && GroomedJet_CA8_mass_pr[0]<=100"
        #Wtagger = "GroomedJet_CA8_tau2tau1[0]< 0.524"
        #Cuts = Cuts + " && " + Wtagger;
        #addition  = "signalregion_withoutjetnumber";
        print "Cuts we apply: " + Cuts
        if options.noX:
           p = subprocess.Popen(["python","BoostedWControlPlots.py","-b","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts,"-n","%s"%CHANNEL,"-a",addition])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )
        else:
           p = subprocess.Popen(["python","BoostedWControlPlots.py","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts,"-n","%s"%CHANNEL,"-a",addition])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )
     
    if options.makeTMVAPlots:
                
        # ---------------------------------------------------
        # make control plots based on the same cuts as the training trees

        if not os.path.isdir("controlPlots"): os.system("mkdir controlPlots");
        #myPlotter.makeControlPlots("controlPlots","nocuts");
        print "Please Check the Cuts used on the BoostedWTMVAPlots.py is reasonable"
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt_pr[0] > 180 && ggdboostedWevt == 1 && event_metMVA_met > 50"
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt_pr[0] > 180 && ggdboostedWevt == 1 && event_metMVA_met > 50 && GroomedJet_numberjets <= 1"
        #Cuts = "W_pt > 180 && GroomedJet_CA8_pt_pr[0] > 180 && ggdboostedWevt == 1 && event_metMVA_met > 50 && GroomedJet_numberbjets == 0"
        Cuts = "jet_mass_pr > 60 && jet_mass_pr < 100"
        #Cuts = "1 > 0"
        print "Cuts we apply: " + Cuts
        if options.noX:
           p = subprocess.Popen(["python","BoostedWTMVAPlots.py","-b","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )
        else:
           p = subprocess.Popen(["python","BoostedWTMVAPlots.py","-f","%s"%sourcefiledirectory,"-t","%s"%treename,"-l","%f"%LUMI,"-s","%s"%lumifile,"-c","%s"%Cuts])
           if(p.wait() != None): raw_input( 'Press ENTER to continue\n ' )

    # ---------------------------------------------------
    # do the training
    # get the training tree names

    if options.doTraining:

     signalTrainingTreeName = ggH600Sample.getTrainingTreeName();
     backgroundTrainingTreeNames = WJets_PythiaSample.getTrainingTreeName();
    
     trainingMethod = options.trainingMethod

     # Trainings
     # --------- #1 -----------
     listOfTrainingVariables1 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1"];
     WWTraining1 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables1, "simple" );
     # --------- #2 -----------
     listOfTrainingVariables2 = ["jet_grsens_tr","jet_grsens_ft","jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore7","jet_planarflow04","jet_planarflow07"];
     WWTraining2 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables2, "no subjets" );
     listOfTrainingVariables2 = ["jet_grsens_tr","jet_grsens_ft","jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_sjdr","jet_rcore4","jet_rcore7","jet_planarflow04","jet_planarflow07"];
     WWTraining2 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables2, "optimal" );

     listOfTrainingVariables3 = ["jet_massdrop_pr"]
     WWTraining3 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables3, trainingMethod, "massdrop")

     listOfTrainingVariables4 = ["jet_massdrop_pr", "jet_tau2tau1"]
     WWTraining4 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables4, "massdroptau2tau1")

     listOfTrainingVariables5 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4"]
     WWTraining5 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables5, "simplercore4")

     listOfTrainingVariables6 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5"]
     WWTraining6 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables6, "simplercore45")
 
     listOfTrainingVariables7 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6"]
     WWTraining7 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables7, "simplercore456")

     listOfTrainingVariables8 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7"]
     WWTraining8 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables8, "simplercore4567")

     listOfTrainingVariables9 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04"]
     WWTraining9 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables9, "simplercore4567planarflow4")

     listOfTrainingVariables10 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05"]
     WWTraining10 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables10, "simplercore4567planarflow45")

     listOfTrainingVariables11 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06"]
     WWTraining11 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables11, "simplercore4567planarflow456")

     listOfTrainingVariables12 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06","jet_planarflow07"]
     WWTraining12 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables12, "simplercore4567planarflow4567")

     listOfTrainingVariables13 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06","jet_planarflow07","jet_pt1frac"]
     WWTraining13 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables13, "simplercore4567planarflow4567subjet1")

     listOfTrainingVariables14 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06","jet_planarflow07","jet_pt2frac"]
     WWTraining14 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables14, "simplercore4567planarflow4567subjet2")

     listOfTrainingVariables15 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06","jet_planarflow07","jet_pt1frac","jet_sjdr"]
     WWTraining15 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables15, "simplercore4567planarflow4567subjet1sjdr")

     listOfTrainingVariables16 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06","jet_planarflow07","jet_pt1frac","jet_sjdr","jet_grsens_ft","jet_grsens_tr"]
     WWTraining16 = trainingClass(signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables16, "simplercore4567planarflow4567subjet1sjdrgrsens_ftgrsens_tr")

     listOfTrainingVariables17 = ["jet_grsens_tr","jet_grsens_ft","jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_rcore7","jet_planarflow07"];
     WWTraining17 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables3, "optimal2" );

     listOfTrainingVariables18 = ["jet_grsens_tr","jet_grsens_ft","jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_pt1frac","jet_sjdr","jet_jetconstituents","jet_rcore4","jet_rcore5","jet_rcore6","jet_rcore7","jet_planarflow04","jet_planarflow05","jet_planarflow06","jet_planarflow07"];
     WWTraining18 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables4, "all" );

     listOfTrainingVariables19 = ["jet_massdrop_pr","jet_qjetvol"];
     WWTraining19 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables17, "massdropqjet" );
    
     listOfTrainingVariables20 = ["jet_qjetvol","jet_tau2tau1"];
     WWTraining20 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables18,trainingMethod,"tau2tau1qjet" );
     
     listOfTrainingVariables21 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_pt1frac"];
     WWTraining21 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables19, "simplesubjet1pt" );

     listOfTrainingVariables22 = ["jet_massdrop_pr","jet_qjetvol","jet_tau2tau1","jet_pt1frac","jet_sjdr"];
     WWTraining22 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables20, "simplesubjet1ptsjdr" );

     listOfTrainingVariables23 = ["jet_massdrop_pr","jet_qjetvol","jet_pt1frac","jet_sjdr"];
     WWTraining213= trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables21, "massdropqjetsubjet1ptsjdr" );

     listOfTrainingVariables23 = ["jet_qjetvol"];
     WWTraining23 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables22,trainingMethod,"qjet" );

     listOfTrainingVariables24 = ["jet_tau2tau1"];
     WWTraining24 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables23,trainingMethod, "tau2tau1" );

     listOfTrainingVariables25 = ["jet_qjetvol","jet_tau2tau1","jet_pt1frac","jet_sjdr"];
     WWTraining25 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables24, trainingMethod, "tau2tau1qjetsubjet1ptsjdr" );
   
     listOfTrainingVariables26 = ["jet_massdrop_pr","jet_qjetvol"];
     WWTraining26 = trainingClass( signalTrainingTreeName, backgroundTrainingTreeNames, listOfTrainingVariables25,trainingMethod,"massdropqjet" );

     PrimaryVertex1 = 0
     PrimaryVertex2 = 15
     PrimaryVertex3 = 30
     PrimaryVertex4 = 40

     WWTraining18.doTraining( 200, 275, PrimaryVertex1, PrimaryVertex2);
     WWTraining18.doTraining( 200, 275, PrimaryVertex2, PrimaryVertex3);
     WWTraining18.doTraining( 200, 275, PrimaryVertex3, PrimaryVertex4);
     WWTraining18.doTraining( 275, 500, PrimaryVertex1, PrimaryVertex2);
     WWTraining18.doTraining( 275, 500, PrimaryVertex2, PrimaryVertex3);
     WWTraining18.doTraining( 275, 500, PrimaryVertex3, PrimaryVertex4);
     WWTraining18.plotTrainingResults( 200, 275, PrimaryVertex1, PrimaryVertex2);
     WWTraining18.plotTrainingResults( 200, 275, PrimaryVertex2, PrimaryVertex3);
     WWTraining18.plotTrainingResults( 200, 275, PrimaryVertex3, PrimaryVertex4);
     WWTraining18.plotTrainingResults( 275, 500, PrimaryVertex1, PrimaryVertex2);
     WWTraining18.plotTrainingResults( 275, 500, PrimaryVertex2, PrimaryVertex3);
     WWTraining18.plotTrainingResults( 275, 500, PrimaryVertex3, PrimaryVertex4);


     WWTraining18.doTraining( 200, 275, 0, 100);
     WWTraining18.doTraining( 275, 500, 0, 100);
        


     WWTraining18.plotTrainingResults( 200, 275, 0, 100);
     WWTraining18.plotTrainingResults( 275, 500, 0, 100);       


     WWTraining1.doTraining( 200, 275 );
     WWTraining1.doTraining( 275, 500 );

     WWTraining1.doTraining( 200, 275, 0,  10);
     WWTraining1.doTraining( 200, 275, 10, 20);
     WWTraining1.doTraining( 200, 275, 20, 30);
     WWTraining1.doTraining( 275, 500, 0,  10);
     WWTraining1.doTraining( 275, 500, 10, 20);
     WWTraining1.doTraining( 275, 500, 20, 30);


     WWTraining22.doTraining( 200, 275, PrimaryVertex1, PrimaryVertex2);
     WWTraining22.doTraining( 200, 275, PrimaryVertex2, PrimaryVertex3);
     WWWTraining22.doTraining( 200, 275, PrimaryVertex3, PrimaryVertex4);
     WWTraining22.doTraining( 275, 500, PrimaryVertex1, PrimaryVertex2);
     WWTraining22.doTraining( 275, 500, PrimaryVertex2, PrimaryVertex3);
     WWTraining22.doTraining( 275, 500, PrimaryVertex3, PrimaryVertex4);

     WWTraining23.doTraining( 200, 275, PrimaryVertex1, PrimaryVertex2);
     WWTraining23.doTraining( 200, 275, PrimaryVertex2, PrimaryVertex3);
     WWTraining23.doTraining( 200, 275, PrimaryVertex3, PrimaryVertex4);
     WWTraining23.doTraining( 275, 500, PrimaryVertex1, PrimaryVertex2);
     WWTraining23.doTraining( 275, 500, PrimaryVertex2, PrimaryVertex3);
     WWTraining23.doTraining( 275, 500, PrimaryVertex3, PrimaryVertex4);

     WWTraining24.doTraining( 200, 275, PrimaryVertex1, PrimaryVertex2);
     WWTraining24.doTraining( 200, 275, PrimaryVertex2, PrimaryVertex3);
     WWTraining24.doTraining( 200, 275, PrimaryVertex3, PrimaryVertex4);
     WWTraining24.doTraining( 275, 500, PrimaryVertex1, PrimaryVertex2);
     WWTraining24.doTraining( 275, 500, PrimaryVertex2, PrimaryVertex3);
     WWTraining24.doTraining( 275, 500, PrimaryVertex3, PrimaryVertex4);

     WWTraining2.doTraining( 200, 275 );
     WWTraining2.doTraining( 275, 500 );
        
     WWTraining3.doTraining( 200, 275, 0, 100);
     WWTraining3.doTraining( 275, 500, 0, 100);

     WWTraining4.doTraining( 200, 275 );
     WWTraining4.doTraining( 275, 500 );
     
     WWTraining5.doTraining( 200, 275 );
     WWTraining5.doTraining( 275, 500 );

     WWTraining6.doTraining( 200, 275 );
     WWTraining6.doTraining( 275, 500 );

     WWTraining7.doTraining( 200, 275 );
     WWTraining7.doTraining( 275, 500 );

     WWTraining8.doTraining( 200, 275 );
     WWTraining8.doTraining( 275, 500 );

     WWTraining9.doTraining( 200, 275 );
     WWTraining9.doTraining( 275, 500 );

     WWTraining10.doTraining( 200, 275 );
     WWTraining10.doTraining( 275, 500 );

     WWTraining11.doTraining( 200, 275 );
     WWTraining11.doTraining( 275, 500 );

     WWTraining12.doTraining( 200, 275 );
     WWTraining12.doTraining( 275, 500 );

     WWTraining13.doTraining( 200, 275 );
     WWTraining13.doTraining( 275, 500 );

     WWTraining14.doTraining( 200, 275 );
     WWTraining14.doTraining( 275, 500 );

     WWTraining15.doTraining( 200, 275 );
     WWTraining15.doTraining( 275, 500 );

     WWTraining16.doTraining( 200, 275 );
     WWTraining16.doTraining( 275, 500 );

    if options.makeFinalTree:

        # ---------------------------------------------------   
        # make the final trees
        print "making final trees"
        if not os.path.isdir("finalPlot"): os.system("mkdir finalPlot");
        if options.noX: gROOT.SetBatch()
        
        # get name from training class (if name is MassDrop, special case)
        # give it the same set of variables
        # give it samples
        Signalefficiency = 0.6
        graphindex = options.graphindex
#        
        rocs18ptv1 = WWTraining18.makeFinalPlotsInternal( 200, 275, PrimaryVertex1, PrimaryVertex2);
        rocs18ptv2 = WWTraining18.makeFinalPlotsInternal( 200, 275, PrimaryVertex2, PrimaryVertex3);
        #rocs18ptv3 = WWTraining18.makeFinalPlotsInternal( 200, 275, PrimaryVertex3, PrimaryVertex4);
        PrintTable(rocs18ptv1[graphindex], rocs18ptv2[graphindex], WWTraining18, 200, 275, PrimaryVertex1, PrimaryVertex2, PrimaryVertex3, Signalefficiency)
        rocs18ptv4 = WWTraining18.makeFinalPlotsInternal( 275, 500, PrimaryVertex1, PrimaryVertex2);
        rocs18ptv5 = WWTraining18.makeFinalPlotsInternal( 275, 500, PrimaryVertex2, PrimaryVertex3);
        #rocs18ptv6 = WWTraining18.makeFinalPlotsInternal( 275, 500, PrimaryVertex3, PrimaryVertex4);
        PrintTable(rocs18ptv4[graphindex], rocs18ptv5[graphindex], WWTraining18, 275, 500, PrimaryVertex1, PrimaryVertex2, PrimaryVertex3, Signalefficiency)
#        
        rocs22ptv1 = WWTraining22.makeFinalPlotsInternal( 200, 275, PrimaryVertex1, PrimaryVertex2);
        rocs22ptv2 = WWTraining22.makeFinalPlotsInternal( 200, 275, PrimaryVertex2, PrimaryVertex3);
        #rocs22ptv3 = WWTraining22.makeFinalPlotsInternal( 200, 275, PrimaryVertex3, PrimaryVertex4);
        PrintTable(rocs22ptv1[graphindex], rocs22ptv2[graphindex], WWTraining22, 200, 275, PrimaryVertex1, PrimaryVertex2, PrimaryVertex3, Signalefficiency)
        rocs22ptv4 = WWTraining22.makeFinalPlotsInternal( 275, 500, PrimaryVertex1, PrimaryVertex2);
        rocs22ptv5 = WWTraining22.makeFinalPlotsInternal( 275, 500, PrimaryVertex2, PrimaryVertex3);
        #rocs22ptv6 = WWTraining22.makeFinalPlotsInternal( 275, 500, PrimaryVertex3, PrimaryVertex4);
        PrintTable(rocs22ptv4[graphindex], rocs22ptv5[graphindex], WWTraining22, 275, 500, PrimaryVertex1, PrimaryVertex2, PrimaryVertex3, Signalefficiency)
        
        rocs23ptv1 = WWTraining23.makeFinalPlotsInternal( 200, 275, PrimaryVertex1, PrimaryVertex2);
        rocs23ptv2 = WWTraining23.makeFinalPlotsInternal( 200, 275, PrimaryVertex2, PrimaryVertex3);
        #rocs23ptv3 = WWTraining23.makeFinalPlotsInternal( 200, 275, PrimaryVertex3, PrimaryVertex4);
        PrintTable(rocs23ptv1[graphindex], rocs23ptv2[graphindex], WWTraining23, 200, 275, PrimaryVertex1, PrimaryVertex2, PrimaryVertex3, Signalefficiency)
        rocs23ptv4 = WWTraining23.makeFinalPlotsInternal( 275, 500, PrimaryVertex1, PrimaryVertex2);
        rocs23ptv5 = WWTraining23.makeFinalPlotsInternal( 275, 500, PrimaryVertex2, PrimaryVertex3);
        #rocs23ptv6 = WWTraining23.makeFinalPlotsInternal( 275, 500, PrimaryVertex3, PrimaryVertex4);
        PrintTable(rocs23ptv4[graphindex], rocs23ptv5[graphindex], WWTraining23, 275, 500, PrimaryVertex1, PrimaryVertex2, PrimaryVertex3, Signalefficiency)
        

        rocs1 = WWTraining1.makeFinalPlotsInternal( 200, 275 );
        rocs2 = WWTraining2.makeFinalPlotsInternal( 200, 275 );
        rocs3 = WWTraining3.makeFinalPlotsInternal( 200, 275 );
        rocs4 = WWTraining4.makeFinalPlotsInternal( 200, 275 );

        rocs1ex = WWTraining1.makeFinalPlots( 200, 275, 0. );
        rocs2ex = WWTraining2.makeFinalPlots( 200, 275, 0. );
        rocs3ex = WWTraining3.makeFinalPlots( 200, 275, 0. );        
        rocs4ex = WWTraining4.makeFinalPlots( 200, 275, 0. );
        
        canBdtRoc = ROOT.TCanvas("canBdtRoc","canBdtRoc",800,800);    
        canBdtRoc.cd();
        hrl = canBdtRoc.DrawFrame(0,0,1.0,1.0);
        hrl.GetXaxis().SetTitle("#epsilon_{sig}");
        hrl.GetYaxis().SetTitle("1 - #epsilon_{bkg}");
        canBdtRoc.SetGrid();
        rocs1[0].Draw();
        rocs1[1].SetLineColor(ROOT.kRed);
        rocs1[1].Draw();
        rocs2[0].SetLineColor(ROOT.kBlue);
        rocs2[0].Draw();
        rocs3[0].SetLineColor(ROOT.kMagenta);
        rocs3[0].Draw();
        rocs4[0].SetLineColor(ROOT.kGreen+2);
        rocs4[0].Draw();

        rocs1[2].SetLineWidth(2);
        rocs1[2].SetLineColor(ROOT.kBlack);
        rocs1[2].Draw();
        rocs2[2].SetLineWidth(2);
        rocs2[2].SetLineColor(ROOT.kBlue);
        rocs2[2].Draw();
        rocs3[2].SetLineWidth(2);
        rocs3[2].SetLineColor(ROOT.kMagenta);
        rocs3[2].Draw();
        rocs4[2].SetLineWidth(2);
        rocs4[2].SetLineColor(ROOT.kGreen+2);
        rocs4[2].Draw();
        
        rocs4[2].SetLineColor(ROOT.kCyan+2);
        rocs4[2].SetLineStyle(2);
        rocs4[2].Draw();
        
        leg = ROOT.TLegend(0.25,0.2,0.55,0.5)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)                
        leg.AddEntry( rocs1[0], "simple", 'l' );
        leg.AddEntry( rocs1[1], "mass drop only", 'l' ); 
        leg.AddEntry( rocs2[0], "optimal", 'l' );
        leg.AddEntry( rocs3[0], "optimal2", 'l' ); 
        leg.AddEntry( rocs4[0], "all", 'l' );         
        leg.AddEntry( rocs4ex[0], "all ex", 'l' );         
        leg.Draw();
        
        canBdtRoc.SaveAs("finalPlot/testROC_compall.eps");
        canBdtRoc.SaveAs("finalPlot/testROC_compall.png");

        hs = rocs2ex[2];
        hs.SetLineColor(4);        
        hs.Scale(1./hs.Integral());                
        hb = rocs2ex[3];
        hb.SetLineColor(2);        
        hb.Scale(1./hb.Integral());                
        canMassPrBdtCut = ROOT.TCanvas("canMassPrBdtCut","canMassPrBdtCut",800,800);    
        hs.Draw("hist");
        hb.Draw("histsames");
        canMassPrBdtCut.SaveAs("finalPlot/testmasspr.eps");
        canMassPrBdtCut.SaveAs("finalPlot/testmasspr.png");
        raw_input( 'Press ENTER to continue\n ' )
