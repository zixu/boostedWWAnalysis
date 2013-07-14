# A class which takes histograms and plots them in a versatile way
# inputs are file names which can be "data" or "MC"

import ROOT
#ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C")
#from ROOT import setTDRStyle
#from ROOT import TTree
#ROOT.setTDRStyle()
#ROOT.gStyle.SetPalette(1)
#ROOT.gStyle.SetOptFit(0)
#
#ROOT.TH1.SetDefaultSumw2()
#ROOT.TH2.SetDefaultSumw2()
#
#ROOT.gStyle.SetPadTopMargin(0.09);
#ROOT.gStyle.SetPadLeftMargin(0.16);

    #    TMVA::Tools::Instance();
ROOT.TMVA.Tools.Instance();    

import os

from array import array
import math
from optparse import OptionParser

def getfractionbelowcut( cut, list ):
    
    ctr = 0;
    for i in range(len(list)): 
        if list[i] < cut: ctr = ctr + 1;
    return (float(ctr)/float(len(list)));
        

def ComputeROCFromList(ls,lb,LtoR):

    lsmax = max(ls);
    lbmax = max(lb);
    lsmin = min(ls);
    lbmin = min(lb);

    allmax = max(lsmax,lbmax);
    allmin = max(lsmin,lbmin);

    xval = array('d', [])
    yval = array('d', [])
    
    nsteps = 1000;
    stepsize = (allmax - allmin)/nsteps;
    for i in range(nsteps+1):
        curCutVal = stepsize*float(i) + allmin;
        effsig = getfractionbelowcut( curCutVal, ls );
        effbkg = getfractionbelowcut( curCutVal, lb );
        
        if LtoR:
            xval.append( 1-effsig );
            yval.append( effbkg );
        else:
            xval.append( effsig );
            yval.append( 1-effbkg );

    tg = ROOT.TGraph( nsteps+1, xval, yval );
    return tg;

def ComputeROC(hsig,hbkg,LtoR):
    
    nbins = hsig.GetNbinsX();
    binsize = hsig.GetBinWidth(1);
    lowedge = hsig.GetBinLowEdge(1);
    
    print "lowedge: ",lowedge
    
    hsigIntegral = hsig.Integral();
    hbkgIntegral = hbkg.Integral();
    
    xval = array('d', [])
    yval = array('d', [])
    ctr = 0;
    for i in range(1,nbins+1):
        
        effBkg = 0;
        effSig = 0;
        
        if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
        else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;
        
        if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
        else: effSig = hsig.Integral( 1, i )/hsigIntegral;
        
        print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", effBkg, ", effSig: ", effSig;
        
        xval.append( effSig );
        yval.append( 1-effBkg );
        ctr = ctr + 1;
    
    print nbins, "and ", ctr
    tg = ROOT.TGraph( nbins, xval, yval );
    return tg;


class trainingClass:

    ### ------------------------------------------------
    def __init__(self, signalFile, backgroundFile, listOfTrainingVariables, label):

        print "Welcome to the training..."
        self.SigFile_ = ROOT.TFile(signalFile);
        self.BkgFile_ = ROOT.TFile(backgroundFile);
        
        self.SigTree_ = self.SigFile_.Get("otree");
        self.BkgTree_ = self.BkgFile_.Get("otree");
    
        self.ListOfTrainingVariables = listOfTrainingVariables;
    
        self.Label_ = label;
        
    ############################
    ############################
    ############################

    def doTraining( self, pTlo, pThi ):
        
        print "pT range: ", pTlo, " - ", pThi;

        self.OutputFileName_ = "classifier/Wtagger_"+str(pTlo)+"to"+str(pThi)+"_"+self.Label_+".root";
        outputFile = ROOT.TFile( self.OutputFileName_, 'RECREATE' )
        factory = ROOT.TMVA.Factory( "Wtagger_"+str(pTlo)+"to"+str(pThi)+"_"+self.Label_, outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )
                
        #training variables
        for i in range(len(self.ListOfTrainingVariables)):
            varString = self.ListOfTrainingVariables[i]+" := "+self.ListOfTrainingVariables[i];
            print varString
            factory.AddVariable( varString, 'F' );
        
        #spectators
        factory.AddSpectator( "jet_pt_pr", 'F' )
        factory.AddSpectator( "jet_mass_pr", 'F' )
        
        # Global event weights (see below for setting event-wise weights)
        signalWeight     = 1.0
        backgroundWeight = 1.0
        
        factory.AddSignalTree    ( self.SigTree_,     signalWeight     )
        factory.AddBackgroundTree( self.BkgTree_,     backgroundWeight )
        
        # what's this for?
        #factory.SetBackgroundWeightExpression( "weight" )
        
        # cuts definition
        mycutSig = ROOT.TCut( "jet_pt_pr < "+str(pThi)+" && jet_pt_pr > "+str(pTlo)+" && jet_mass_pr > 60 && jet_mass_pr < 100" );
        mycutBkg = ROOT.TCut( "jet_pt_pr < "+str(pThi)+" && jet_pt_pr > "+str(pTlo)+" && jet_mass_pr > 60 && jet_mass_pr < 100" );
        
        # Here, the relevant variables are copied over in new, slim trees that are
        # used for TMVA training and testing
        # "SplitMode=Random" means that the input events are randomly shuffled before
        # splitting them into training and test samples
        print "PrepareTrainingAndTestTree ... "
        factory.PrepareTrainingAndTestTree( mycutSig, mycutBkg,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )
        print "BookMethod ... "
#        factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT","!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" )
        factory.BookMethod( ROOT.TMVA.Types.kBDT, "BDT","!H:!V:NTrees=1000:nEventsMin=40:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=4" )
        
        # Train MVAs
        print "TrainAllMethods ... "        
        factory.TrainAllMethods()
        
        # Test MVAs
        print "TestAllMethods ... "        
        factory.TestAllMethods()
        
        # Evaluate MVAs
        print "EvaluateAllMethods ... "        
        factory.EvaluateAllMethods()    
        
        # Save the output.
        outputFile.Close()

    ############################
    ############################
    ############################

    def plotTrainingResults( self, pTlo, pThi ):

        rootname = "Wtagger_"+str(pTlo)+"to"+str(pThi)+"_"+self.Label_;
        
        if not os.path.isdir("plots"): os.system("mkdir plots");
        if not os.path.isdir("plots/classifier"): os.system("mkdir plots/classifier");
        os.system("bash TMVAscripts/runAllScripts.sh classifier/"+rootname+".root");

        if not os.path.isdir("classifier/trainingPlots"): os.system("mkdir classifier/trainingPlots");        
        if not os.path.isdir("classifier/trainingPlots/"+rootname): os.system("mkdir classifier/trainingPlots/"+rootname);        
        os.system("mv plots/classifier/"+rootname+".root* classifier/trainingPlots/"+rootname+"/.");
        
    def getWeightFile( self, pTlo, pThi ):
        self.WeightFile_ = "weights/"+"Wtagger_"+str(pTlo)+"to"+str(pThi)+"_"+self.Label_+"_BDT.weights.xml";
        return self.WeightFile_
    
    ############################
    ############################
    ############################
    def makeFinalPlots( self, pTlo, pThi, bdtcut ):

        print self.ListOfTrainingVariables;
        print self.SigTree_.GetEntries();
        print self.BkgTree_.GetEntries();
                
#        h_sig = [];
#        h_bkg = [];
        
        nameOfWeightFile = self.getWeightFile( pTlo, pThi );
        
        # here we have a flag to either use TMVA file or original files
        
        # make a list of arrays
        reader = ROOT.TMVA.Reader("!Color:!Silent")
        listOfVarArray = [];
        for i in range(len(self.ListOfTrainingVariables)):
#            curVar = array('f',[0.]);
            listOfVarArray.append( array('f',[0.]) );
            varString = self.ListOfTrainingVariables[i]+" := "+self.ListOfTrainingVariables[i];
            reader.AddVariable( varString, listOfVarArray[i] );
        
        #spectators
        spec1 = array('f',[0.]);
        spec2 = array('f',[0.]);
        reader.AddSpectator( "jet_pt_pr", spec1 )
        reader.AddSpectator( "jet_mass_pr", spec2 )

        reader.BookMVA("BDT",nameOfWeightFile);

        mdVal_sig = []; mdVal_bkg = [];
        discrVal_sig = []; discrVal_bkg = [];
        
        h_mass_pr_sig = ROOT.TH1F("h_mass_pr_sig","; pruned mass;", 50, 0, 150 );
        h_mass_pr_bkg = ROOT.TH1F("h_mass_pr_bkg","; pruned mass;", 50, 0, 150 );

        ## ---------------
        for i in range(self.SigTree_.GetEntries()): 
            
            self.SigTree_.GetEntry(i);
            
            for j in range(len(self.ListOfTrainingVariables)):
                listOfVarArray[j][0] = getattr( self.SigTree_,self.ListOfTrainingVariables[j] );
            
            bdtv = reader.EvaluateMVA("BDT");
            if bdtv > bdtcut: h_mass_pr_sig.Fill( getattr( self.SigTree_, "jet_mass_pr" ) );

#            print "mass = ", getattr( self.SigTree_, "jet_mass_pr" ), " and discriminant = ", bdtv

            if getattr(self.SigTree_,"jet_pt_pr") > pTlo and getattr(self.SigTree_,"jet_pt_pr") < pThi and getattr(self.SigTree_,"jet_mass_pr") > 60. and getattr(self.SigTree_,"jet_mass_pr") < 100:
                
                mdVal_sig.append( getattr( self.SigTree_, "jet_massdrop_pr") );
                discrVal_sig.append( bdtv );
            

        ## ---------------
        for i in range(self.BkgTree_.GetEntries()): 
            
            self.BkgTree_.GetEntry(i);
            
            for j in range(len(self.ListOfTrainingVariables)):
                listOfVarArray[j][0] = getattr( self.BkgTree_,self.ListOfTrainingVariables[j] );
            
            bdtv = reader.EvaluateMVA("BDT");
            if bdtv > bdtcut: h_mass_pr_bkg.Fill( getattr( self.BkgTree_, "jet_mass_pr" ) );                                                      

            if getattr(self.BkgTree_,"jet_pt_pr") > pTlo and getattr(self.BkgTree_,"jet_pt_pr") < pThi and getattr(self.BkgTree_,"jet_mass_pr") > 60. and getattr(self.BkgTree_,"jet_mass_pr") < 100:

                mdVal_bkg.append( getattr( self.BkgTree_, "jet_massdrop_pr") );
                discrVal_bkg.append( reader.EvaluateMVA("BDT") );
                                                      
        ## ---------------
                                                      
        hmd_sig = ROOT.TH1F("hmd_sig","hmd_sig",100,min(mdVal_sig),max(mdVal_sig));
        for i in range(len(mdVal_sig)): hmd_sig.Fill( mdVal_sig[i] );
        hmd_bkg = ROOT.TH1F("hmd_bkg","hmd_bkg",100,min(mdVal_bkg),max(mdVal_bkg)); 
        for i in range(len(mdVal_bkg)): hmd_bkg.Fill( mdVal_bkg[i] );

        hdiscr_sig = ROOT.TH1F("hdiscr_sig","hdiscr_sig",100,min(discrVal_sig),max(discrVal_sig));
        for i in range(len(discrVal_sig)): hdiscr_sig.Fill( discrVal_sig[i] );
        hdiscr_bkg = ROOT.TH1F("hdiscr_bkg","hdiscr_bkg",100,min(discrVal_bkg),max(discrVal_bkg));
        for i in range(len(discrVal_bkg)): hdiscr_bkg.Fill( discrVal_bkg[i] );

        can1ex = ROOT.TCanvas("can1ex","can1ex",800,800);
        can1ex.cd();
        hmd_sig.Draw("hist");
        hmd_bkg.SetLineColor(2);
        hmd_bkg.Draw("histsames");
        can1ex.SaveAs("finalPlot/extestmd_"+self.Label_+".eps");
        can1ex.SaveAs("finalPlot/extestmd_"+self.Label_+".png");

        can2ex = ROOT.TCanvas("can2ex","can2ex",800,800);
        can2ex.cd();
        hdiscr_sig.Draw("hist");
        hdiscr_bkg.SetLineColor(2);
        hdiscr_bkg.Draw("histsames");
        can2ex.SaveAs("finalPlot/extestdiscr_"+self.Label_+".eps");
        can2ex.SaveAs("finalPlot/extestdiscr_"+self.Label_+".png");

        tgs = [];
        tgs.append( ComputeROCFromList(discrVal_sig,discrVal_bkg, True) );
        tgs.append( ComputeROCFromList(mdVal_sig,mdVal_bkg, False) );
        tgs.append( h_mass_pr_sig );
        tgs.append( h_mass_pr_bkg );
        return tgs;

    ############################
    ############################
    ############################
    def makeFinalPlotsInternal(self, pTlo, pThi):

        self.OutputFileName_ = "classifier/Wtagger_"+str(pTlo)+"to"+str(pThi)+"_"+self.Label_+".root";

        internalTreeName = ROOT.TFile(self.OutputFileName_);
        internalTree = internalTreeName.Get("TestTree");
        internalTrainTree = internalTreeName.Get("TrainTree");
        
        mdVal_sig = []; mdVal_bkg = [];
        discrVal_sig = []; discrVal_bkg = [];
        discrVal_sig_train = []; discrVal_bkg_train = [];        
        ## ---------------
        for i in range(internalTree.GetEntries()): 
            
            internalTree.GetEntry(i);
            
            if getattr(internalTree,"jet_pt_pr") > pTlo and getattr(internalTree,"jet_pt_pr") < pThi and getattr(internalTree,"jet_mass_pr") > 60. and getattr(internalTree,"jet_mass_pr") < 100 and internalTree.classID == 0:
                                    
                mdVal_sig.append( getattr( internalTree, "jet_massdrop_pr") );
                discrVal_sig.append( getattr( internalTree, "BDT") );

        ## ---------------
        for i in range(internalTree.GetEntries()): 
            
            internalTree.GetEntry(i);
            
            if getattr(internalTree,"jet_pt_pr") > pTlo and getattr(internalTree,"jet_pt_pr") < pThi and getattr(internalTree,"jet_mass_pr") > 60. and getattr(internalTree,"jet_mass_pr") < 100 and internalTree.classID == 1:
                
                mdVal_bkg.append( getattr( internalTree, "jet_massdrop_pr") );
                discrVal_bkg.append( getattr( internalTree, "BDT") );
        ## ---------------
        for i in range(internalTrainTree.GetEntries()): 
            
            internalTrainTree.GetEntry(i);
            
            if getattr(internalTrainTree,"jet_pt_pr") > pTlo and getattr(internalTrainTree,"jet_pt_pr") < pThi and getattr(internalTrainTree,"jet_mass_pr") > 60. and getattr(internalTrainTree,"jet_mass_pr") < 100 and internalTrainTree.classID == 0:
                
                mdVal_sig.append( getattr( internalTrainTree, "jet_massdrop_pr") );
                discrVal_sig_train.append( getattr( internalTrainTree, "BDT") );
        
        ## ---------------
        for i in range(internalTrainTree.GetEntries()): 
            
            internalTrainTree.GetEntry(i);
            
            if getattr(internalTrainTree,"jet_pt_pr") > pTlo and getattr(internalTrainTree,"jet_pt_pr") < pThi and getattr(internalTrainTree,"jet_mass_pr") > 60. and getattr(internalTrainTree,"jet_mass_pr") < 100 and internalTrainTree.classID == 1:
                
                mdVal_bkg.append( getattr( internalTrainTree, "jet_massdrop_pr") );
                discrVal_bkg_train.append( getattr( internalTrainTree, "BDT") );

            
                    
        hmd_sig = ROOT.TH1F("hmd_sig",";mass drop ("+self.Label_+");",30,min(mdVal_sig),max(mdVal_sig));
        for i in range(len(mdVal_sig)): hmd_sig.Fill( mdVal_sig[i] );
        hmd_bkg = ROOT.TH1F("hmd_bkg",";mass drop ("+self.Label_+");",30,min(mdVal_bkg),max(mdVal_bkg)); 
        for i in range(len(mdVal_bkg)): hmd_bkg.Fill( mdVal_bkg[i] );
        
        hdiscr_sig = ROOT.TH1F("hdiscr_sig",";BDT discr ("+self.Label_+");",30,min(discrVal_sig),max(discrVal_sig));
        for i in range(len(discrVal_sig)): hdiscr_sig.Fill( discrVal_sig[i] );
        hdiscr_bkg = ROOT.TH1F("hdiscr_bkg",";BDT discr ("+self.Label_+");",30,min(discrVal_bkg),max(discrVal_bkg));
        for i in range(len(discrVal_bkg)): hdiscr_bkg.Fill( discrVal_bkg[i] );

        hdiscr_sig_train = ROOT.TH1F("hdiscr_sig_train",";BDT discr ("+self.Label_+");",30,min(discrVal_sig_train),max(discrVal_sig_train));
        for i in range(len(discrVal_sig_train)): hdiscr_sig_train.Fill( discrVal_sig_train[i] );
        hdiscr_bkg_train = ROOT.TH1F("hdiscr_bkg_train",";BDT discr ("+self.Label_+");",30,min(discrVal_bkg_train),max(discrVal_bkg_train));
        for i in range(len(discrVal_bkg_train)): hdiscr_bkg_train.Fill( discrVal_bkg_train[i] );
                    
                    
        can1 = ROOT.TCanvas("can1"+self.Label_,"can1"+self.Label_,800,800);
        can1.cd();
        hmd_sig.Draw("hist");
        hmd_bkg.SetLineColor(2);
        hmd_bkg.Draw("histsames");
        can1.SaveAs("finalPlot/testmd_"+self.Label_+".eps");
        can1.SaveAs("finalPlot/testmd_"+self.Label_+".png");
        
        can2 = ROOT.TCanvas("can2"+self.Label_,"can2"+self.Label_,800,800);
        can2.cd();
        hdiscr_sig.Draw("hist");
        hdiscr_bkg.SetLineColor(2);
        hdiscr_bkg.Draw("histsames");
                    
        hdiscr_sig_train.SetLineStyle(2);
        hdiscr_sig_train.Draw("histsames");
                    
        hdiscr_bkg_train.SetLineColor(2);
        hdiscr_bkg_train.SetLineStyle(2);
                    
        hdiscr_bkg_train.Draw("histsames");
        can2.SaveAs("finalPlot/testdiscr_"+self.Label_+".eps");
        can2.SaveAs("finalPlot/testdiscr_"+self.Label_+".png");
        
        tgs = [];
#        tgs.append( ComputeROC(hdiscr_sig,hdiscr_bkg, True) );
#        tgs.append( ComputeROC(hmd_sig,hmd_bkg, False) );
        tgs.append( ComputeROCFromList(discrVal_sig,discrVal_bkg, True) );
        tgs.append( ComputeROCFromList(mdVal_sig,mdVal_bkg, False) );
        tgs.append( ComputeROCFromList(discrVal_sig_train,discrVal_bkg_train, True) );

        return tgs;
