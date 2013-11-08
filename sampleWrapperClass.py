# A class which takes histograms and plots them in a versatile way
# inputs are file names which can be "data" or "MC"

import ROOT
import os
import sys
import math

from array import array
from optparse import OptionParser

from ROOT import RooTrace
from ROOT import TTree
from ROOT import gROOT

from bsmReweighter import *

import sys
from DataFormats.FWLite import Events, Handle

ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(0)

ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);


# FWLITE stuff
ROOT.gSystem.Load('libCondFormatsJetMETObjects')
ROOT.gSystem.Load('libFWCoreFWLite');
ROOT.gSystem.Load('libFWCoreUtilities');  

### ------------ h e l p e r s --------------------

def getListRMS(list):
    mean = sum(list)/float(len(list));
    return math.sqrt(sum((n-mean)*(n-mean) for n in list)/len(list));
def getListMean(list):
    return sum(list)/float(len(list));
### ----------- class implementation -----------

class sampleWrapperClass:

    ### ------------------------------------------------
    def __init__(self, label, file, channel, sampleEffLumi, lumi, treename, isData,outputfiledirectory):

        
        self.IsData_ = isData; ### flag when running on data
        self.FileName_ = file; ### input file name 
        self.File_ = ROOT.TFile(file); ## get the root file
        self.InputTree_ = self.File_.Get(treename); ## get the input tree
        
        self.SampleWeight_ = lumi/sampleEffLumi; ## take the lumi re-weight factor
        
        self.JetPrefix_ = "GroomedJet_CA8";
        self.Label_ = label;
        self.Channel_ = channel
        self.OFileName_ = outputfiledirectory+"trainingtrees_"+channel+"/ofile_"+label+".root";
    
        # initialization for doing BSM reweighting
        self.SignalMass_ = -1;
        
        if file.find("HWW") > 0 and file.find("600") > 0: self.SignalMass_ = 600;
        if file.find("HWW") > 0 and file.find("700") > 0: self.SignalMass_ = 700;
        if file.find("HWW") > 0 and file.find("800") > 0: self.SignalMass_ = 800;
        if file.find("HWW") > 0 and file.find("900") > 0: self.SignalMass_ = 900;    
        if file.find("HWW") > 0 and file.find("1000") > 0: self.SignalMass_ = 1000;   

        self.FitSMSignal = False;
        self.FitSMSignal_mean = -1;
        self.FitSMSignal_gamma = -1;
        self.isVBF_ = False;

        if file.find("VBF") > 0: self.isVBF_ = True;

        ###### reweight file
        if self.SignalMass_ > 0:
            self.rwName = "H_CPSrw_%03d"%(self.SignalMass_);
            self.rwF = ROOT.TFile("CPSrw/"+self.rwName+".root");
            self.h_rwCPS = self.rwF.Get(self.rwName);
            self.x_rwCPS = self.h_rwCPS.GetXaxis();

        # ---------- Set up jet corrections on the fly of R >= 0.7 jets
        fDir = "JECs/"      
        jecUncStr       = ROOT.std.string(fDir + "GR_R_53_V10_Uncertainty_AK7PFchs.txt")
        self.jecUnc_    = ROOT.JetCorrectionUncertainty(jecUncStr)
        jecUncStrAK5    = ROOT.std.string(fDir + "GR_R_53_V10_Uncertainty_AK5PFchs.txt")
        self.jecUncAK5_ = ROOT.JetCorrectionUncertainty(jecUncStr)

            
    def getLabel(self):
        return self.Label_

    def getTrainingTreeName(self):
        return self.OFileName_
    
    def getSampleLabel(self):
        return self.Label_


    def createTrainingTree(self):
        
        print self.FileName_
        self.NTree_ = self.InputTree_.GetEntries();
        
        print "Turning off branches...", self.FileName_

        self.InputTree_.SetBranchStatus("vbf*Charged*",0); 
        self.InputTree_.SetBranchStatus("vbf*Neutral*",0); 
        self.InputTree_.SetBranchStatus("vbf*Neutral*",0); 
        self.InputTree_.SetBranchStatus("vbf*Photon*",0); 
        self.InputTree_.SetBranchStatus("vbf*Electron*",0); 
        self.InputTree_.SetBranchStatus("vbf*HF*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Charged*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Neutral*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Neutral*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Electron*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*HF*",0); 
        self.InputTree_.SetBranchStatus("Hadronic*",0); 
        self.InputTree_.SetBranchStatus("boosted*ang*",0); 
        self.InputTree_.SetBranchStatus("fit_*",0); 
        self.InputTree_.SetBranchStatus("W_tb*",0); 
        self.InputTree_.SetBranchStatus("W_Parton*",0); 
        self.InputTree_.SetBranchStatus("JetGen*",0); 

        print "Initializing sample: ", self.FileName_
        print "Nentries = ", self.NTree_

        # fill histograms
        self.createBRDTree();
        self.File_.Close();


    ### produce weights for alternative models
    def GetInteferenceWeights(self, mass, Cprime, BRnew=0):
        
        massin = self.SignalMass_;
        if massin < 0: return;
        
        massmin = {600:200,700:200,800:400,900:400,1000:400};
        massmax = {600:1200,700:1200,800:1500,900:1600,1000:1800};

        # read in original file and get lineshape from fit
        if not self.FitSMSignal: 
            fitSM = FitMassPoint(self.FileName_, massin, massmin[massin], massmax[massin]);
            self.FitSMSignal_mean  = fitSM[0];
            self.FitSMSignal_gamma = fitSM[1];
            self.FitSMSignal = True;
    
        # create a weight for the given event after switch is turned
        weight_width = lineshapeWidthReweight(mass, self.FitSMSignal_mean, self.FitSMSignal_gamma, Cprime/(1.-BRnew), massmin[massin], massmax[massin]);
        weight_xs = Cprime*(1-BRnew);   
        return (weight_width*weight_xs);

    ### ------------------------------------------------
    def createBRDTree(self):

        reader1_ = ROOT.TMVA.Reader("!Color:!Silent")
        reader2_ = ROOT.TMVA.Reader("!Color:!Silent")
        ungroomed_reader1_ = ROOT.TMVA.Reader("!Color:!Silent")
        ungroomed_reader2_ = ROOT.TMVA.Reader("!Color:!Silent")

        fname = self.OFileName_;        
        self.OFile_ = ROOT.TFile(fname,"RECREATE");        
        self.otree = ROOT.TTree("otree","otree");

        ########## Initialize Variables
        self.InitializeVariables();

        ########## bsm branches for ewk signlet 

        self.cprimeVals   = [1,2,3,4,5,6,7,8,9,10] ;
        self.brnewVals    = [00, 01, 02, 03, 04, 05] ;
        self.bsmReweights = [];
        
        for i in range(len(self.cprimeVals)): 
          col_bsmReweights = [];
          for j in range(len(self.brnewVals)): 
            col_bsmReweights.append( array( 'f', [ 0. ] ) );
            self.bsmReweights.append( col_bsmReweights );

        ########## Create Output branches
        self.createBranches();
        
        ###### loop preparation
        prefix = self.JetPrefix_;        
        NLoop = min(self.NTree_,1e9);
        NLoopWeight = self.NTree_/NLoop;
        
        wSampleWeight = NLoopWeight*self.SampleWeight_;

        wjettaggervariables = []
        wjettaggervariables.append("jet_qjetvol")
        wjettaggervariables.append("jet_tau2tau1")
        listOfVarArray1 = []
        listOfVarArray2 = []
        for i in range(len(wjettaggervariables)):
            listOfVarArray1.append( array('f',[0.]) )
            listOfVarArray2.append( array('f',[0.]) )
            varString = wjettaggervariables[i] + " := " + wjettaggervariables[i]
            reader1_.AddVariable(varString, listOfVarArray1[i])
            reader2_.AddVariable(varString, listOfVarArray2[i])
            ungroomed_reader1_.AddVariable(varString, listOfVarArray1[i])
            ungroomed_reader2_.AddVariable(varString, listOfVarArray2[i])
        
        spec11 = array('f',[0.])
        spec12 = array('f',[0.])
        spec13 = array('f',[0.])
        spec14 = array('f',[0.])
        spec21 = array('f',[0.])
        spec22 = array('f',[0.])
        spec23 = array('f',[0.])
        spec24 = array('f',[0.])
         
        ungroomed_spec11 = array('f',[0.])
        ungroomed_spec12 = array('f',[0.])
        ungroomed_spec13 = array('f',[0.])
        ungroomed_spec14 = array('f',[0.])
        ungroomed_spec21 = array('f',[0.])
        ungroomed_spec22 = array('f',[0.])
        ungroomed_spec23 = array('f',[0.])
        ungroomed_spec24 = array('f',[0.])
        reader1_.AddSpectator( "jet_pt_pr", spec11 )
        reader1_.AddSpectator( "jet_mass_pr", spec12 )
        reader1_.AddSpectator( "nPV", spec13 )
        reader1_.AddSpectator( "jet_massdrop_pr", spec14 )
        reader2_.AddSpectator( "jet_pt_pr", spec21 )
        reader2_.AddSpectator( "jet_mass_pr", spec22 )
        reader2_.AddSpectator( "nPV", spec23 )
        reader2_.AddSpectator( "jet_massdrop_pr", spec24 )

        ungroomed_reader1_.AddSpectator( "ungroomed_jet_pt", ungroomed_spec11 )
        ungroomed_reader1_.AddSpectator( "jet_mass_pr", ungroomed_spec12 )
        ungroomed_reader1_.AddSpectator( "nPV", ungroomed_spec13 )
        ungroomed_reader1_.AddSpectator( "jet_massdrop_pr", ungroomed_spec14 )
        ungroomed_reader2_.AddSpectator( "ungroomed_jet_pt", ungroomed_spec21 )
        ungroomed_reader2_.AddSpectator( "jet_mass_pr", ungroomed_spec22 )
        ungroomed_reader2_.AddSpectator( "nPV", ungroomed_spec23 )
        ungroomed_reader2_.AddSpectator( "jet_massdrop_pr", ungroomed_spec24 )

        TMVAMethod = "Likelihood" 

        RooTrace.active(ROOT.kTRUE);
        RooTrace.mark();

        ################################################
        ########## Start Loop On the Evvents ###########
        ################################################
        
        for i in range(NLoop) :
            
            if i % 100000 == 0: print "iEvent = ", i

            self.InputTree_.GetEntry(i);

            ########## Fill TTbar Control Region Info

            ttbarlike  = 0;

            if self.Channel_ == 'mu': lepLabel = "muon";
            if self.Channel_ == 'el': lepLabel = "electron";

            index_ca8_in_oppoHemi = [];
            
            for i in range(6):
                if getattr( self.InputTree_, "GroomedJet_CA8_pt" )[i] > 200 and math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_eta" )[i])<2.4:
                    j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[i];
                    j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[i];
                    l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                    l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );
                    l_charge = getattr( self.InputTree_, "W_"+lepLabel+"_charge" );
                    dR_lj = math.sqrt( (l_eta - j_ca8_eta)**2 + (l_phi - j_ca8_phi)**2 ); ## delta R between lepton and hadronic W candidate over 200 GeV
                    if dR_lj > ROOT.TMath.Pi()/2.: index_ca8_in_oppoHemi.append(i); ## opposite hemishpere

            minMass = -1;
            theca8Index = -1;
            for i in range(len(index_ca8_in_oppoHemi)):
                curmass = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[index_ca8_in_oppoHemi[i]]; ## search for the jet with mass closer to the W
                if math.fabs(curmass-80.385) < math.fabs(minMass-80.385): 
                    minMass = curmass;
                    theca8Index = index_ca8_in_oppoHemi[i];

            index_ak5_in_sameHemi = [];
            index_ak5_in_oppoHemi = [];
            index_ak5_in_sameHemi_vetoca8 = [];
            index_ak5_in_oppoHemi_vetoca8 = [];
            
            index_ak5_in_sameHemi_csvl = [];
            index_ak5_in_oppoHemi_csvl = [];
            index_ak5_in_sameHemi_vetoca8_csvl = [];
            index_ak5_in_oppoHemi_vetoca8_csvl = [];
            
            index_ak5_in_sameHemi_csvm = [];
            index_ak5_in_oppoHemi_csvm = [];
            index_ak5_in_sameHemi_vetoca8_csvm = [];
            index_ak5_in_oppoHemi_vetoca8_csvm = [];
            
            index_ak5_in_sameHemi_csvt = [];
            index_ak5_in_oppoHemi_csvt = [];
            index_ak5_in_sameHemi_vetoca8_csvt = [];
            index_ak5_in_oppoHemi_vetoca8_csvt = [];

            ttb_ht = getattr( self.InputTree_, "W_"+lepLabel+"_pt" );
            ttb_ht += getattr( self.InputTree_, "event_met_pfmet" ); 
        
            if theca8Index >= 0: ## if a W candidate is found 
                for i in range(6):
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[i] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[i] > 30: ## loop on the ak5 with pt > 30 only central jets 
                        ttb_ht += getattr( self.InputTree_, "JetPFCor_Pt" )[i];
                        j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];
                        j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[i];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[i];
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );               
                        dR_jj = math.sqrt( (j_ak5_eta - j_ca8_eta)**2 + (j_ak5_phi - j_ca8_phi)**2 ); ## delta R jet jet 
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + (l_phi - j_ak5_phi)**2 ); ## delta R jet lep
                        
                        if dR_lj < ROOT.TMath.Pi()/2. :  ## same hemisphere wrt to jet and cleaned wrt to the fat jet
                            index_ak5_in_sameHemi.append( i );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_sameHemi_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_sameHemi_csvm.append(i);                        
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_sameHemi_csvt.append(i);                        
                        elif dR_lj > ROOT.TMath.Pi()/2. : 
                            index_ak5_in_oppoHemi.append( i );    
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_oppoHemi_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_oppoHemi_csvm.append(i);                                            
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_oppoHemi_csvt.append(i);                                            

                        if dR_lj > ROOT.TMath.Pi()/2. and dR_jj > 0.8: 
                            index_ak5_in_oppoHemi_vetoca8.append( i );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_oppoHemi_vetoca8_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_oppoHemi_vetoca8_csvm.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_oppoHemi_vetoca8_csvt.append(i);
                        elif dR_lj < ROOT.TMath.Pi()/2. and  dR_jj > 0.8:
                            index_ak5_in_sameHemi_vetoca8.append( i );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_sameHemi_vetoca8_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_sameHemi_vetoca8_csvm.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_sameHemi_vetoca8_csvt.append(i);


                for i in range(6):
                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[i] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[i] > 30: ## loop on the ak5 with pt > 30 only central jets 
                        ttb_ht += getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[i];
                        j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];
                        j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCorVBFTag_Eta" )[i];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCorVBFTag_Phi" )[i];
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );               
                        dR_jj = math.sqrt( (j_ak5_eta - j_ca8_eta)**2 + (j_ak5_phi - j_ca8_phi)**2 ); ## delta R jet jet 
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + (l_phi - j_ak5_phi)**2 ); ## delta R jet lep
                        
                        if dR_lj < ROOT.TMath.Pi()/2. :  ## same hemisphere wrt to jet and cleaned wrt to the fat jet
                            index_ak5_in_sameHemi.append( i );
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_sameHemi_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_sameHemi_csvm.append(i);                        
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_sameHemi_csvt.append(i);                        
                        elif dR_lj > ROOT.TMath.Pi()/2. : 
                            index_ak5_in_oppoHemi.append( i );    
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_oppoHemi_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_oppoHemi_csvm.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_oppoHemi_csvt.append(i);
                        if dR_lj > ROOT.TMath.Pi()/2. and dR_jj > 0.8: 
                            index_ak5_in_oppoHemi_vetoca8.append( i );
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_oppoHemi_vetoca8_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_oppoHemi_vetoca8_csvm.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_oppoHemi_vetoca8_csvt.append(i);
                        elif dR_lj < ROOT.TMath.Pi()/2. and  dR_jj > 0.8:
                            index_ak5_in_sameHemi_vetoca8.append( i );
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_sameHemi_vetoca8_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_sameHemi_vetoca8_csvm.append(i);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_sameHemi_vetoca8_csvt.append(i);

                ### number of jets in the same emisphere of the lepton 
                self.ttb_nak5_same_[0] = int(len(index_ak5_in_sameHemi));
                ### number of jets in the same emisphere of the lepton passing csvl                 
                self.ttb_nak5_same_csvl_[0] = int(len(index_ak5_in_sameHemi_csvl));
                ### number of jets in the same emisphere of the lepton passing csvm                                 
                self.ttb_nak5_same_csvm_[0] = int(len(index_ak5_in_sameHemi_csvm));
                ### number of jets in the same emisphere of the lepton passing csvt                                                 
                self.ttb_nak5_same_csvt_[0] = int(len(index_ak5_in_sameHemi_csvt));
                ### number of jets in the opposite emisphere of the lepton                 
                self.ttb_nak5_oppo_[0] = int(len(index_ak5_in_oppoHemi));
                self.ttb_nak5_oppo_csvl_[0] = int(len(index_ak5_in_oppoHemi_csvl));
                self.ttb_nak5_oppo_csvm_[0] = int(len(index_ak5_in_oppoHemi_csvm));
                self.ttb_nak5_oppo_csvt_[0] = int(len(index_ak5_in_oppoHemi_csvt));
                ### number of jets in the opposite emisphere of the lepton cleaned wrt to the hadronic candidate                                
                self.ttb_nak5_oppoveto_[0] = int(len(index_ak5_in_oppoHemi_vetoca8));
                self.ttb_nak5_oppoveto_csvl_[0] = int(len(index_ak5_in_oppoHemi_vetoca8_csvl));
                self.ttb_nak5_oppoveto_csvm_[0] = int(len(index_ak5_in_oppoHemi_vetoca8_csvm));
                self.ttb_nak5_oppoveto_csvt_[0] = int(len(index_ak5_in_oppoHemi_vetoca8_csvt));
                ### number of jets in the same emisphere of the lepton cleaned wrt to the hadronic candidate                                
                self.ttb_nak5_sameveto_[0] = int(len(index_ak5_in_sameHemi_vetoca8));
                self.ttb_nak5_sameveto_csvl_[0] = int(len(index_ak5_in_sameHemi_vetoca8_csvl));
                self.ttb_nak5_sameveto_csvm_[0] = int(len(index_ak5_in_sameHemi_vetoca8_csvm));
                self.ttb_nak5_sameveto_csvt_[0] = int(len(index_ak5_in_sameHemi_vetoca8_csvt));
                
                self.ttb_ht_[0] = ttb_ht;
                self.ttb_ca8_mass_pr_[0]       = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[theca8Index];
                self.ttb_ca8_ungroomed_pt_[0]  = getattr( self.InputTree_, "GroomedJet_CA8_pt" )[theca8Index];
                self.ttb_ca8_ungroomed_eta_[0] = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];
                self.ttb_ca8_ungroomed_phi_[0] = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];
                self.ttb_ca8_ungroomed_e_[0]   = getattr( self.InputTree_, "GroomedJet_CA8_e" )[theca8Index];
                
                if not self.IsData_ :
                 self.ttb_ca8_ungroomed_gen_pt_[0]  = getattr( self.InputTree_, "GenGroomedJet_CA8_pt" )[theca8Index];
                 self.ttb_ca8_ungroomed_gen_eta_[0] = getattr( self.InputTree_, "GenGroomedJet_CA8_eta" )[theca8Index];
                 self.ttb_ca8_ungroomed_gen_phi_[0] = getattr( self.InputTree_, "GenGroomedJet_CA8_phi" )[theca8Index];
                 self.ttb_ca8_ungroomed_gen_e_[0]   = getattr( self.InputTree_, "GenGroomedJet_CA8_e" )[theca8Index];

                self.ttb_ca8_charge_[0]     = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge" )[theca8Index];
#                self.ttb_ca8_charge_k05_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge_k05" )[theca8Index];
#                self.ttb_ca8_charge_k07_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge_k07" )[theca8Index];
#                self.ttb_ca8_charge_k10_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge_k10" )[theca8Index];
                 
                self.ttb_ca8_tau2tau1_[0]       = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1" )[theca8Index];
#                self.ttb_ca8_tau2tau1_exkT_[0]  = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1_exkT" )[theca8Index];
#                self.ttb_ca8_tau2tau1_pr_[0]    = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1_pr" )[theca8Index];
#                self.ttb_ca8_GeneralizedECF_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetGeneralizedECF" )[theca8Index];
                self.ttb_ca8_mu_[0]             = getattr( self.InputTree_, "GroomedJet_CA8_massdrop_pr" )[theca8Index];

                ttb_ca8J_p4  = ROOT.TLorentzVector();        
                ttb_ca8J_pt  = getattr( self.InputTree_, "GroomedJet_CA8_pt" )[theca8Index];    
                ttb_ca8J_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];    
                ttb_ca8J_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];    
                ttb_ca8J_e   = getattr( self.InputTree_, "GroomedJet_CA8_e" )[theca8Index];                        

                ttb_ca8J_p4.SetPtEtaPhiE(ttb_ca8J_pt, ttb_ca8J_eta, ttb_ca8J_phi, ttb_ca8J_e)

                ttb_V_p4 = ROOT.TLorentzVector(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type0"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type0_[0] = (ttb_V_p4+ttb_ca8J_p4).M();
                ttb_V_p4.SetPtEtaPhiE(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type2"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type2_[0] = (ttb_V_p4+ttb_ca8J_p4).M();

                ttb_V_p4.SetPtEtaPhiE(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type0_met"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type0_met_[0] = (ttb_V_p4+ttb_ca8J_p4).M();
                ttb_V_p4.SetPtEtaPhiE(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type2_met"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type2_met_[0] = (ttb_V_p4+ttb_ca8J_p4).M();

                self.ttb_ca8_px_[0] = ttb_ca8J_p4.Px();
                self.ttb_ca8_py_[0] = ttb_ca8J_p4.Py();
                self.ttb_ca8_pz_[0] = ttb_ca8J_p4.Pz();
                self.ttb_ca8_e_[0] =  ttb_ca8J_p4.E();
                ## preselection for the ttbar control region selection: at least one btag loose same or opposite
                oppo1same1 = (self.ttb_nak5_same_csvl_[0] > 0  or self.ttb_nak5_oppo_csvl_[0] > 0) or (self.ttb_nak5_sameveto_csvl_[0]>0 or self.ttb_nak5_oppoveto_csvl_[0]>0);
                oppo2same0 = (self.ttb_nak5_same_csvl_[0] == 0 and self.ttb_nak5_oppo_csvl_[0] > 1) or (self.ttb_nak5_sameveto_csvl_[0]==0 or self.ttb_nak5_oppoveto_csvl_[0] > 1);

                if oppo1same1 or oppo2same0:
                    self.isttbar_[0] = 1;
            else:

                self.ttb_nak5_same_[0] = -1;
                self.ttb_nak5_same_csvl_[0] = -1;
                self.ttb_nak5_same_csvm_[0] = -1;
                self.ttb_nak5_same_csvt_[0] = -1;

                self.ttb_nak5_oppo_[0] = -1;
                self.ttb_nak5_oppo_csvl_[0] = -1;
                self.ttb_nak5_oppo_csvm_[0] = -1;
                self.ttb_nak5_oppo_csvt_[0] = -1;
                self.ttb_nak5_oppoveto_[0] = -1;

                self.ttb_nak5_oppoveto_csvl_[0] = -1;
                self.ttb_nak5_oppoveto_csvm_[0] = -1;
                self.ttb_nak5_oppoveto_csvt_[0] = -1;

                self.ttb_nak5_sameveto_csvl_[0] = -1;
                self.ttb_nak5_sameveto_csvm_[0] = -1;
                self.ttb_nak5_sameveto_csvt_[0] = -1;

                self.ttb_ca8_mass_pr_[0]       = -1;
                self.ttb_ca8_ungroomed_pt_[0]  = -1;
                self.ttb_ca8_ungroomed_eta_[0] = -999;
                self.ttb_ca8_ungroomed_phi_[0] = -999;
                self.ttb_ca8_ungroomed_e_[0]   = -1;
                
                self.ttb_ca8_ungroomed_gen_pt_[0]  = -1; 
                self.ttb_ca8_ungroomed_gen_eta_[0] = -999; 
                self.ttb_ca8_ungroomed_gen_phi_[0] = -999; 
                self.ttb_ca8_ungroomed_gen_e_[0]   = -1;

                self.ttb_ca8_charge_[0]     = -999; 
                self.ttb_ca8_charge_k05_[0] = -999; 
                self.ttb_ca8_charge_k07_[0] = -999; 
                self.ttb_ca8_charge_k10_[0] = -999; 
                 
                self.ttb_ca8_tau2tau1_[0]       = -999; 
                self.ttb_ca8_tau2tau1_exkT_[0]  = -999; 
                self.ttb_ca8_tau2tau1_pr_[0]    = -999;
                self.ttb_ca8_GeneralizedECF_[0] = -999;
                self.ttb_ca8_mu_[0]             = -999;

                self.ttb_ca8_mlvj_type0_[0] = -999;
                self.ttb_ca8_mlvj_type2_[0] = -999;

                self.ttb_ca8_mlvj_type0_met_[0] = -999;
                self.ttb_ca8_mlvj_type2_met_[0] = -999;

                self.ttb_ca8_px_[0] = -999; 
                self.ttb_ca8_py_[0] = -999;
                self.ttb_ca8_pz_[0] = -999;
                self.ttb_ca8_e_[0] =  -999;
                               
            ##################### Cuts And divinding event in SR and TTbar Control Region
            
            leptonCut = 30; ## basic lepton thresolds
            leptonCutString = "W_muon_pt";
            metCut = 40;
            if self.Channel_ == "el": 
                leptonCut = 35;    
                leptonCutString = "W_electron_pt";
                metCut = 40;

            signallike = 0;
            
            if ( getattr( self.InputTree_, "W_pt" ) > 200 and
                 getattr( self.InputTree_, "GroomedJet_CA8_pt" )[0] > 200 and
                 math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_eta" )[0]) < 2.4 and
                 getattr( self.InputTree_, "event_met_pfmet" ) > metCut   and
                 getattr( self.InputTree_, leptonCutString ) > leptonCut  and
                 getattr( self.InputTree_, "GroomedJet_CA8_deltaR_lca8jet") > 1.57) :

                if ( self.Channel_ == "mu" and
                     math.fabs(getattr( self.InputTree_, "W_muon_dz000")) < 0.02 and
                     math.fabs(getattr( self.InputTree_, "W_muon_dzPV")) < 0.5   and
                     math.fabs( getattr( self.InputTree_, "W_muon_eta" )) < 2.1 ) :
                    signallike = 1 ; 
                elif self.Channel_ == "el" : signallike = 1 ; 

            ttbarlike = 0;
            if self.isttbar_[0] == 1 and getattr( self.InputTree_, "event_met_pfmet" ) > metCut and getattr( self.InputTree_, leptonCutString ) > leptonCut:
                ttbarlike = 1;

            ### store the info just for ttbar like or signal like events
            if ttbarlike == 1 or signallike == 1 :
              
             ############### Fill Event property only for interesting events
                             
             self.event_[0]            = getattr( self.InputTree_, "event_evtNo");
             self.event_runNo_[0]      = getattr( self.InputTree_, "event_runNo" );
             self.event_lumi_[0]       = getattr( self.InputTree_, "event_lumi" );
             
             effwt = getattr( self.InputTree_, "effwt" );
             puwt  = getattr( self.InputTree_, "puwt" );
             
             totSampleWeight = 1.;
             if self.IsData_: totSampleWeight = wSampleWeight;
             else: totSampleWeight = wSampleWeight*effwt*puwt; ## total weight takes into account also pileUp and efficiency (not the btag)

             self.totalEventWeight_[0]  = totSampleWeight;
             self.eff_and_pu_Weight_[0] = effwt*puwt;
             self.wSampleWeight_[0]     = wSampleWeight;

             self.btag_weight_[0]    = getattr( self.InputTree_, "eff_btag" );
             self.btag_weight_up_[0] = getattr( self.InputTree_, "eff_btag_up" );
             self.btag_weight_dn_[0] = getattr( self.InputTree_, "eff_btag_dw" );
             self.btag_weight_up_dn_[0] = getattr( self.InputTree_, "eff_btag_up_dw" );
             self.btag_weight_dn_up_[0] = getattr( self.InputTree_, "eff_btag_dw_up" );

             self.nPV_[0]      = getattr( self.InputTree_, "event_nPV" );

             self.issignal_[0] = signallike;
             self.isttbar_[0]  = ttbarlike;

             self.numberJetBin_[0] = getattr(self.InputTree_, "numberJetBin")[0];
             self.numberJetBin2_[0] = getattr(self.InputTree_, "numberJetBin")[1];
             self.numberJetBin3_[0] = getattr(self.InputTree_, "numberJetBin")[2];
             self.numberJetBin4_[0] = getattr(self.InputTree_, "numberJetBin")[3];

             if not self.IsData_:
              self.numberJetBinGen_[0] = getattr(self.InputTree_, "numberJetBinGen")[0];
              self.numberJetBinGen2_[0] = getattr(self.InputTree_, "numberJetBinGen")[1];
              self.numberJetBinGen3_[0] = getattr(self.InputTree_, "numberJetBinGen")[2];
              self.numberJetBinGen4_[0] = getattr(self.InputTree_, "numberJetBinGen")[3];

             if self.Label_ == "TTbar_mcatnlo" :
              self.event_weight_[0] = getattr( self.InputTree_, "event_weight" )/math.fabs(getattr( self.InputTree_, "event_weight" )) ;
             else: self.event_weight_[0] = 1.;


             #### CPS part -> take value from histos
             rwCPS = 1;
             if self.SignalMass_ > 0:
                 binVal = self.x_rwCPS.FindBin(getattr(self.InputTree_,"W_H_mass_gen"));
                 if binVal > self.h_rwCPS.GetNbinsX(): binVal = self.h_rwCPS.GetNbinsX();
                 if binVal < 1: binVal = 1;
                 rwCPS = self.h_rwCPS.GetBinContent( binVal );

             #interference weight
             self.complexpolewtggH600    = getattr(self.InputTree_,"complexpolewtggH600")*rwCPS;
             self.interferencewtggH600   = getattr(self.InputTree_,"interferencewtggH600");
             self.avecomplexpolewtggH600 = getattr(self.InputTree_,"avecomplexpolewtggH600"); 
             self.interference_Weight_H600_[0] = self.complexpolewtggH600*self.interferencewtggH600/self.avecomplexpolewtggH600;  ## complete weight for standard higgs

             self.complexpolewtggH700    = getattr(self.InputTree_,"complexpolewtggH700")*rwCPS; 
             self.interferencewtggH700   = getattr(self.InputTree_,"interferencewtggH700");
             self.avecomplexpolewtggH700 = getattr(self.InputTree_,"avecomplexpolewtggH700"); 
             self.interference_Weight_H700_[0] = self.complexpolewtggH700*self.interferencewtggH700/self.avecomplexpolewtggH700;  ## complete weight for standard higgs

             self.complexpolewtggH800    = getattr(self.InputTree_,"complexpolewtggH800")*rwCPS; 
             self.interferencewtggH800   = getattr(self.InputTree_,"interferencewtggH800");
             self.avecomplexpolewtggH800 = getattr(self.InputTree_,"avecomplexpolewtggH800"); 
             self.interference_Weight_H800_[0] = self.complexpolewtggH800*self.interferencewtggH800/self.avecomplexpolewtggH800;  ## complete weight for standard higgs

             self.complexpolewtggH900    = getattr(self.InputTree_,"complexpolewtggH900")*rwCPS; 
             self.interferencewtggH900   = getattr(self.InputTree_,"interferencewtggH900");
             self.avecomplexpolewtggH900 = getattr(self.InputTree_,"avecomplexpolewtggH900"); 
             self.interference_Weight_H900_[0] = self.complexpolewtggH900*self.interferencewtggH900/self.avecomplexpolewtggH900;  ## complete weight for standard higgs

             self.complexpolewtggH1000    = getattr(self.InputTree_,"complexpolewtggH1000")*rwCPS; 
             self.interferencewtggH1000   = getattr(self.InputTree_,"interferencewtggH1000");
             self.avecomplexpolewtggH1000 = getattr(self.InputTree_,"avecomplexpolewtggH1000"); 
             self.interference_Weight_H1000_[0] = self.complexpolewtggH1000*self.interferencewtggH1000/self.avecomplexpolewtggH1000;  ## complete weight for standard higgs

             self.cps_Weight_H600_[0] = self.complexpolewtggH600/self.avecomplexpolewtggH600;
             self.cps_Weight_H700_[0] = self.complexpolewtggH700/self.avecomplexpolewtggH700;
             self.cps_Weight_H800_[0] = self.complexpolewtggH800/self.avecomplexpolewtggH800;
             self.cps_Weight_H900_[0] = self.complexpolewtggH900/self.avecomplexpolewtggH900;
             self.cps_Weight_H1000_[0] = self.complexpolewtggH1000/self.avecomplexpolewtggH1000;

             #### produce weights for alternative models
        
             if self.SignalMass_ > 0:

                curIntfRw = getattr(self.InputTree_,"interferencewtggH%03d"%(self.SignalMass_));  ## take the interference value
                
                self.genHMass_[0] = getattr(self.InputTree_,"W_H_mass_gen");
                self.genHphi_[0]  = getattr(self.InputTree_,"W_H_phi_gen");
                self.genHeta_[0]  = getattr(self.InputTree_,"W_H_eta_gen");
                self.genHpt_[0]   = getattr(self.InputTree_,"W_H_pt_gen");
               
                if self.isVBF_:
                 self.genTagQuark1E_[0]    = getattr(self.InputTree_,"W_TagQuark_E")[0];
                 self.genTagQuark1eta_[0]  = getattr(self.InputTree_,"W_TagQuark_eta")[0];
                 self.genTagQuark1phi_[0]  = getattr(self.InputTree_,"W_TagQuark_phi")[0];
                 self.genTagQuark1pt_[0]   = getattr(self.InputTree_,"W_TagQuark_pt")[0];
                 self.genTagQuark2E_[0]    = getattr(self.InputTree_,"W_TagQuark_E")[1];
                 self.genTagQuark2eta_[0]  = getattr(self.InputTree_,"W_TagQuark_eta")[1];
                 self.genTagQuark2phi_[0]  = getattr(self.InputTree_,"W_TagQuark_phi")[1];
                 self.genTagQuark2pt_[0]   = getattr(self.InputTree_,"W_TagQuark_pt")[1];
                                                             

                for i in range(len(self.cprimeVals)): ## run over the possible value of the new couplig constant and BR
                 for j in range(len(self.brnewVals)): 
                  curCprime = float(self.cprimeVals[i])/10.;
                  curBRnew = float(self.brnewVals[j])/10.;
                  if self.isVBF_: 
                    self.bsmReweights[i][j][0] = self.GetInteferenceWeights( getattr(self.InputTree_,"W_H_mass_gen"), curCprime, curBRnew );  
                  else:
                   self.bsmReweights[i][j][0] = self.GetInteferenceWeights( getattr(self.InputTree_,"W_H_mass_gen"), curCprime, curBRnew )*IntfRescale(curIntfRw,curCprime,curBRnew);
                
             else:   
                    self.genHMass_[0] = -1;
                    self.genHphi_[0]  = -999;
                    self.genHeta_[0]  = -999;
                    self.genHpt_[0]   = -1;
                    self.genTagQuark1E_ = -1;
                    self.genTagQuark1eta_  = -999.;
                    self.genTagQuark1phi_  = -999.;
                    self.genTagQuark1pt_   = -1.;
                    self.genTagQuark2E_ = -1;
                    self.genTagQuark2eta_  = -999.;
                    self.genTagQuark2phi_  = -999.;
                    self.genTagQuark2pt_   = -1.;
                                                             
                    for i in range(len(self.cprimeVals)): 
                        for j in range(len(self.brnewVals)): 
                            self.bsmReweights[i][j][0] = -1;
            
                
             ################# Other basic Observables
                            
             self.mass_lvj_type0_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type0" );
             self.mass_lvj_type2_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type2" );

             self.mass_lvj_type0_met_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type0_met" );
             self.mass_lvj_type2_met_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type2_met" );

             self.mass_lv_subj_type0_[0] = getattr( self.InputTree_, "boosted_lvj_m_type0" );
             self.mass_lv_subj_type2_[0] = getattr( self.InputTree_, "boosted_lvj_m_type2" );

             self.mass_lv_subj_type0_met_[0] = getattr( self.InputTree_, "boosted_lvj_m_type0_met" );
             self.mass_lv_subj_type2_met_[0] = getattr( self.InputTree_, "boosted_lvj_m_type2_met" );

             self.v_pt_[0]   = getattr( self.InputTree_, "W_pt" );
             self.v_mt_[0]   = getattr( self.InputTree_, "W_mt" );
             self.v_eta_[0]  = getattr( self.InputTree_, "W_eta" );
             self.v_phi_[0]  = getattr( self.InputTree_, "W_phi" );

             self.nu_pz_type0_[0] = getattr( self.InputTree_, "W_nu1_pz_type0" );
             self.nu_pz_type2_[0] = getattr( self.InputTree_, "W_nu1_pz_type2" );

             self.nu_pz_type0_met_[0] = getattr( self.InputTree_, "W_nu1_pz_type0_met" );
             self.nu_pz_type2_met_[0] = getattr( self.InputTree_, "W_nu1_pz_type2_met" );

             self.W_pz_type0_[0] = getattr( self.InputTree_, "W_pz_type0" );
             self.W_pz_type2_[0] = getattr( self.InputTree_, "W_pz_type2" );

             self.W_pz_type0_met_[0] = getattr( self.InputTree_, "W_pz_type0_met" );
             self.W_pz_type2_met_[0] = getattr( self.InputTree_, "W_pz_type2_met" );

             W_Lepton_gen = ROOT.TLorentzVector() ;
             
             if self.IsData_ == False  :
                 self.nu_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" );                 
                 if self.Channel_ == "mu" :
                     W_Lepton_gen.SetPxPyPzE(getattr(self.InputTree_, "W_muon_px_gen" )+getattr(self.InputTree_, "W_neutrino_px_gen"),
                                             getattr(self.InputTree_, "W_muon_py_gen" )+getattr(self.InputTree_, "W_neutrino_py_gen" ),
                                             getattr(self.InputTree_, "W_muon_pz_gen" )+getattr(self.InputTree_, "W_neutrino_pz_gen" ),
                                             getattr(self.InputTree_, "W_muon_e_gen" )+getattr(self.InputTree_, "W_neutrino_e_gen" ));                     

                     self.W_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" ) + getattr( self.InputTree_, "W_muon_pz_gen" );
                     self.W_pt_gen_[0]  =  W_Lepton_gen.Pt();
                 
                 elif self.Channel_ == "el" :
                     W_Lepton_gen.SetPxPyPzE(getattr(self.InputTree_, "W_electron_px_gen" )+getattr(self.InputTree_, "W_neutrino_px_gen"),
                                             getattr(self.InputTree_, "W_electron_py_gen" )+getattr(self.InputTree_, "W_neutrino_py_gen" ),
                                             getattr(self.InputTree_, "W_electron_pz_gen" )+getattr(self.InputTree_, "W_neutrino_pz_gen" ),
                                             getattr(self.InputTree_, "W_electron_e_gen" )+getattr(self.InputTree_, "W_neutrino_e_gen" ));

                     self.W_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" ) + getattr( self.InputTree_, "W_electron_pz_gen" );
                     self.W_pt_gen_[0]  =  W_Lepton_gen.Pt();

                 
                 
             ######## jet stuff
             self.ungroomed_jet_pt_[0]  = getattr( self.InputTree_, prefix+"_pt" )[0];
             self.ungroomed_jet_eta_[0] = getattr( self.InputTree_, prefix+"_eta" )[0];
             self.ungroomed_jet_phi_[0] = getattr( self.InputTree_, prefix+"_phi" )[0];
             self.ungroomed_jet_e_[0]   = getattr( self.InputTree_, prefix+"_e" )[0];

             self.jet_mass_pr_[0] = getattr( self.InputTree_, prefix + "_mass_pr" )[0];
             self.jet_pt_pr_[0]   = getattr( self.InputTree_, prefix + "_pt_pr" )[0];

             self.jet_charge_[0]     = getattr( self.InputTree_, prefix + "_jetcharge" )[0];
#             self.jet_charge_k05_[0] = getattr( self.InputTree_, prefix + "_jetcharge_k05" )[0];
#             self.jet_charge_k07_[0] = getattr( self.InputTree_, prefix + "_jetcharge_k07" )[0];
#             self.jet_charge_k10_[0] = getattr( self.InputTree_, prefix + "_jetcharge_k10" )[0];

             self.jet_grsens_ft_[0]   = getattr( self.InputTree_, prefix + "_mass_ft" )[0] / getattr( self.InputTree_, prefix + "_mass" )[0];
             self.jet_grsens_tr_[0]   = getattr( self.InputTree_, prefix + "_mass_tr" )[0] / getattr( self.InputTree_, prefix + "_mass" )[0];
             self.jet_massdrop_pr_[0] = getattr( self.InputTree_, prefix + "_massdrop_pr" )[0];    

             ### CA8 jet collection
             self.jecUnc_.setJetEta( getattr( self.InputTree_, prefix+"_eta" )[0]);
             self.jecUnc_.setJetPt( getattr( self.InputTree_, prefix+"_pt" )[0]);                        
             self.j_jecfactor_up_[0] = self.jecUnc_.getUncertainty( True );

             self.jecUnc_.setJetEta( getattr( self.InputTree_, prefix+"_eta" )[0]);
             self.jecUnc_.setJetPt( getattr( self.InputTree_, prefix+"_pt" )[0]);               
             self.j_jecfactor_dn_[0] = self.jecUnc_.getUncertainty( False ) ;

             self.j_jecfactor_up_[0] = math.sqrt( self.j_jecfactor_up_[0]**2 + 0.02**2 );
             self.j_jecfactor_dn_[0] = math.sqrt( self.j_jecfactor_dn_[0]**2 + 0.02**2 );
             self.curjes_up = 1 + self.j_jecfactor_up_[0]
             self.curjes_dn = 1 - self.j_jecfactor_dn_[0]

             jorig_pt  = getattr( self.InputTree_, prefix + "_pt_pr" )[0];    
             jorig_eta = getattr( self.InputTree_, prefix + "_eta_pr" )[0];    
             jorig_phi = getattr( self.InputTree_, prefix + "_phi_pr" )[0];    
             jorig_e   = getattr( self.InputTree_, prefix + "_e_pr" )[0];                        
             jdef_ptetaphie = ROOT.TLorentzVector();
             jdef_ptetaphie.SetPtEtaPhiE(jorig_pt, jorig_eta, jorig_phi, jorig_e)

             jdef_up = ROOT.TLorentzVector(jdef_ptetaphie.Px() * self.curjes_up, jdef_ptetaphie.Py() * self.curjes_up,
                                           jdef_ptetaphie.Pz() * self.curjes_up, jdef_ptetaphie.E() * self.curjes_up);
             jdef_dn = ROOT.TLorentzVector(jdef_ptetaphie.Px() * self.curjes_dn, jdef_ptetaphie.Py() * self.curjes_dn,
                                           jdef_ptetaphie.Pz() * self.curjes_dn, jdef_ptetaphie.E() * self.curjes_dn);

             self.jet_mass_pr_up_[0] = jdef_up.M();  
             self.jet_mass_pr_dn_[0] = jdef_dn.M();
              
             self.jet_ungroomed_jet_pt_up_[0] = jdef_up.Pt();
             self.jet_ungroomed_jet_pt_dn_[0] = jdef_dn.Pt();

             qjetmassdistribution = getattr( self.InputTree_, prefix+"_qjetmass" );
             qjetvol              = getListRMS(qjetmassdistribution)/getListMean(qjetmassdistribution);
             self.jet_qjetvol_[0] = qjetvol;

             self.jet_tau2tau1_[0]        = getattr( self.InputTree_, prefix + "_tau2tau1" )[0];     
#             self.jet_tau2tau1_exkT_[0]   = getattr( self.InputTree_, prefix + "_tau2tau1_exkT" )[0];     
#             self.jet_tau2tau1_pr_[0]     = getattr( self.InputTree_, prefix + "_tau2tau1_pr" )[0];     
#             self.jet_GeneralizedECF_[0]  = getattr( self.InputTree_, prefix + "_jetGeneralizedECF" )[0];     
             self.jet_jetconstituents_[0] = getattr( self.InputTree_, prefix + "_jetconstituents" )[0];     

             self.jet_rcore4_[0] = getattr( self.InputTree_, prefix + "_rcores")[3*6 + 0];
             self.jet_rcore5_[0] = getattr( self.InputTree_, prefix + "_rcores")[4*6 + 0];
             self.jet_rcore6_[0] = getattr( self.InputTree_, prefix + "_rcores")[5*6 + 0];
             self.jet_rcore7_[0] = getattr( self.InputTree_, prefix + "_rcores")[6*6 + 0];

             self.jet_planarlow04_[0] = getattr( self.InputTree_, prefix + "_planarflow04");
             self.jet_planarlow05_[0] = getattr( self.InputTree_, prefix + "_planarflow05");
             self.jet_planarlow06_[0] = getattr( self.InputTree_, prefix + "_planarflow06");
             self.jet_planarlow07_[0] = getattr( self.InputTree_, prefix + "_planarflow07");

             self.pt1FracVal = max( getattr( self.InputTree_, prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, prefix + "_prsubjet2ptoverjetpt" ) );
             self.pt2FracVal = min( getattr( self.InputTree_, prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, prefix + "_prsubjet2ptoverjetpt" ) );

             self.jet_pt1frac_[0] = self.pt1FracVal;
             self.jet_pt2frac_[0] = self.pt2FracVal;

             self.jet_sjdr_[0] = getattr( self.InputTree_, prefix + "_prsubjet1subjet2_deltaR" );       
                
             self.deltaR_lca8jet_[0] = getattr( self.InputTree_, prefix + "_deltaR_lca8jet" );       

             self.deltaphi_METca8jet_[0]     = getattr( self.InputTree_, prefix + "_deltaphi_METca8jet_type2" );       
             self.deltaphi_Vca8jet_[0]       = getattr( self.InputTree_, prefix + "_deltaphi_Vca8jet_type2" );       
             self.deltaphi_METca8jet_met_[0] = getattr( self.InputTree_, prefix + "_deltaphi_METca8jet_type2_met" );       
             self.deltaphi_Vca8jet_met_[0]   = getattr( self.InputTree_, prefix + "_deltaphi_Vca8jet_type2_met" );       

             if not self.IsData_ :

              self.ungroomed_gen_jet_pt_[0]  = getattr( self.InputTree_, "Gen"+prefix+"_pt" )[0];
              self.ungroomed_gen_jet_eta_[0] = getattr( self.InputTree_, "Gen"+prefix+"_eta" )[0];
              self.ungroomed_gen_jet_phi_[0] = getattr( self.InputTree_, "Gen"+prefix+"_phi" )[0];
              self.ungroomed_gen_jet_e_[0]   = getattr( self.InputTree_, "Gen"+prefix+"_e" )[0];
 
              self.gen_jet_mass_pr_       = getattr( self.InputTree_, "Gen"+prefix + "_mass_pr" )[0];
              self.gen_jet_pt_pr_         = getattr( self.InputTree_, "Gen"+prefix + "_pt_pr" )[0];
#              self.gen_jet_charge_        = getattr( self.InputTree_, "Gen"+prefix + "_charge" )[0];
#              self.gen_jet_charge_k05_    = getattr( self.InputTree_, "Gen"+prefix + "_charge_k05" )[0];
#              self.gen_jet_charge_k07_    = getattr( self.InputTree_, "Gen"+prefix + "_charge_k07" )[0];
#              self.gen_jet_charge_k10_    = getattr( self.InputTree_, "Gen"+prefix + "_charge_k10" )[0];

              self.gen_jet_grsens_ft_[0]   = getattr( self.InputTree_, "Gen"+prefix + "_mass_ft" )[0] / getattr( self.InputTree_, "Gen"+prefix + "_mass" )[0];
              self.gen_jet_grsens_tr_[0]   = getattr( self.InputTree_, "Gen"+prefix + "_mass_tr" )[0] / getattr( self.InputTree_, "Gen"+prefix + "_mass" )[0];
              self.gen_jet_massdrop_pr_[0] = getattr( self.InputTree_, "Gen"+prefix + "_massdrop_pr" )[0];    

#              qjetmassdistribution     = getattr( self.InputTree_, "Gen"+prefix+"_qjetmass" );
#              qjetvol                  = getListRMS(qjetmassdistribution)/getListMean(qjetmassdistribution);
#              self.gen_jet_qjetvol_[0] = qjetvol;

              self.gen_jet_tau2tau1_[0]        = getattr( self.InputTree_, "Gen"+prefix + "_tau2tau1" )[0];     
#              self.gen_jet_tau2tau1_exkT_[0]   = getattr( self.InputTree_, "Gen"+prefix + "_tau2tau1_exkT" )[0];     
#              self.gen_jet_tau2tau1_pr_[0]     = getattr( self.InputTree_, "Gen"+prefix + "_tau2tau1_pr" )[0];     
#              self.gen_jet_jetconstituents_[0] = getattr( self.InputTree_, "Gen"+prefix + "_jetconstituents" )[0];     

              self.gen_jet_rcore4_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[3*6 + 0];
              self.gen_jet_rcore5_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[4*6 + 0];
              self.gen_jet_rcore6_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[5*6 + 0];
              self.gen_jet_rcore7_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[6*6 + 0];
 
#              self.gen_jet_planarlow04_[0] = getattr( self.InputTree_, "Gen"+prefix + "_planarflow04");
#              self.gen_jet_planarlow05_[0] = getattr( self.InputTree_, "Gen"+prefix + "_planarflow05");
#              self.gen_jet_planarlow06_[0] = getattr( self.InputTree_, "Gen"+prefix + "_planarflow06");
#              self.gen_jet_planarlow07_[0] = getattr( self.InputTree_, "Gen"+prefix + "_planarflow07");

#              self.pt1FracVal = max( getattr( self.InputTree_, "Gen"+prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, "Gen"+prefix + "_prsubjet2ptoverjetpt" ) );
#              self.pt2FracVal = min( getattr( self.InputTree_, "Gen"+prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, "Gen"+prefix + "_prsubjet2ptoverjetpt" ) );

#              self.gen_jet_pt1frac_[0] = pt1FracVal;
#              self.gen_jet_pt2frac_[0] = pt2FracVal;

#              self.gen_jet_sjdr_[0] = getattr( self.InputTree_,"Gen"+prefix + "_prsubjet1subjet2_deltaR" );       

             ### lepton and met side                
             if self.Channel_ == "mu" :
                    self.l_pt_[0]     = getattr( self.InputTree_, "W_muon_pt" );
                    self.l_eta_[0]    = getattr( self.InputTree_, "W_muon_eta" );
                    self.l_phi_[0]    = getattr( self.InputTree_, "W_muon_phi" );
                    self.l_charge_[0] = getattr( self.InputTree_, "W_muon_charge" );
                    
             elif self.Channel_ == "el":
                    self.l_pt_[0]     = getattr( self.InputTree_, "W_electron_pt" );
                    self.l_eta_[0]    = getattr( self.InputTree_, "W_electron_eta" );
                    self.l_phi_[0]    = getattr( self.InputTree_, "W_electron_phi" );
                    self.l_charge_[0] = getattr( self.InputTree_, "W_electron_charge" );

             self.pfMET_[0]     = getattr( self.InputTree_, "event_met_pfmet" );        
             self.pfMET_Phi_[0] = getattr( self.InputTree_, "event_met_pfmetPhi" );        

             ######## VBF jet Stuff
             if self.numberJetBin_[0] ==1 :

 
              self.vbf_maxpt_j1_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m");                               
              self.vbf_maxpt_j1_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt");                               
              self.vbf_maxpt_j1_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta");                               
              self.vbf_maxpt_j1_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi");                               
 
              self.vbf_maxpt_j1_QGLikelihood_[0] = getattr(self.InputTree_, "vbf_maxpt_j1_QGLikelihood");                                
 
              self.vbf_maxpt_j1_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpMedium");                               
              self.vbf_maxpt_j1_isPileUpTight_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpTight");                               

              self.vbf_maxpt_j1_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV");

              if self.numberJetBinGen_[0] ==1 and not self.IsData_:

               self.vbf_maxpt_j1_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m_gen");                               
               self.vbf_maxpt_j1_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt_gen");                               
               self.vbf_maxpt_j1_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta_gen");                               
               self.vbf_maxpt_j1_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi_gen");                               

               self.vbf_maxpt_j1_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV_gen");

                  
             if self.numberJetBin_[0] >=2 :

              self.vbf_maxpt_jj_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_jj_m");                               
              self.vbf_maxpt_jj_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_jj_pt");                               
              self.vbf_maxpt_jj_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_eta");                               
              self.vbf_maxpt_jj_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_phi");                               

              self.vbf_maxpt_j1_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m");                               
              self.vbf_maxpt_j1_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt");                               
              self.vbf_maxpt_j1_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta");                               
              self.vbf_maxpt_j1_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi");                               
 
              self.vbf_maxpt_j2_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_j2_m");                               
              self.vbf_maxpt_j2_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_j2_pt");                               
              self.vbf_maxpt_j2_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_eta");                               
              self.vbf_maxpt_j2_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_phi");                               

              self.vbf_maxpt_j1_QGLikelihood_[0] = getattr(self.InputTree_, "vbf_maxpt_j1_QGLikelihood");                                
              self.vbf_maxpt_j2_QGLikelihood_[0] = getattr(self.InputTree_, "vbf_maxpt_j2_QGLikelihood");                               
 
              self.vbf_maxpt_j1_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpMedium");                               
              self.vbf_maxpt_j1_isPileUpTight_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpTight");                               

              self.vbf_maxpt_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_isPileUpMedium");                               
              self.vbf_maxpt_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_isPileUpMedium");                               

              self.vbf_maxpt_j1_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV");
              self.vbf_maxpt_j2_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_bDiscriminatorCSV");                                

              if not self.IsData_ and self.numberJetBinGen_[0]>= 2:
                  
               self.vbf_maxpt_jj_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_jj_m_gen");                               
               self.vbf_maxpt_jj_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_jj_pt_gen");                               
               self.vbf_maxpt_jj_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_eta_gen");                               
               self.vbf_maxpt_jj_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_phi_gen");                               

               self.vbf_maxpt_j1_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m_gen");                               
               self.vbf_maxpt_j1_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt_gen");                               
               self.vbf_maxpt_j1_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta_gen");                               
               self.vbf_maxpt_j1_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi_gen");                               
 
               self.vbf_maxpt_j2_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_j2_m_gen");                               
               self.vbf_maxpt_j2_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_j2_pt_gen");                               
               self.vbf_maxpt_j2_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_eta_gen");                               
               self.vbf_maxpt_j2_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_phi_gen");                               
 
               self.vbf_maxpt_j1_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV_gen");
               self.vbf_maxpt_j2_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_bDiscriminatorCSV_gen");                                

              ### scaling jets
              j_jecfactorAK5_up_ = array('f',[0.]);
              j_jecfactorAK5_dn_ = array('f',[0.]);
              self.jecUncAK5_.setJetEta(getattr(self.InputTree_, "vbf_maxpt_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr(self.InputTree_, "vbf_maxpt_j1_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta(getattr( self.InputTree_, "vbf_maxpt_j1_eta" ));
              self.jecUncAK5_.setJetPt(getattr( self.InputTree_, "vbf_maxpt_j1_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              self.vbf_maxpt_j1 = ROOT.TLorentzVector();
              self.vbf_maxpt_j1.SetPtEtaPhiM(self.vbf_maxpt_j1_pt_[0],self.vbf_maxpt_j1_eta_[0],self.vbf_maxpt_j1_phi_[0],self.vbf_maxpt_j1_m_[0]); 
              vbf_maxpt_j1_up = ROOT.TLorentzVector(self.vbf_maxpt_j1.Px()*(1+j_jecfactorAK5_up_[0]),
                                                    self.vbf_maxpt_j1.Py()*(1+j_jecfactorAK5_up_[0]),
                                                    self.vbf_maxpt_j1.Pz()*(1+j_jecfactorAK5_up_[0]),
                                                    self.vbf_maxpt_j1.E()*(1+j_jecfactorAK5_up_[0]));
                                                   

              self.vbf_maxpt_j1_m_up_[0]   = vbf_maxpt_j1_up.M();                                
              self.vbf_maxpt_j1_pt_up_[0]  = vbf_maxpt_j1_up.Pt();                              
              self.vbf_maxpt_j1_eta_up_[0] = vbf_maxpt_j1_up.Eta();                                
              self.vbf_maxpt_j1_phi_up_[0] = vbf_maxpt_j1_up.Phi();                                

              vbf_maxpt_j1_dn = ROOT.TLorentzVector(self.vbf_maxpt_j1.Px()*(1+j_jecfactorAK5_dn_[0]),
                                                    self.vbf_maxpt_j1.Py()*(1+j_jecfactorAK5_dn_[0]),
                                                    self.vbf_maxpt_j1.Pz()*(1+j_jecfactorAK5_dn_[0]),
                                                    self.vbf_maxpt_j1.E()*(1+j_jecfactorAK5_dn_[0]));
                                                    

              self.vbf_maxpt_j1_m_dn_[0]   = vbf_maxpt_j1_dn.M();                                
              self.vbf_maxpt_j1_pt_dn_[0]  = vbf_maxpt_j1_dn.Pt();                              
              self.vbf_maxpt_j1_eta_dn_[0] = vbf_maxpt_j1_dn.Eta();                                
 
              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxpt_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxpt_j2_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );
  
              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxpt_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxpt_j2_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              self.vbf_maxpt_j2 = ROOT.TLorentzVector();
              self.vbf_maxpt_j2.SetPtEtaPhiM(self.vbf_maxpt_j2_pt_[0],self.vbf_maxpt_j2_eta_[0],self.vbf_maxpt_j2_phi_[0],self.vbf_maxpt_j2_m_[0]); 
              vbf_maxpt_j2_up = ROOT.TLorentzVector(self.vbf_maxpt_j2.Px()*(1+j_jecfactorAK5_up_[0]),
                                                    self.vbf_maxpt_j2.Py()*(1+j_jecfactorAK5_up_[0]),
                                                    self.vbf_maxpt_j2.Pz()*(1+j_jecfactorAK5_up_[0]),
                                                    self.vbf_maxpt_j2.E()*(1+j_jecfactorAK5_up_[0]));
                                                   

              self.vbf_maxpt_j2_m_up_[0]   = vbf_maxpt_j2_up.M();                                
              self.vbf_maxpt_j2_pt_up_[0]  = vbf_maxpt_j2_up.Pt();                              
              self.vbf_maxpt_j2_eta_up_[0] = vbf_maxpt_j2_up.Eta();                                
              self.vbf_maxpt_j2_phi_up_[0] = vbf_maxpt_j2_up.Phi();                                

              vbf_maxpt_j2_dn = ROOT.TLorentzVector(self.vbf_maxpt_j2.Px()*(1+j_jecfactorAK5_dn_[0]),
                                                    self.vbf_maxpt_j2.Py()*(1+j_jecfactorAK5_dn_[0]),
                                                    self.vbf_maxpt_j2.Pz()*(1+j_jecfactorAK5_dn_[0]),
                                                    self.vbf_maxpt_j2.E()*(1+j_jecfactorAK5_dn_[0]));

              self.vbf_maxpt_j2_m_dn_[0]   = vbf_maxpt_j2_dn.M();                                
              self.vbf_maxpt_j2_pt_dn_[0]  = vbf_maxpt_j2_dn.Pt();                              
              self.vbf_maxpt_j2_eta_dn_[0] = vbf_maxpt_j2_dn.Eta();                                


              self.vbf_maxDeta_jj_m_[0]     = getattr(self.InputTree_, "vbf_maxDeta_jj_m");                               
              self.vbf_maxDeta_jj_pt_[0]    = getattr(self.InputTree_, "vbf_maxDeta_jj_pt");                               
              self.vbf_maxDeta_jj_eta_[0]   = getattr(self.InputTree_, "vbf_maxDeta_jj_eta");                               
              self.vbf_maxDeta_jj_phi_[0]   = getattr(self.InputTree_, "vbf_maxDeta_jj_phi");                               

              self.vbf_maxDeta_j1_m_[0]     = getattr(self.InputTree_, "vbf_maxDeta_j1_m");                               
              self.vbf_maxDeta_j1_pt_[0]    = getattr(self.InputTree_, "vbf_maxDeta_j1_pt");                               
              self.vbf_maxDeta_j1_eta_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j1_eta");                               
              self.vbf_maxDeta_j1_phi_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j1_phi");                               

              self.vbf_maxDeta_j2_m_[0]     = getattr(self.InputTree_, "vbf_maxDeta_j2_m");                               
              self.vbf_maxDeta_j2_pt_[0]    = getattr(self.InputTree_, "vbf_maxDeta_j2_pt");                               
              self.vbf_maxDeta_j2_eta_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j2_eta");                               
              self.vbf_maxDeta_j2_phi_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j2_phi");                               

              self.vbf_maxDeta_j1_QGLikelihood_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j1_QGLikelihood");                                
              self.vbf_maxDeta_j2_QGLikelihood_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j2_QGLikelihood");                               

              self.vbf_maxDeta_j1_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j1_isPileUpMedium");                               
              self.vbf_maxDeta_j1_isPileUpTight_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j1_isPileUpTight");                               

              self.vbf_maxDeta_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j2_isPileUpMedium");                               
              self.vbf_maxDeta_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j2_isPileUpMedium");                               

              self.vbf_maxDeta_j1_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j1_bDiscriminatorCSV");
              if self.Channel_ == "mu":
               self.vbf_maxDeta_j2_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j2_bDiscriminatorCSV");                                

              if not self.IsData_ and self.numberJetBinGen_[0]>= 2:
                  
               self.vbf_maxDeta_jj_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxDeta_jj_m_gen");                               
               self.vbf_maxDeta_jj_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxDeta_jj_pt_gen");                               
               self.vbf_maxDeta_jj_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxDeta_jj_eta_gen");                               
               self.vbf_maxDeta_jj_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxDeta_jj_phi_gen");                               

               self.vbf_maxDeta_j1_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxDeta_j1_m_gen");                               
               self.vbf_maxDeta_j1_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxDeta_j1_pt_gen");                               
               self.vbf_maxDeta_j1_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j1_eta_gen");                               
               self.vbf_maxDeta_j1_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j1_phi_gen");                               
 
               self.vbf_maxDeta_j2_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxDeta_j2_m_gen");                               
               self.vbf_maxDeta_j2_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxDeta_j2_pt_gen");                               
               self.vbf_maxDeta_j2_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j2_eta_gen");                               
               self.vbf_maxDeta_j2_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxDeta_j2_phi_gen");                               
 
               self.vbf_maxDeta_j1_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j1_bDiscriminatorCSV_gen");
               self.vbf_maxDeta_j2_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxDeta_j2_bDiscriminatorCSV_gen");                                

              ### scaling jets
              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxDeta_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxDeta_j1_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxDeta_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxDeta_j1_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              self.vbf_maxDeta_j1 = ROOT.TLorentzVector();
              self.vbf_maxDeta_j1.SetPtEtaPhiM(self.vbf_maxDeta_j1_pt_[0],self.vbf_maxDeta_j1_eta_[0],self.vbf_maxDeta_j1_phi_[0],self.vbf_maxDeta_j1_m_[0]); 
              vbf_maxDeta_j1_up = ROOT.TLorentzVector(self.vbf_maxDeta_j1.Px()*(1+j_jecfactorAK5_up_[0]),
                                                      self.vbf_maxDeta_j1.Py()*(1+j_jecfactorAK5_up_[0]),
                                                      self.vbf_maxDeta_j1.Pz()*(1+j_jecfactorAK5_up_[0]),
                                                      self.vbf_maxDeta_j1.E()*(1+j_jecfactorAK5_up_[0]));
                                                   

              self.vbf_maxDeta_j1_m_up_[0]   = vbf_maxDeta_j1_up.M();                                
              self.vbf_maxDeta_j1_pt_up_[0]  = vbf_maxDeta_j1_up.Pt();                              
              self.vbf_maxDeta_j1_eta_up_[0] = vbf_maxDeta_j1_up.Eta();                                
              self.vbf_maxDeta_j1_phi_up_[0] = vbf_maxDeta_j1_up.Phi();                                

              vbf_maxDeta_j1_dn = ROOT.TLorentzVector(self.vbf_maxDeta_j1.Px()*(1+j_jecfactorAK5_dn_[0]),
                                                      self.vbf_maxDeta_j1.Py()*(1+j_jecfactorAK5_dn_[0]),
                                                      self.vbf_maxDeta_j1.Pz()*(1+j_jecfactorAK5_dn_[0]),
                                                      self.vbf_maxDeta_j1.E()*(1+j_jecfactorAK5_dn_[0]));
                                                    
              self.vbf_maxDeta_j1_m_dn_[0]   = vbf_maxDeta_j1_dn.M();                                
              self.vbf_maxDeta_j1_pt_dn_[0]  = vbf_maxDeta_j1_dn.Pt();                              
              self.vbf_maxDeta_j1_eta_dn_[0] = vbf_maxDeta_j1_dn.Eta();                                

              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxDeta_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxDeta_j2_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxDeta_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxDeta_j2_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              self.vbf_maxDeta_j2 = ROOT.TLorentzVector();
              self.vbf_maxDeta_j2.SetPtEtaPhiM(self.vbf_maxDeta_j2_pt_[0],self.vbf_maxDeta_j2_eta_[0],self.vbf_maxDeta_j2_phi_[0],self.vbf_maxDeta_j2_m_[0]); 
              vbf_maxDeta_j2_up = ROOT.TLorentzVector(self.vbf_maxDeta_j2.Px()*(1+j_jecfactorAK5_up_[0]),
                                                      self.vbf_maxDeta_j2.Py()*(1+j_jecfactorAK5_up_[0]),
                                                      self.vbf_maxDeta_j2.Pz()*(1+j_jecfactorAK5_up_[0]),
                                                      self.vbf_maxDeta_j2.E()*(1+j_jecfactorAK5_up_[0]));
                                                   

              self.vbf_maxDeta_j2_m_up_[0]   = vbf_maxDeta_j2_up.M();                                
              self.vbf_maxDeta_j2_pt_up_[0]  = vbf_maxDeta_j2_up.Pt();                              
              self.vbf_maxDeta_j2_eta_up_[0] = vbf_maxDeta_j2_up.Eta();                                
              self.vbf_maxDeta_j2_phi_up_[0] = vbf_maxDeta_j2_up.Phi();                                

              vbf_maxDeta_j2_dn = ROOT.TLorentzVector(self.vbf_maxDeta_j2.Px()*(1+j_jecfactorAK5_dn_[0]),
                                                      self.vbf_maxDeta_j2.Py()*(1+j_jecfactorAK5_dn_[0]),
                                                      self.vbf_maxDeta_j2.Pz()*(1+j_jecfactorAK5_dn_[0]),
                                                      self.vbf_maxDeta_j2.E()*(1+j_jecfactorAK5_dn_[0]));
                                                   

              self.vbf_maxDeta_j2_m_dn_[0]   = vbf_maxDeta_j2_dn.M();                                
              self.vbf_maxDeta_j2_pt_dn_[0]  = vbf_maxDeta_j2_dn.Pt();                              
              self.vbf_maxDeta_j2_eta_dn_[0] = vbf_maxDeta_j2_dn.Eta();                                


              self.vbf_maxMjj_jj_m_[0]     = getattr(self.InputTree_, "vbf_maxMjj_jj_m");                               
              self.vbf_maxMjj_jj_pt_[0]    = getattr(self.InputTree_, "vbf_maxMjj_jj_pt");                               
              self.vbf_maxMjj_jj_eta_[0]   = getattr(self.InputTree_, "vbf_maxMjj_jj_eta");                               
              self.vbf_maxMjj_jj_phi_[0]   = getattr(self.InputTree_, "vbf_maxMjj_jj_phi");                               

              self.vbf_maxMjj_j1_m_[0]     = getattr(self.InputTree_, "vbf_maxMjj_j1_m");                               
              self.vbf_maxMjj_j1_pt_[0]    = getattr(self.InputTree_, "vbf_maxMjj_j1_pt");                               
              self.vbf_maxMjj_j1_eta_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j1_eta");                               
              self.vbf_maxMjj_j1_phi_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j1_phi");                               

              self.vbf_maxMjj_j2_m_[0]     = getattr(self.InputTree_, "vbf_maxMjj_j2_m");                               
              self.vbf_maxMjj_j2_pt_[0]    = getattr(self.InputTree_, "vbf_maxMjj_j2_pt");                               
              self.vbf_maxMjj_j2_eta_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j2_eta");                               
              self.vbf_maxMjj_j2_phi_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j2_phi");                               
 
              self.vbf_maxMjj_j1_QGLikelihood_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j1_QGLikelihood");                                
              self.vbf_maxMjj_j2_QGLikelihood_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j2_QGLikelihood");                               

              self.vbf_maxMjj_j1_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j1_isPileUpMedium");                               
              self.vbf_maxMjj_j1_isPileUpTight_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j1_isPileUpTight");                               

              self.vbf_maxMjj_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j2_isPileUpMedium");                               
              self.vbf_maxMjj_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j2_isPileUpMedium");                               

              self.vbf_maxMjj_j1_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j1_bDiscriminatorCSV");                                
              self.vbf_maxMjj_j2_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j2_bDiscriminatorCSV");                                

              if not self.IsData_ and self.numberJetBinGen_[0]>= 2:
                  
               self.vbf_maxMjj_jj_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxMjj_jj_m_gen");                               
               self.vbf_maxMjj_jj_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxMjj_jj_pt_gen");                               
               self.vbf_maxMjj_jj_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxMjj_jj_eta_gen");                               
               self.vbf_maxMjj_jj_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxMjj_jj_phi_gen");                               

               self.vbf_maxMjj_j1_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxMjj_j1_m_gen");                               
               self.vbf_maxMjj_j1_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxMjj_j1_pt_gen");                               
               self.vbf_maxMjj_j1_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j1_eta_gen");                               
               self.vbf_maxMjj_j1_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j1_phi_gen");                               
 
               self.vbf_maxMjj_j2_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxMjj_j2_m_gen");                               
               self.vbf_maxMjj_j2_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxMjj_j2_pt_gen");                               
               self.vbf_maxMjj_j2_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j2_eta_gen");                               
               self.vbf_maxMjj_j2_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxMjj_j2_phi_gen");                               
 
               self.vbf_maxMjj_j1_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j1_bDiscriminatorCSV_gen");
               self.vbf_maxMjj_j2_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxMjj_j2_bDiscriminatorCSV_gen");                                
              ### scaling jets
              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxMjj_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxMjj_j1_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxMjj_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxMjj_j1_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ; 

              self.vbf_maxMjj_j1 = ROOT.TLorentzVector();
              self.vbf_maxMjj_j1.SetPtEtaPhiM(self.vbf_maxMjj_j1_pt_[0],self.vbf_maxMjj_j1_eta_[0],self.vbf_maxMjj_j1_phi_[0],self.vbf_maxMjj_j1_m_[0]); 
              vbf_maxMjj_j1_up = ROOT.TLorentzVector(self.vbf_maxMjj_j1.Px()*(1+j_jecfactorAK5_up_[0]),
                                                     self.vbf_maxMjj_j1.Py()*(1+j_jecfactorAK5_up_[0]),
                                                     self.vbf_maxMjj_j1.Pz()*(1+j_jecfactorAK5_up_[0]),
                                                     self.vbf_maxMjj_j1.E()*(1+j_jecfactorAK5_up_[0]));
                                                   

              self.vbf_maxMjj_j1_m_up_[0]   = vbf_maxMjj_j1_up.M();                                
              self.vbf_maxMjj_j1_pt_up_[0]  = vbf_maxMjj_j1_up.Pt();                              
              self.vbf_maxMjj_j1_eta_up_[0] = vbf_maxMjj_j1_up.Eta();                                
              self.vbf_maxMjj_j1_phi_up_[0] = vbf_maxMjj_j1_up.Phi();                                
 
              vbf_maxMjj_j1_dn = ROOT.TLorentzVector(self.vbf_maxMjj_j1.Px()*(1+j_jecfactorAK5_dn_[0]),
                                                     self.vbf_maxMjj_j1.Py()*(1+j_jecfactorAK5_dn_[0]),
                                                     self.vbf_maxMjj_j1.Pz()*(1+j_jecfactorAK5_dn_[0]),
                                                     self.vbf_maxMjj_j1.E()*(1+j_jecfactorAK5_dn_[0]));
                                                   
              self.vbf_maxMjj_j1_m_dn_[0]   = vbf_maxMjj_j1_dn.M();                                
              self.vbf_maxMjj_j1_pt_dn_[0]  = vbf_maxMjj_j1_dn.Pt();                              
              self.vbf_maxMjj_j1_eta_dn_[0] = vbf_maxMjj_j1_dn.Eta();                                

              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxMjj_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxMjj_j2_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxMjj_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxMjj_j2_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              self.vbf_maxMjj_j2 = ROOT.TLorentzVector();
              self.vbf_maxMjj_j2.SetPtEtaPhiM(self.vbf_maxMjj_j2_pt_[0],self.vbf_maxMjj_j2_eta_[0],self.vbf_maxMjj_j2_phi_[0],self.vbf_maxMjj_j2_m_[0]); 
              vbf_maxMjj_j2_up = ROOT.TLorentzVector(self.vbf_maxMjj_j2.Px()*(1+j_jecfactorAK5_up_[0]),
                                                     self.vbf_maxMjj_j2.Py()*(1+j_jecfactorAK5_up_[0]),
                                                     self.vbf_maxMjj_j2.Pz()*(1+j_jecfactorAK5_up_[0]),
                                                     self.vbf_maxMjj_j2.E()*(1+j_jecfactorAK5_up_[0]));
                                                   

              self.vbf_maxMjj_j2_m_up_[0]   = vbf_maxMjj_j2_up.M();                                
              self.vbf_maxMjj_j2_pt_up_[0]  = vbf_maxMjj_j2_up.Pt();                              
              self.vbf_maxMjj_j2_eta_up_[0] = vbf_maxMjj_j2_up.Eta();                                
              self.vbf_maxMjj_j2_phi_up_[0] = vbf_maxMjj_j2_up.Phi();                                

              vbf_maxMjj_j2_dn = ROOT.TLorentzVector(self.vbf_maxMjj_j2.Px()*(1+j_jecfactorAK5_dn_[0]),
                                                     self.vbf_maxMjj_j2.Py()*(1+j_jecfactorAK5_dn_[0]),
                                                     self.vbf_maxMjj_j2.Pz()*(1+j_jecfactorAK5_dn_[0]),
                                                     self.vbf_maxMjj_j2.E()*(1+j_jecfactorAK5_dn_[0]));
                                                   

              self.vbf_maxMjj_j2_m_dn_[0]   = vbf_maxMjj_j2_dn.M();                                
              self.vbf_maxMjj_j2_pt_dn_[0]  = vbf_maxMjj_j2_dn.Pt();                              
              self.vbf_maxMjj_j2_eta_dn_[0] = vbf_maxMjj_j2_dn.Eta();                                


             ######### btag stuff  

             dR_lj = 0. ;
             index_ak5_cvst = array( 'f', [0.] );
             index_ak5_cvsm = array( 'f', [0.] );
             index_ak5_cvsl = array( 'f', [0.] );

             self.nbjets_csvt_veto_cleaned_[0] = 0. ;
             self.nbjets_csvm_veto_cleaned_[0] = 0. ;
             self.nbjets_csvl_veto_cleaned_[0] = 0. ;
             self.nbjets_csvl_veto_[0] = 0. ;
             self.nbjets_csvm_veto_[0] = 0. ;
             self.nbjets_csvt_veto_[0] = 0. ;

                 

             ## take the delta R between leading CA8 jet and the lepton 
             for i in range(6):
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[i] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[i] > 30:
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[i]
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[i]
                                                                                
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" )                
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" )                
                        dR_jj = math.sqrt( (j_ak5_eta - getattr(self.InputTree_, prefix+"_eta")[0])**2 + (j_ak5_phi - getattr(self.InputTree_, prefix+"_phi")[0])**2 );
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + (l_phi - j_ak5_phi)**2 );
                        if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.244: self.nbjets_csvl_veto_[0] = self.nbjets_csvl_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.679: self.nbjets_csvm_veto_[0] = self.nbjets_csvm_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.898: self.nbjets_csvt_veto_[0] = self.nbjets_csvt_veto_[0]+1;
                        
                        if dR_jj > 0.8: 
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.244: index_ak5_cvsl.append(i)
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.679: index_ak5_cvsm.append(i)
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_cvst.append(i);

             for i in range(6):
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[i] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[i] > 30:
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCorVBFTag_Eta" )[i]
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCorVBFTag_Phi" )[i]
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" )                
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" )                
                        dR_jj = math.sqrt( (j_ak5_eta - getattr(self.InputTree_, prefix+"_eta")[0])**2 + (j_ak5_phi - getattr(self.InputTree_, prefix+"_phi")[0])**2 );
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + (l_phi - j_ak5_phi)**2 );
                        if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] >=0.244: self.nbjets_csvl_veto_[0] = self.nbjets_csvl_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] >=0.679: self.nbjets_csvm_veto_[0] = self.nbjets_csvm_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[i] >=0.898: self.nbjets_csvt_veto_[0] = self.nbjets_csvt_veto_[0]+1;
                        if dR_jj > 0.8: 
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.244: index_ak5_cvsl.append(i)
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.679: index_ak5_cvsm.append(i)
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_cvst.append(i);

               
             self.nbjets_csvl_veto_cleaned_[0] = len(index_ak5_in_oppoHemi_csvl);
             self.nbjets_csvm_veto_cleaned_[0] = len(index_ak5_in_oppoHemi_csvm);
             self.nbjets_csvt_veto_cleaned_[0] = len(index_ak5_in_oppoHemi_csvt);

             self.nbjets_ssvhem_veto_  = getattr( self.InputTree_, "numPFCorJetBTags");

             self.njets_[0] = getattr( self.InputTree_, "numPFCorJets" ) + getattr( self.InputTree_, "numPFCorVBFTagJets" ) ;
               
             listOfVarArray1[0][0] = self.jet_qjetvol_[0];
             listOfVarArray1[1][0] = self.jet_tau2tau1_[0];

             listOfVarArray2[0][0] = self.jet_qjetvol_[0];
             listOfVarArray2[1][0] = self.jet_tau2tau1_[0];

             ### Fill the output tree

             self.otree.Fill();

        ### close the file 
        self.OFile_.cd();
        self.otree.Write();
        self.OFile_.Close();

    def InitializeVariables(self):

        ### event Property and weights

        self.nPV_               = array( 'f', [ 0. ] );                        
        self.event_runNo_       = array( 'i', [ 0 ] );
        self.event_lumi_        = array( 'i', [ 0 ] );
        self.event_             = array( 'i', [0] );

        self.totalEventWeight_  = array( 'f', [ 0. ] );                                
        self.eff_and_pu_Weight_ = array( 'f', [ 0. ] );                                
        self.wSampleWeight_     = array( 'f', [ 0. ] );
        self.event_weight_      = array( 'f', [ 0. ] );
        
        self.btag_weight_ = array( 'f', [0.]); 
        self.btag_weight_up_ = array( 'f', [0.]); 
        self.btag_weight_dn_ = array( 'f', [0.]); 
        self.btag_weight_up_dn_ = array( 'f', [0.]); 
        self.btag_weight_dn_up_ = array( 'f', [0.]); 

        self.interference_Weight_H600_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H700_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H800_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H900_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H1000_ = array( 'f', [ 0. ] );                                

        self.cps_Weight_H600_ = array( 'f', [ 0. ] );                                
        self.cps_Weight_H700_ = array( 'f', [ 0. ] );                                
        self.cps_Weight_H800_ = array( 'f', [ 0. ] );                                
        self.cps_Weight_H900_ = array( 'f', [ 0. ] );                                
        self.cps_Weight_H1000_ = array( 'f', [ 0. ] );                                

        ########## event topology variables

        self.issignal_     = array( 'i', [ 0 ] );             
        self.numberJetBin_ = array( 'i', [ 0 ] );
        self.numberJetBin2_ = array( 'i', [ 0 ] );
        self.numberJetBin2_ = array( 'i', [ 0 ] );
        self.numberJetBin3_ = array( 'i', [ 0 ] );
        self.numberJetBin4_ = array( 'i', [ 0 ] );

        self.numberJetBinGen_ = array( 'i', [ 0 ] );
        self.numberJetBinGen2_ = array( 'i', [ 0 ] );
        self.numberJetBinGen2_ = array( 'i', [ 0 ] );
        self.numberJetBinGen3_ = array( 'i', [ 0 ] );
        self.numberJetBinGen4_ = array( 'i', [ 0 ] );

        ### Leptonic W, Lepton and Nueutrino
        
        self.l_pt_     = array( 'f', [ 0. ] );
        self.l_eta_    = array( 'f', [ 0. ] );
        self.l_phi_    = array( 'f', [ 0. ] );
        self.l_charge_ = array( 'f', [ 0. ] );
        
        self.pfMET_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_ = array( 'f', [ 0. ] );

        self.nu_pz_type0_ = array( 'f', [ 0. ] );
        self.nu_pz_type2_ = array( 'f', [ 0. ] );

        self.nu_pz_type0_met_ = array( 'f', [ 0. ] );
        self.nu_pz_type2_met_ = array( 'f', [ 0. ] );

        self.W_pz_type0_ = array( 'f', [ 0. ] );
        self.W_pz_type2_ = array( 'f', [ 0. ] );

        self.W_pz_type0_met_ = array( 'f', [ 0. ] );
        self.W_pz_type2_met_ = array( 'f', [ 0. ] );

        self.nu_pz_gen_   = array( 'f', [ 0. ] );
        self.W_pz_gen_    = array( 'f', [ 0. ] );
        self.W_pt_gen_    = array( 'f', [ 0. ] );


        ###### new variables for tree branches         
        self.mass_lvj_type0_met_    = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_    = array( 'f', [ 0. ] );

        self.mass_lvj_type0_    = array( 'f', [ 0. ] );
        self.mass_lvj_type2_    = array( 'f', [ 0. ] );

        self.mass_lv_subj_type0_met_    = array( 'f', [ 0. ] );
        self.mass_lv_subj_type2_met_    = array( 'f', [ 0. ] );

        self.mass_lv_subj_type0_    = array( 'f', [ 0. ] );
        self.mass_lv_subj_type2_    = array( 'f', [ 0. ] );

        ### leptonic W transverse mass and pt 
        self.v_pt_             = array( 'f', [ 0. ] );
        self.v_mt_             = array( 'f', [ 0. ] );
        self.v_eta_            = array( 'f', [ 0. ] );
        self.v_phi_            = array( 'f', [ 0. ] );

        ######### Hadronic W Variables ###############

        self.ungroomed_jet_eta_ = array( 'f', [ 0. ] );
        self.ungroomed_jet_pt_  = array( 'f', [ 0. ] );
        self.ungroomed_jet_phi_ = array( 'f', [ 0. ] );
        self.ungroomed_jet_e_   = array( 'f', [ 0. ] );

        self.ungroomed_gen_jet_pt_  = array( 'f', [ 0. ] );
        self.ungroomed_gen_jet_eta_ = array( 'f', [ 0. ] );
        self.ungroomed_gen_jet_phi_ = array( 'f', [ 0. ] );
        self.ungroomed_gen_jet_e_   = array( 'f', [ 0. ] );

        self.jet_mass_pr_       = array( 'f', [ 0. ] );
        self.jet_pt_pr_         = array( 'f', [ 0. ] );
        self.jet_charge_        = array( 'f', [ 0. ] );
        self.jet_charge_k05_    = array( 'f', [ 0. ] );
        self.jet_charge_k07_    = array( 'f', [ 0. ] );
        self.jet_charge_k10_    = array( 'f', [ 0. ] );

        self.gen_jet_mass_pr_       = array( 'f', [ 0. ] );
        self.gen_jet_pt_pr_         = array( 'f', [ 0. ] );
        self.gen_jet_charge_        = array( 'f', [ 0. ] );
        self.gen_jet_charge_k05_    = array( 'f', [ 0. ] );
        self.gen_jet_charge_k07_    = array( 'f', [ 0. ] );
        self.gen_jet_charge_k10_    = array( 'f', [ 0. ] );

        self.jet_grsens_ft_   = array( 'f', [ 0. ] );
        self.jet_grsens_tr_   = array( 'f', [ 0. ] );
        self.jet_massdrop_pr_ = array( 'f', [ 0. ] );    
        self.jet_qjetvol_     = array( 'f', [ 0. ] ); 

        self.gen_jet_grsens_ft_   = array( 'f', [ 0. ] );
        self.gen_jet_grsens_tr_   = array( 'f', [ 0. ] );
        self.gen_jet_massdrop_pr_ = array( 'f', [ 0. ] );    
        self.gen_jet_qjetvol_     = array( 'f', [ 0. ] ); 

        self.jet_tau2tau1_ = array( 'f', [ 0. ] );     
        self.jet_tau2tau1_exkT_ = array( 'f', [ 0. ] );     
        self.jet_tau2tau1_pr_ = array( 'f', [ 0. ] );
        self.jet_GeneralizedECF_ = array( 'f', [ 0. ] );

        self.gen_jet_tau2tau1_ = array( 'f', [ 0. ] );     
        self.gen_jet_tau2tau1_exkT_ = array( 'f', [ 0. ] );     
        self.gen_jet_tau2tau1_pr_ = array( 'f', [ 0. ] );
        self.gen_jet_GeneralizedECF_ = array( 'f', [ 0. ] );

        self.jet_jetconstituents_  = array( 'f', [ 0. ] );     
        self.gen_jet_jetconstituents_  = array( 'f', [ 0. ] );     
        
        self.jet_rcore4_ = array( 'f', [ 0. ] );
        self.jet_rcore5_ = array( 'f', [ 0. ] );
        self.jet_rcore6_ = array( 'f', [ 0. ] );
        self.jet_rcore7_ = array( 'f', [ 0. ] );

        self.gen_jet_rcore4_ = array( 'f', [ 0. ] );
        self.gen_jet_rcore5_ = array( 'f', [ 0. ] );
        self.gen_jet_rcore6_ = array( 'f', [ 0. ] );
        self.gen_jet_rcore7_ = array( 'f', [ 0. ] );


        self.jet_pt1frac_  = array ('f',[ 0. ]);
        self.jet_pt2frac_  = array ('f',[ 0. ]);
        self.jet_sjdr_     = array ('f',[ 0. ]);
        

        self.gen_jet_pt1frac_  = array ('f',[ 0. ]);
        self.gen_jet_pt2frac_  = array ('f',[ 0. ]);
        self.gen_jet_sjdr_     = array ('f',[ 0. ]);
        

        self.j_jecfactor_up_ = array( 'f', [ 0. ] );
        self.j_jecfactor_dn_ = array( 'f', [ 0. ] );

        self.jet_mass_pr_up_ = array( 'f', [ 0. ] );
        self.jet_mass_pr_dn_ = array( 'f', [ 0. ] );

        self.jet_ungroomed_jet_pt_up_ = array( 'f', [ 0. ] );
        self.jet_ungroomed_jet_pt_dn_ = array( 'f', [ 0. ] );

        self.jet_planarlow04_ = array( 'f', [0.] );
        self.jet_planarlow05_ = array( 'f', [0.] );
        self.jet_planarlow06_ = array( 'f', [0.] );
        self.jet_planarlow07_ = array( 'f', [0.] );

        self.gen_jet_planarlow04_ = array( 'f', [0.] );
        self.gen_jet_planarlow05_ = array( 'f', [0.] );
        self.gen_jet_planarlow06_ = array( 'f', [0.] );
        self.gen_jet_planarlow07_ = array( 'f', [0.] );

        ###### vbf variables
        self.vbf_maxpt_jj_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_jj_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_QGLikelihood_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_QGLikelihood_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxpt_j1_isPileUpTight_  = array( 'i', [ 0 ] );                                
        self.vbf_maxpt_j2_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxpt_j2_isPileUpTight_  = array( 'i', [ 0 ] );                                

        self.vbf_maxpt_j1_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_jj_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_jj_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_jj_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_jj_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_jj_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_jj_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_jj_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_jj_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_m_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_pt_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_eta_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_phi_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_m_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_pt_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_eta_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j1_phi_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j2_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j2_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j2_m_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_pt_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_eta_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_phi_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j2_m_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_pt_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_eta_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_phi_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_QGLikelihood_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_QGLikelihood_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxDeta_j1_isPileUpTight_  = array( 'i', [ 0 ] );                                
        self.vbf_maxDeta_j2_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxDeta_j2_isPileUpTight_  = array( 'i', [ 0 ] );                                

        self.vbf_maxDeta_j1_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                

        self.vbf_maxDeta_j1_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxDeta_j2_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_jj_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_jj_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_jj_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_jj_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_jj_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_jj_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_jj_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_jj_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_m_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_pt_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_eta_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_phi_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_m_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_pt_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_eta_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j1_phi_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j2_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j2_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j2_m_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_pt_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_eta_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_phi_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j2_m_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_pt_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_eta_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_phi_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_QGLikelihood_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_QGLikelihood_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxMjj_j1_isPileUpTight_  = array( 'i', [ 0 ] );                                
        self.vbf_maxMjj_j2_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxMjj_j2_isPileUpTight_  = array( 'i', [ 0 ] );                                

        self.vbf_maxMjj_j1_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                

        self.vbf_maxMjj_j1_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxMjj_j2_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                

        ###### btag counters
        self.nbjets_csvl_veto_ = array( 'f', [ 0. ] );
        self.nbjets_csvm_veto_ = array( 'f', [ 0. ] );
        self.nbjets_csvt_veto_ = array( 'f', [ 0. ] );
        self.nbjets_ssvhem_veto_ = array( 'f', [ 0. ] );

        self.nbjets_csvl_veto_cleaned_   = array( 'f', [ 0. ] );
        self.nbjets_csvm_veto_cleaned_   = array( 'f', [ 0. ] );
        self.nbjets_csvt_veto_cleaned_   = array( 'f', [ 0. ] );
        self.nbjets_ssvhem_veto_cleaned_ = array( 'f', [ 0. ] );

        self.njets_ = array( 'f', [ 0. ] );

        self.deltaR_lca8jet_          = array( 'f', [0.] );
        self.deltaphi_METca8jet_      = array( 'f', [0.] );
        self.deltaphi_Vca8jet_        = array( 'f', [0.] );
        self.deltaphi_METca8jet_met_  = array( 'f', [0.] );
        self.deltaphi_Vca8jet_met_    = array( 'f', [0.] );

        ### some generator level quantites
        self.genHMass_ = array( 'f', [ 0. ] );              
        self.genHphi_  = array( 'f', [ 0. ] );              
        self.genHeta_  = array( 'f', [ 0. ] );              
        self.genHpt_   = array( 'f', [ 0. ] );              

        self.genTagQuark1E_    = array( 'f', [ 0. ] );              
        self.genTagQuark1phi_  = array( 'f', [ 0. ] );              
        self.genTagQuark1eta_  = array( 'f', [ 0. ] );              
        self.genTagQuark1pt_   = array( 'f', [ 0. ] );              

        self.genTagQuark2E_    = array( 'f', [ 0. ] );              
        self.genTagQuark2phi_  = array( 'f', [ 0. ] );              
        self.genTagQuark2eta_  = array( 'f', [ 0. ] );              
        self.genTagQuark2pt_   = array( 'f', [ 0. ] );              
                                                                        
        # variables for ttbar control region
        self.ttb_nak5_same_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_same_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_same_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_same_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_nak5_oppo_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppo_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppo_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppo_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_nak5_oppoveto_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppoveto_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppoveto_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppoveto_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_nak5_sameveto_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_sameveto_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_sameveto_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_sameveto_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_ht_           = array( 'f', [ 0. ] );                
        self.ttb_ca8_mass_pr_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_charge_   = array( 'f', [ 0. ] );
        
        self.ttb_ca8_charge_k05_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_charge_k07_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_charge_k10_  = array( 'f', [ 0. ] );                        

        self.ttb_ca8_px_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_py_  = array( 'f', [ 0. ] );                                
        self.ttb_ca8_pz_  = array( 'f', [ 0. ] );                                
        self.ttb_ca8_e_   = array( 'f', [ 0. ] );                                
        
        self.ttb_ca8_ungroomed_pt_   = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_eta_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_phi_  = array( 'f', [ 0. ] );
        self.ttb_ca8_ungroomed_e_    = array( 'f', [ 0. ] );

        self.ttb_ca8_ungroomed_gen_pt_   = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_gen_eta_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_gen_phi_  = array( 'f', [ 0. ] );
        self.ttb_ca8_ungroomed_gen_e_    = array( 'f', [ 0. ] );

        self.ttb_ca8_tau2tau1_        = array( 'f', [ 0. ] );                        
        self.ttb_ca8_tau2tau1_exkT_   = array( 'f', [ 0. ] );                        
        self.ttb_ca8_tau2tau1_pr_     = array( 'f', [ 0. ] );                        
        self.ttb_ca8_GeneralizedECF_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_mu_  = array( 'f', [ 0. ] );                        

        self.ttb_ca8_mlvj_type0_    = array( 'f', [ 0. ] );                        
        self.ttb_ca8_mlvj_type2_    = array( 'f', [ 0. ] );                        

        self.ttb_ca8_mlvj_type0_met_    = array( 'f', [ 0. ] );                        
        self.ttb_ca8_mlvj_type2_met_    = array( 'f', [ 0. ] );                        

        self.isttbar_ = array( 'i', [ 0 ] );            

        self.gen_parton1_px_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton1_py_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton1_pz_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton1_e_fromttbar_  = array( 'f', [ 0. ] );
        self.gen_parton1_id_fromttbar_ = array( 'i', [ 0 ] );
        self.gen_parton2_px_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton2_py_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton2_pz_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton2_e_fromttbar_  = array( 'f', [ 0. ] );
        self.gen_parton2_id_fromttbar_ = array( 'i', [ 0 ] );

    def createBranches(self):

        ################ Branch output Tree

        self.otree.Branch("nPV", self.nPV_,"nPV/F");
        self.otree.Branch("event", self.event_ , "event/I");
        self.otree.Branch("event_lumi", self.event_lumi_ , "event_lumi/I");
        self.otree.Branch("event_runNo", self.event_runNo_ , "event_runNo/I");

        self.otree.Branch("issignal", self.issignal_ , "issignal/I");
        self.otree.Branch("numberJetBin", self.numberJetBin_ , "numberJetBin/I");
        self.otree.Branch("numberJetBin2", self.numberJetBin2_ , "numberJetBin2/I");
        self.otree.Branch("numberJetBin3", self.numberJetBin3_ , "numberJetBin3/I");
        self.otree.Branch("numberJetBin4", self.numberJetBin4_ , "numberJetBin4/I");

        self.otree.Branch("numberJetBinGen", self.numberJetBinGen_ , "numberJetBinGen/I");
        self.otree.Branch("numberJetBinGen2", self.numberJetBinGen2_ , "numberJetBinGen2/I");
        self.otree.Branch("numberJetBinGen3", self.numberJetBinGen3_ , "numberJetBinGen3/I");
        self.otree.Branch("numberJetBinGen4", self.numberJetBinGen4_ , "numberJetBinGen4/I");

        self.otree.Branch("totalEventWeight",  self.totalEventWeight_ , "totalEventWeight/F"); ## total xs * pile Up * trigger * lepton ID
        self.otree.Branch("eff_and_pu_Weight", self.eff_and_pu_Weight_ , "eff_and_pu_Weight/F"); ## product of pileUp, trigger and lepton ID
        self.otree.Branch("wSampleWeight",     self.wSampleWeight_ , "wSampleWeight/F"); ## only cross section weight
        self.otree.Branch("event_weight",      self.event_weight_ , "event_weight/F"); ## in case the MC produce weighted events at LHE level

        self.otree.Branch("btag_weight", self.btag_weight_ , "btag_weight/F"); ## btag weight according to btag pog recipe
        self.otree.Branch("btag_weight_up", self.btag_weight_up_ , "btag_weight_up/F"); 
        self.otree.Branch("btag_weight_dn", self.btag_weight_dn_ , "btag_weight_dn/F");
        self.otree.Branch("btag_weight_dn_up", self.btag_weight_dn_up_ , "btag_weight_dn_up/F");
        self.otree.Branch("btag_weight_up_dn", self.btag_weight_up_dn_ , "btag_weight_up_dn/F");

        self.otree.Branch("interference_Weight_H600", self.interference_Weight_H600_ , "interference_Weight_H600/F");
        self.otree.Branch("interference_Weight_H700", self.interference_Weight_H700_ , "interference_Weight_H700/F");
        self.otree.Branch("interference_Weight_H800", self.interference_Weight_H800_ , "interference_Weight_H800/F");
        self.otree.Branch("interference_Weight_H900", self.interference_Weight_H900_ , "interference_Weight_H900/F");
        self.otree.Branch("interference_Weight_H1000", self.interference_Weight_H1000_ , "interference_Weight_H1000/F");

        self.otree.Branch("cps_Weight_H600", self.cps_Weight_H600_ , "cps_Weight_H600/F");
        self.otree.Branch("cps_Weight_H700", self.cps_Weight_H700_ , "cps_Weight_H700/F");
        self.otree.Branch("cps_Weight_H800", self.cps_Weight_H800_ , "cps_Weight_H800/F");
        self.otree.Branch("cps_Weight_H900", self.cps_Weight_H900_ , "cps_Weight_H900/F");
        self.otree.Branch("cps_Weight_H1000", self.cps_Weight_H1000_ , "cps_Weight_H1000/F");

        ##############
        
        self.otree.Branch("mass_lvj_type0_met", self.mass_lvj_type0_met_ , "mass_lvj_type0_met/F");
        self.otree.Branch("mass_lvj_type2_met", self.mass_lvj_type2_met_ , "mass_lvj_type2_met/F");

        self.otree.Branch("mass_lvj_type0", self.mass_lvj_type0_ , "mass_lvj_type0/F");
        self.otree.Branch("mass_lvj_type2", self.mass_lvj_type2_ , "mass_lvj_type2/F");

        self.otree.Branch("mass_lv_subj_type0_met", self.mass_lv_subj_type0_met_ , "mass_lv_subj_type0_met/F");
        self.otree.Branch("mass_lv_subj_type2_met", self.mass_lv_subj_type2_met_ , "mass_lv_subj_type2_met/F");

        self.otree.Branch("mass_lv_subj_type0", self.mass_lv_subj_type0_ , "mass_lv_subj_type0/F");
        self.otree.Branch("mass_lv_subj_type2", self.mass_lv_subj_type2_ , "mass_lv_subj_type2/F");

        ########### Leptonic W, Lepton and Nueutrino
        
        self.otree.Branch("l_pt", self.l_pt_ , "l_pt/F");
        self.otree.Branch("l_eta", self.l_eta_ , "l_eta/F");
        self.otree.Branch("l_charge", self.l_charge_ , "l_charge/F");
        self.otree.Branch("l_phi", self.l_phi_ , "l_phi/F");

        self.otree.Branch("pfMET", self.pfMET_ , "pfMET/F");
        self.otree.Branch("pfMET_Phi", self.pfMET_Phi_ , "pfMET_Phi/F");

        self.otree.Branch("nu_pz_type0", self.nu_pz_type0_ , "nu_pz_type0/F");
        self.otree.Branch("nu_pz_type2", self.nu_pz_type2_ , "nu_pz_type2/F");

        self.otree.Branch("nu_pz_type0_met", self.nu_pz_type0_met_ , "nu_pz_type0_met/F");
        self.otree.Branch("nu_pz_type2_met", self.nu_pz_type2_met_ , "nu_pz_type2_met/F");

        self.otree.Branch("W_pz_type0", self.W_pz_type0_ , "W_pz_type0/F");
        self.otree.Branch("W_pz_type2", self.W_pz_type2_ , "W_pz_type2/F");

        self.otree.Branch("W_pz_type0_met", self.W_pz_type0_met_ , "W_pz_type0_met/F");
        self.otree.Branch("W_pz_type2_met", self.W_pz_type2_met_ , "W_pz_type2_met/F");

        self.otree.Branch("nu_pz_gen", self.nu_pz_gen_ , "nu_pz_gen/F");
        self.otree.Branch("W_pz_gen", self.W_pz_gen_ , "W_pz_gen/F");
        self.otree.Branch("W_pt_gen", self.W_pt_gen_ , "W_pt_gen/F");

        self.otree.Branch("v_pt", self.v_pt_ , "v_pt/F");
        self.otree.Branch("v_mt", self.v_mt_ , "v_mt/F");
        self.otree.Branch("v_eta", self. v_eta_ , "v_eta/F");
        self.otree.Branch("v_phi", self.v_phi_ , "v_phi/F");

        ###################

        self.otree.Branch("ungroomed_jet_eta", self.ungroomed_jet_eta_ , "ungroomed_jet_eta/F");
        self.otree.Branch("ungroomed_jet_phi", self.ungroomed_jet_phi_ , "ungroomed_jet_phi/F");
        self.otree.Branch("ungroomed_jet_pt", self.ungroomed_jet_pt_ , "ungroomed_jet_pt/F");
        self.otree.Branch("ungroomed_jet_e", self.ungroomed_jet_e_ , "ungroomed_jet_e/F");

        self.otree.Branch("ungroomed_gen_jet_eta", self.ungroomed_gen_jet_eta_ , "ungroomed_gen_jet_eta/F");
        self.otree.Branch("ungroomed_gen_jet_phi", self.ungroomed_gen_jet_phi_ , "ungroomed_gen_jet_phi/F");
        self.otree.Branch("ungroomed_gen_jet_pt", self.ungroomed_gen_jet_pt_ , "ungroomed_gen_jet_pt/F");
        self.otree.Branch("ungroomed_gen_jet_e", self.ungroomed_gen_jet_e_ , "ungroomed_gen_jet_e/F");

        self.otree.Branch("jet_mass_pr", self.jet_mass_pr_ , "jet_mass_pr/F");
        self.otree.Branch("jet_pt_pr", self.jet_pt_pr_ , "jet_pt_pr/F");
        self.otree.Branch("jet_charge", self.jet_charge_ , "jet_charge/F");
        self.otree.Branch("jet_charge_k05", self.jet_charge_k05_ , "jet_charge_k05/F");
        self.otree.Branch("jet_charge_k07", self.jet_charge_k07_ , "jet_charge_k07/F");
        self.otree.Branch("jet_charge_k10", self.jet_charge_k10_ , "jet_charge_k10/F");

        self.otree.Branch("gen_jet_mass_pr", self.gen_jet_mass_pr_ , "gen_jet_mass_pr/F");
        self.otree.Branch("gen_jet_pt_pr", self.gen_jet_pt_pr_ , "gen_jet_pt_pr/F");
        self.otree.Branch("gen_jet_charge", self.gen_jet_charge_ , "gen_jet_charge/F");
        self.otree.Branch("gen_jet_charge_k05", self.gen_jet_charge_k05_ , "gen_jet_charge_k05/F");
        self.otree.Branch("gen_jet_charge_k07", self.gen_jet_charge_k07_ , "gen_jet_charge_k07/F");
        self.otree.Branch("gen_jet_charge_k10", self.gen_jet_charge_k10_ , "gen_jet_charge_k10/F");

        self.otree.Branch("jet_grsens_ft", self.jet_grsens_ft_ , "jet_grsens_ft/F");
        self.otree.Branch("jet_grsens_tr", self.jet_grsens_tr_ , "jet_grsens_tr/F");
        self.otree.Branch("jet_massdrop_pr", self.jet_massdrop_pr_ , "jet_massdrop_pr/F");
        self.otree.Branch("jet_qjetvol", self.jet_qjetvol_ , "jet_qjetvol/F");

        self.otree.Branch("gen_jet_grsens_ft", self.gen_jet_grsens_ft_ , "gen_jet_grsens_ft/F");
        self.otree.Branch("gen_jet_grsens_tr", self.gen_jet_grsens_tr_ , "gen_jet_grsens_tr/F");
        self.otree.Branch("gen_jet_massdrop_pr", self.gen_jet_massdrop_pr_ , "gen_jet_massdrop_pr/F");
        self.otree.Branch("gen_jet_qjetvol", self.gen_jet_qjetvol_ , "gen_jet_qjetvol/F");

        self.otree.Branch("jet_tau2tau1", self.jet_tau2tau1_ , "jet_tau2tau1/F");
        self.otree.Branch("jet_tau2tau1_exkT", self.jet_tau2tau1_exkT_ , "jet_tau2tau1_exkT/F");
        self.otree.Branch("jet_tau2tau1_pr", self.jet_tau2tau1_pr_ , "jet_tau2tau1_pr/F");
        self.otree.Branch("jet_GeneralizedECF", self.jet_GeneralizedECF_ , "jet_GeneralizedECF/F");

        self.otree.Branch("gen_jet_tau2tau1", self.gen_jet_tau2tau1_ , "gen_jet_tau2tau1/F");
        self.otree.Branch("gen_jet_tau2tau1_exkT", self.gen_jet_tau2tau1_exkT_ , "gen_jet_tau2tau1_exkT/F");
        self.otree.Branch("gen_jet_tau2tau1_pr", self.gen_jet_tau2tau1_pr_ , "gen_jet_tau2tau1_pr/F");
        self.otree.Branch("gen_jet_GeneralizedECF", self.gen_jet_GeneralizedECF_ , "gen_jet_GeneralizedECF/F");

        self.otree.Branch("jet_jetconstituents", self.jet_jetconstituents_ , "jet_jetconstituents/F");
        self.otree.Branch("gen_jet_jetconstituents", self.gen_jet_jetconstituents_ , "gen_jet_jetconstituents/F");

        self.otree.Branch("jet_rcore4", self.jet_rcore4_ , "jet_rcore4/F");
        self.otree.Branch("jet_rcore5", self.jet_rcore5_ , "jet_rcore5/F");
        self.otree.Branch("jet_rcore6", self.jet_rcore6_ , "jet_rcore6/F");
        self.otree.Branch("jet_rcore7", self.jet_rcore7_ , "jet_rcore7/F");

        self.otree.Branch("gen_jet_rcore4", self.gen_jet_rcore4_ , "gen_jet_rcore4/F");
        self.otree.Branch("gen_jet_rcore5", self.gen_jet_rcore5_ , "gen_jet_rcore5/F");
        self.otree.Branch("gen_jet_rcore6", self.gen_jet_rcore6_ , "gen_jet_rcore6/F");
        self.otree.Branch("gen_jet_rcore7", self.gen_jet_rcore7_ , "gen_jet_rcore7/F");

        self.otree.Branch("jet_pt1frac", self.jet_pt1frac_ , "jet_pt1frac/F");
        self.otree.Branch("jet_pt2frac", self.jet_pt2frac_ , "jet_pt2frac/F");
        self.otree.Branch("jet_sjdr", self.jet_sjdr_ , "jet_sjdr/F");

        self.otree.Branch("gen_jet_pt1frac", self.gen_jet_pt1frac_ , "gen_jet_pt1frac/F");
        self.otree.Branch("gen_jet_pt2frac", self.gen_jet_pt2frac_ , "gen_jet_pt2frac/F");
        self.otree.Branch("gen_jet_sjdr", self.gen_jet_sjdr_ , "gen_jet_sjdr/F");


        self.otree.Branch("j_jecfactor_up", self.j_jecfactor_up_ , "j_jecfactor_up/F");
        self.otree.Branch("j_jecfactor_dn", self.j_jecfactor_dn_ , "j_jecfactor_dn/F");

        self.otree.Branch("jet_mass_pr_up", self.jet_mass_pr_up_ , "jet_mass_pr_up/F");
        self.otree.Branch("jet_mass_pr_dn", self.jet_mass_pr_dn_ , "jet_mass_pr_dn/F");

        self.otree.Branch("jet_ungroomed_jet_pt_dn", self.jet_ungroomed_jet_pt_dn_ , "jet_ungroomed_jet_pt_dn/F");
        self.otree.Branch("jet_ungroomed_jet_pt_up", self.jet_ungroomed_jet_pt_up_ , "jet_ungroomed_jet_pt_up/F");

        self.otree.Branch("jet_planarlow04", self.jet_planarlow04_ , "jet_planarlow04/F");
        self.otree.Branch("jet_planarlow05", self.jet_planarlow05_ , "jet_planarlow05/F");
        self.otree.Branch("jet_planarlow06", self.jet_planarlow06_ , "jet_planarlow06/F");
        self.otree.Branch("jet_planarlow07", self.jet_planarlow07_ , "jet_planarlow07/F");

        self.otree.Branch("gen_jet_planarlow04", self.gen_jet_planarlow04_ , "gen_jet_planarlow04/F");
        self.otree.Branch("gen_jet_planarlow05", self.gen_jet_planarlow05_ , "gen_jet_planarlow05/F");
        self.otree.Branch("gen_jet_planarlow06", self.gen_jet_planarlow06_ , "gen_jet_planarlow06/F");
        self.otree.Branch("gen_jet_planarlow07", self.gen_jet_planarlow07_ , "gen_jet_planarlow07/F");

        ########## vbf jets 
        
        self.otree.Branch("vbf_maxpt_jj_m", self.vbf_maxpt_jj_m_ , "vbf_maxpt_jj_m/F");
        self.otree.Branch("vbf_maxpt_jj_pt", self.vbf_maxpt_jj_pt_ , "vbf_maxpt_jj_pt/F");
        self.otree.Branch("vbf_maxpt_jj_eta", self.vbf_maxpt_jj_eta_ , "vbf_maxpt_jj_eta/F");
        self.otree.Branch("vbf_maxpt_jj_phi", self.vbf_maxpt_jj_phi_ , "vbf_maxpt_jj_phi/F");

        self.otree.Branch("vbf_maxpt_j1_m", self.vbf_maxpt_j1_m_ , "vbf_maxpt_j1_m/F");
        self.otree.Branch("vbf_maxpt_j1_pt", self.vbf_maxpt_j1_pt_ , "vbf_maxpt_j1_pt/F");
        self.otree.Branch("vbf_maxpt_j1_eta", self.vbf_maxpt_j1_eta_ , "vbf_maxpt_j1_eta/F");
        self.otree.Branch("vbf_maxpt_j1_phi", self.vbf_maxpt_j1_phi_ , "vbf_maxpt_j1_phi/F");

        self.otree.Branch("vbf_maxpt_j2_m", self.vbf_maxpt_j2_m_ , "vbf_maxpt_j2_m/F");
        self.otree.Branch("vbf_maxpt_j2_pt", self.vbf_maxpt_j2_pt_ , "vbf_maxpt_j2_pt/F");
        self.otree.Branch("vbf_maxpt_j2_eta", self.vbf_maxpt_j2_eta_ , "vbf_maxpt_j2_eta/F");
        self.otree.Branch("vbf_maxpt_j2_phi", self.vbf_maxpt_j2_phi_ , "vbf_maxpt_j2_phi/F");

        self.otree.Branch("vbf_maxpt_jj_m_gen", self.vbf_maxpt_jj_m_gen_ , "vbf_maxpt_jj_m_gen/F");
        self.otree.Branch("vbf_maxpt_jj_pt_gen", self.vbf_maxpt_jj_pt_gen_ , "vbf_maxpt_jj_pt_gen/F");
        self.otree.Branch("vbf_maxpt_jj_eta_gen", self.vbf_maxpt_jj_eta_gen_ , "vbf_maxpt_jj_eta_gen/F");
        self.otree.Branch("vbf_maxpt_jj_phi_gen", self.vbf_maxpt_jj_phi_gen_ , "vbf_maxpt_jj_phi_gen/F");

        self.otree.Branch("vbf_maxpt_j1_m_gen", self.vbf_maxpt_j1_m_gen_ , "vbf_maxpt_j1_m_gen/F");
        self.otree.Branch("vbf_maxpt_j1_pt_gen", self.vbf_maxpt_j1_pt_gen_ , "vbf_maxpt_j1_pt_gen/F");
        self.otree.Branch("vbf_maxpt_j1_eta_gen", self.vbf_maxpt_j1_eta_gen_ , "vbf_maxpt_j1_eta_gen/F");
        self.otree.Branch("vbf_maxpt_j1_phi_gen", self.vbf_maxpt_j1_phi_gen_ , "vbf_maxpt_j1_phi_gen/F");

        self.otree.Branch("vbf_maxpt_j2_m_gen", self.vbf_maxpt_j2_m_gen_ , "vbf_maxpt_j2_m_gen/F");
        self.otree.Branch("vbf_maxpt_j2_pt_gen", self.vbf_maxpt_j2_pt_gen_ , "vbf_maxpt_j2_pt_gen/F");
        self.otree.Branch("vbf_maxpt_j2_eta_gen", self.vbf_maxpt_j2_eta_gen_ , "vbf_maxpt_j2_eta_gen/F");
        self.otree.Branch("vbf_maxpt_j2_phi_gen", self.vbf_maxpt_j2_phi_gen_ , "vbf_maxpt_j2_phi_gen/F");

        self.otree.Branch("vbf_maxpt_j1_m_up", self.vbf_maxpt_j1_m_up_ , "vbf_maxpt_j1_m_up/F");
        self.otree.Branch("vbf_maxpt_j1_pt_up", self.vbf_maxpt_j1_pt_up_ , "vbf_maxpt_j1_pt_up/F");
        self.otree.Branch("vbf_maxpt_j1_eta_up", self.vbf_maxpt_j1_eta_up_ , "vbf_maxpt_j1_eta_up/F");
        self.otree.Branch("vbf_maxpt_j1_phi_up", self.vbf_maxpt_j1_phi_up_ , "vbf_maxpt_j1_phi_up/F");

        self.otree.Branch("vbf_maxpt_j1_m_dn", self.vbf_maxpt_j1_m_dn_ , "vbf_maxpt_j1_m_dn/F");
        self.otree.Branch("vbf_maxpt_j1_pt_dn", self.vbf_maxpt_j1_pt_dn_ , "vbf_maxpt_j1_pt_dn/F");
        self.otree.Branch("vbf_maxpt_j1_eta_dn", self.vbf_maxpt_j1_eta_dn_ , "vbf_maxpt_j1_eta_dn/F");
        self.otree.Branch("vbf_maxpt_j1_phi_dn", self.vbf_maxpt_j1_phi_dn_ , "vbf_maxpt_j1_phi_dn/F");

        self.otree.Branch("vbf_maxpt_j2_m_up", self.vbf_maxpt_j2_m_up_ , "vbf_maxpt_j2_m_up/F");
        self.otree.Branch("vbf_maxpt_j2_pt_up", self.vbf_maxpt_j2_pt_up_ , "vbf_maxpt_j2_pt_up/F");
        self.otree.Branch("vbf_maxpt_j2_eta_up", self.vbf_maxpt_j2_eta_up_ , "vbf_maxpt_j2_eta_up/F");
        self.otree.Branch("vbf_maxpt_j2_phi_up", self.vbf_maxpt_j2_phi_up_ , "vbf_maxpt_j2_phi_up/F");

        self.otree.Branch("vbf_maxpt_j2_m_dn", self.vbf_maxpt_j2_m_dn_ , "vbf_maxpt_j2_m_dn/F");
        self.otree.Branch("vbf_maxpt_j2_pt_dn", self.vbf_maxpt_j2_pt_dn_ , "vbf_maxpt_j2_pt_dn/F");
        self.otree.Branch("vbf_maxpt_j2_eta_dn", self.vbf_maxpt_j2_eta_dn_ , "vbf_maxpt_j2_eta_dn/F");
        self.otree.Branch("vbf_maxpt_j2_phi_dn", self.vbf_maxpt_j2_phi_dn_ , "vbf_maxpt_j2_phi_dn/F");

        self.otree.Branch("vbf_maxpt_j1_QGLikelihood", self.vbf_maxpt_j1_QGLikelihood_ , "vbf_maxpt_j1_QGLikelihood/F");
        self.otree.Branch("vbf_maxpt_j2_QGLikelihood", self.vbf_maxpt_j2_QGLikelihood_ , "vbf_maxpt_j2_QGLikelihood/F");

        self.otree.Branch("vbf_maxpt_j1_isPileUpMedium", self.vbf_maxpt_j1_isPileUpMedium_ , "vbf_maxpt_j1_isPileUpMedium/I");
        self.otree.Branch("vbf_maxpt_j2_isPileUpMedium", self.vbf_maxpt_j2_isPileUpMedium_ , "vbf_maxpt_j2_isPileUpMedium/I");

        self.otree.Branch("vbf_maxpt_j1_isPileUpTight", self.vbf_maxpt_j1_isPileUpTight_ , "vbf_maxpt_j1_isPileUpTight/I");
        self.otree.Branch("vbf_maxpt_j2_isPileUpTight", self.vbf_maxpt_j2_isPileUpTight_ , "vbf_maxpt_j2_isPileUpTight/I");

        self.otree.Branch("vbf_maxpt_j1_bDiscriminatorCSV", self.vbf_maxpt_j1_bDiscriminatorCSV_ , "vbf_maxpt_j1_bDiscriminatorCSV/F");
        self.otree.Branch("vbf_maxpt_j2_bDiscriminatorCSV", self.vbf_maxpt_j2_bDiscriminatorCSV_ , "vbf_maxpt_j2_bDiscriminatorCSV/F");

        self.otree.Branch("vbf_maxpt_j1_bDiscriminatorCSV_gen", self.vbf_maxpt_j1_bDiscriminatorCSV_gen_ , "vbf_maxpt_j1_bDiscriminatorCSV_gen/F");
        self.otree.Branch("vbf_maxpt_j2_bDiscriminatorCSV_gen", self.vbf_maxpt_j2_bDiscriminatorCSV_gen_ , "vbf_maxpt_j2_bDiscriminatorCSV_gen/F");

        self.otree.Branch("vbf_maxDeta_jj_m", self.vbf_maxDeta_jj_m_ , "vbf_maxDeta_jj_m/F");
        self.otree.Branch("vbf_maxDeta_jj_pt", self.vbf_maxDeta_jj_pt_ , "vbf_maxDeta_jj_pt/F");
        self.otree.Branch("vbf_maxDeta_jj_eta", self.vbf_maxDeta_jj_eta_ , "vbf_maxDeta_jj_eta/F");
        self.otree.Branch("vbf_maxDeta_jj_phi", self.vbf_maxDeta_jj_phi_ , "vbf_maxDeta_jj_phi/F");

        self.otree.Branch("vbf_maxDeta_j1_m", self.vbf_maxDeta_j1_m_ , "vbf_maxDeta_j1_m/F");
        self.otree.Branch("vbf_maxDeta_j1_pt", self.vbf_maxDeta_j1_pt_ , "vbf_maxDeta_j1_pt/F");
        self.otree.Branch("vbf_maxDeta_j1_eta", self.vbf_maxDeta_j1_eta_ , "vbf_maxDeta_j1_eta/F");
        self.otree.Branch("vbf_maxDeta_j1_phi", self.vbf_maxDeta_j1_phi_ , "vbf_maxDeta_j1_phi/F");

        self.otree.Branch("vbf_maxDeta_j2_m", self.vbf_maxDeta_j2_m_ , "vbf_maxDeta_j2_m/F");
        self.otree.Branch("vbf_maxDeta_j2_pt", self.vbf_maxDeta_j2_pt_ , "vbf_maxDeta_j2_pt/F");
        self.otree.Branch("vbf_maxDeta_j2_eta", self.vbf_maxDeta_j2_eta_ , "vbf_maxDeta_j2_eta/F");
        self.otree.Branch("vbf_maxDeta_j2_phi", self.vbf_maxDeta_j2_phi_ , "vbf_maxDeta_j2_phi/F");


        self.otree.Branch("vbf_maxDeta_jj_m_gen", self.vbf_maxDeta_jj_m_gen_ , "vbf_maxDeta_jj_m_gen/F");
        self.otree.Branch("vbf_maxDeta_jj_pt_gen", self.vbf_maxDeta_jj_pt_gen_ , "vbf_maxDeta_jj_pt_gen/F");
        self.otree.Branch("vbf_maxDeta_jj_eta_gen", self.vbf_maxDeta_jj_eta_gen_ , "vbf_maxDeta_jj_eta_gen/F");
        self.otree.Branch("vbf_maxDeta_jj_phi_gen", self.vbf_maxDeta_jj_phi_gen_ , "vbf_maxDeta_jj_phi_gen/F");

        self.otree.Branch("vbf_maxDeta_j1_m_gen", self.vbf_maxDeta_j1_m_gen_ , "vbf_maxDeta_j1_m_gen/F");
        self.otree.Branch("vbf_maxDeta_j1_pt_gen", self.vbf_maxDeta_j1_pt_gen_ , "vbf_maxDeta_j1_pt_gen/F");
        self.otree.Branch("vbf_maxDeta_j1_eta_gen", self.vbf_maxDeta_j1_eta_gen_ , "vbf_maxDeta_j1_eta_gen/F");
        self.otree.Branch("vbf_maxDeta_j1_phi_gen", self.vbf_maxDeta_j1_phi_gen_ , "vbf_maxDeta_j1_phi_gen/F");

        self.otree.Branch("vbf_maxDeta_j2_m_gen", self.vbf_maxDeta_j2_m_gen_ , "vbf_maxDeta_j2_m_gen/F");
        self.otree.Branch("vbf_maxDeta_j2_pt_gen", self.vbf_maxDeta_j2_pt_gen_ , "vbf_maxDeta_j2_pt_gen/F");
        self.otree.Branch("vbf_maxDeta_j2_eta_gen", self.vbf_maxDeta_j2_eta_gen_ , "vbf_maxDeta_j2_eta_gen/F");
        self.otree.Branch("vbf_maxDeta_j2_phi_gen", self.vbf_maxDeta_j2_phi_gen_ , "vbf_maxDeta_j2_phi_gen/F");

        self.otree.Branch("vbf_maxDeta_j1_m_up", self.vbf_maxDeta_j1_m_up_ , "vbf_maxDeta_j1_m_up/F");
        self.otree.Branch("vbf_maxDeta_j1_pt_up", self.vbf_maxDeta_j1_pt_up_ , "vbf_maxDeta_j1_pt_up/F");
        self.otree.Branch("vbf_maxDeta_j1_eta_up", self.vbf_maxDeta_j1_eta_up_ , "vbf_maxDeta_j1_eta_up/F");
        self.otree.Branch("vbf_maxDeta_j1_phi_up", self.vbf_maxDeta_j1_phi_up_ , "vbf_maxDeta_j1_phi_up/F");

        self.otree.Branch("vbf_maxDeta_j1_m_dn", self.vbf_maxDeta_j1_m_dn_ , "vbf_maxDeta_j1_m_dn/F");
        self.otree.Branch("vbf_maxDeta_j1_pt_dn", self.vbf_maxDeta_j1_pt_dn_ , "vbf_maxDeta_j1_pt_dn/F");
        self.otree.Branch("vbf_maxDeta_j1_eta_dn", self.vbf_maxDeta_j1_eta_dn_ , "vbf_maxDeta_j1_eta_dn/F");
        self.otree.Branch("vbf_maxDeta_j1_phi_dn", self.vbf_maxDeta_j1_phi_dn_ , "vbf_maxDeta_j1_phi_dn/F");

        self.otree.Branch("vbf_maxDeta_j2_m_up", self.vbf_maxDeta_j2_m_up_ , "vbf_maxDeta_j2_m_up/F");
        self.otree.Branch("vbf_maxDeta_j2_pt_up", self.vbf_maxDeta_j2_pt_up_ , "vbf_maxDeta_j2_pt_up/F");
        self.otree.Branch("vbf_maxDeta_j2_eta_up", self.vbf_maxDeta_j2_eta_up_ , "vbf_maxDeta_j2_eta_up/F");
        self.otree.Branch("vbf_maxDeta_j2_phi_up", self.vbf_maxDeta_j2_phi_up_ , "vbf_maxDeta_j2_phi_up/F");

        self.otree.Branch("vbf_maxDeta_j2_m_dn", self.vbf_maxDeta_j2_m_dn_ , "vbf_maxDeta_j2_m_dn/F");
        self.otree.Branch("vbf_maxDeta_j2_pt_dn", self.vbf_maxDeta_j2_pt_dn_ , "vbf_maxDeta_j2_pt_dn/F");
        self.otree.Branch("vbf_maxDeta_j2_eta_dn", self.vbf_maxDeta_j2_eta_dn_ , "vbf_maxDeta_j2_eta_dn/F");
        self.otree.Branch("vbf_maxDeta_j2_phi_dn", self.vbf_maxDeta_j2_phi_dn_ , "vbf_maxDeta_j2_phi_dn/F");

        self.otree.Branch("vbf_maxDeta_j1_QGLikelihood", self.vbf_maxDeta_j1_QGLikelihood_ , "vbf_maxDeta_j1_QGLikelihood/F");
        self.otree.Branch("vbf_maxDeta_j2_QGLikelihood", self.vbf_maxDeta_j2_QGLikelihood_ , "vbf_maxDeta_j2_QGLikelihood/F");

        self.otree.Branch("vbf_maxDeta_j1_isPileUpMedium", self.vbf_maxDeta_j1_isPileUpMedium_ , "vbf_maxDeta_j1_isPileUpMedium/I");
        self.otree.Branch("vbf_maxDeta_j2_isPileUpMedium", self.vbf_maxDeta_j2_isPileUpMedium_ , "vbf_maxDeta_j2_isPileUpMedium/I");

        self.otree.Branch("vbf_maxDeta_j1_isPileUpTight", self.vbf_maxDeta_j1_isPileUpTight_ , "vbf_maxDeta_j1_isPileUpTight/I");
        self.otree.Branch("vbf_maxDeta_j2_isPileUpTight", self.vbf_maxDeta_j2_isPileUpTight_ , "vbf_maxDeta_j2_isPileUpTight/I");

        self.otree.Branch("vbf_maxDeta_j1_bDiscriminatorCSV", self.vbf_maxDeta_j1_bDiscriminatorCSV_ , "vbf_maxDeta_j1_bDiscriminatorCSV/F");
        self.otree.Branch("vbf_maxDeta_j2_bDiscriminatorCSV", self.vbf_maxDeta_j2_bDiscriminatorCSV_ , "vbf_maxDeta_j2_bDiscriminatorCSV/F");

        self.otree.Branch("vbf_maxDeta_j1_bDiscriminatorCSV_gen", self.vbf_maxDeta_j1_bDiscriminatorCSV_gen_ , "vbf_maxDeta_j1_bDiscriminatorCSV_gen/F");
        self.otree.Branch("vbf_maxDeta_j2_bDiscriminatorCSV_gen", self.vbf_maxDeta_j2_bDiscriminatorCSV_gen_ , "vbf_maxDeta_j2_bDiscriminatorCSV_gen/F");

        self.otree.Branch("vbf_maxMjj_jj_m", self.vbf_maxMjj_jj_m_ , "vbf_maxMjj_jj_m/F");
        self.otree.Branch("vbf_maxMjj_jj_pt", self.vbf_maxMjj_jj_pt_ , "vbf_maxMjj_jj_pt/F");
        self.otree.Branch("vbf_maxMjj_jj_eta", self.vbf_maxMjj_jj_eta_ , "vbf_maxMjj_jj_eta/F");
        self.otree.Branch("vbf_maxMjj_jj_phi", self.vbf_maxMjj_jj_phi_ , "vbf_maxMjj_jj_phi/F");

        self.otree.Branch("vbf_maxMjj_j1_m", self.vbf_maxMjj_j1_m_ , "vbf_maxMjj_j1_m/F");
        self.otree.Branch("vbf_maxMjj_j1_pt", self.vbf_maxMjj_j1_pt_ , "vbf_maxMjj_j1_pt/F");
        self.otree.Branch("vbf_maxMjj_j1_eta", self.vbf_maxMjj_j1_eta_ , "vbf_maxMjj_j1_eta/F");
        self.otree.Branch("vbf_maxMjj_j1_phi", self.vbf_maxMjj_j1_phi_ , "vbf_maxMjj_j1_phi/F");

        self.otree.Branch("vbf_maxMjj_j2_m", self.vbf_maxMjj_j2_m_ , "vbf_maxMjj_j2_m/F");
        self.otree.Branch("vbf_maxMjj_j2_pt", self.vbf_maxMjj_j2_pt_ , "vbf_maxMjj_j2_pt/F");
        self.otree.Branch("vbf_maxMjj_j2_eta", self.vbf_maxMjj_j2_eta_ , "vbf_maxMjj_j2_eta/F");
        self.otree.Branch("vbf_maxMjj_j2_phi", self.vbf_maxMjj_j2_phi_ , "vbf_maxMjj_j2_phi/F");

        self.otree.Branch("vbf_maxMjj_jj_m_gen", self.vbf_maxMjj_jj_m_gen_ , "vbf_maxMjj_jj_m_gen/F");
        self.otree.Branch("vbf_maxMjj_jj_pt_gen", self.vbf_maxMjj_jj_pt_gen_ , "vbf_maxMjj_jj_pt_gen/F");
        self.otree.Branch("vbf_maxMjj_jj_eta_gen", self.vbf_maxMjj_jj_eta_gen_ , "vbf_maxMjj_jj_eta_gen/F");
        self.otree.Branch("vbf_maxMjj_jj_phi_gen", self.vbf_maxMjj_jj_phi_gen_ , "vbf_maxMjj_jj_phi_gen/F");

        self.otree.Branch("vbf_maxMjj_j1_m_gen", self.vbf_maxMjj_j1_m_gen_ , "vbf_maxMjj_j1_m_gen/F");
        self.otree.Branch("vbf_maxMjj_j1_pt_gen", self.vbf_maxMjj_j1_pt_gen_ , "vbf_maxMjj_j1_pt_gen/F");
        self.otree.Branch("vbf_maxMjj_j1_eta_gen", self.vbf_maxMjj_j1_eta_gen_ , "vbf_maxMjj_j1_eta_gen/F");
        self.otree.Branch("vbf_maxMjj_j1_phi_gen", self.vbf_maxMjj_j1_phi_gen_ , "vbf_maxMjj_j1_phi_gen/F");

        self.otree.Branch("vbf_maxMjj_j2_m_gen", self.vbf_maxMjj_j2_m_gen_ , "vbf_maxMjj_j2_m_gen/F");
        self.otree.Branch("vbf_maxMjj_j2_pt_gen", self.vbf_maxMjj_j2_pt_gen_ , "vbf_maxMjj_j2_pt_gen/F");
        self.otree.Branch("vbf_maxMjj_j2_eta_gen", self.vbf_maxMjj_j2_eta_gen_ , "vbf_maxMjj_j2_eta_gen/F");
        self.otree.Branch("vbf_maxMjj_j2_phi_gen", self.vbf_maxMjj_j2_phi_gen_ , "vbf_maxMjj_j2_phi_gen/F");

        self.otree.Branch("vbf_maxMjj_j1_m_up", self.vbf_maxMjj_j1_m_up_ , "vbf_maxMjj_j1_m_up/F");
        self.otree.Branch("vbf_maxMjj_j1_pt_up", self.vbf_maxMjj_j1_pt_up_ , "vbf_maxMjj_j1_pt_up/F");
        self.otree.Branch("vbf_maxMjj_j1_eta_up", self.vbf_maxMjj_j1_eta_up_ , "vbf_maxMjj_j1_eta_up/F");
        self.otree.Branch("vbf_maxMjj_j1_phi_up", self.vbf_maxMjj_j1_phi_up_ , "vbf_maxMjj_j1_phi_up/F");

        self.otree.Branch("vbf_maxMjj_j1_m_dn", self.vbf_maxMjj_j1_m_dn_ , "vbf_maxMjj_j1_m_dn/F");
        self.otree.Branch("vbf_maxMjj_j1_pt_dn", self.vbf_maxMjj_j1_pt_dn_ , "vbf_maxMjj_j1_pt_dn/F");
        self.otree.Branch("vbf_maxMjj_j1_eta_dn", self.vbf_maxMjj_j1_eta_dn_ , "vbf_maxMjj_j1_eta_dn/F");
        self.otree.Branch("vbf_maxMjj_j1_phi_dn", self.vbf_maxMjj_j1_phi_dn_ , "vbf_maxMjj_j1_phi_dn/F");

        self.otree.Branch("vbf_maxMjj_j2_m_up", self.vbf_maxMjj_j2_m_up_ , "vbf_maxMjj_j2_m_up/F");
        self.otree.Branch("vbf_maxMjj_j2_pt_up", self.vbf_maxMjj_j2_pt_up_ , "vbf_maxMjj_j2_pt_up/F");
        self.otree.Branch("vbf_maxMjj_j2_eta_up", self.vbf_maxMjj_j2_eta_up_ , "vbf_maxMjj_j2_eta_up/F");
        self.otree.Branch("vbf_maxMjj_j2_phi_up", self.vbf_maxMjj_j2_phi_up_ , "vbf_maxMjj_j2_phi_up/F");

        self.otree.Branch("vbf_maxMjj_j2_m_dn", self.vbf_maxMjj_j2_m_dn_ , "vbf_maxMjj_j2_m_dn/F");
        self.otree.Branch("vbf_maxMjj_j2_pt_dn", self.vbf_maxMjj_j2_pt_dn_ , "vbf_maxMjj_j2_pt_dn/F");
        self.otree.Branch("vbf_maxMjj_j2_eta_dn", self.vbf_maxMjj_j2_eta_dn_ , "vbf_maxMjj_j2_eta_dn/F");
        self.otree.Branch("vbf_maxMjj_j2_phi_dn", self.vbf_maxMjj_j2_phi_dn_ , "vbf_maxMjj_j2_phi_dn/F");

        self.otree.Branch("vbf_maxMjj_j1_QGLikelihood", self.vbf_maxMjj_j1_QGLikelihood_ , "vbf_maxMjj_j1_QGLikelihood/F");
        self.otree.Branch("vbf_maxMjj_j2_QGLikelihood", self.vbf_maxMjj_j2_QGLikelihood_ , "vbf_maxMjj_j2_QGLikelihood/F");

        self.otree.Branch("vbf_maxMjj_j1_isPileUpMedium", self.vbf_maxMjj_j1_isPileUpMedium_ , "vbf_maxMjj_j1_isPileUpMedium/I");
        self.otree.Branch("vbf_maxMjj_j2_isPileUpMedium", self.vbf_maxMjj_j2_isPileUpMedium_ , "vbf_maxMjj_j2_isPileUpMedium/I");

        self.otree.Branch("vbf_maxMjj_j1_isPileUpTight", self.vbf_maxMjj_j1_isPileUpTight_ , "vbf_maxMjj_j1_isPileUpTight/I");
        self.otree.Branch("vbf_maxMjj_j2_isPileUpTight", self.vbf_maxMjj_j2_isPileUpTight_ , "vbf_maxMjj_j2_isPileUpTight/I");

        self.otree.Branch("vbf_maxMjj_j1_bDiscriminatorCSV", self.vbf_maxMjj_j1_bDiscriminatorCSV_ , "vbf_maxMjj_j1_bDiscriminatorCSV/F");
        self.otree.Branch("vbf_maxMjj_j2_bDiscriminatorCSV", self.vbf_maxMjj_j2_bDiscriminatorCSV_ , "vbf_maxMjj_j2_bDiscriminatorCSV/F");

        self.otree.Branch("vbf_maxMjj_j1_bDiscriminatorCSV_gen", self.vbf_maxMjj_j1_bDiscriminatorCSV_gen_ , "vbf_maxMjj_j1_bDiscriminatorCSV_gen/F");
        self.otree.Branch("vbf_maxMjj_j2_bDiscriminatorCSV_gen", self.vbf_maxMjj_j2_bDiscriminatorCSV_gen_ , "vbf_maxMjj_j2_bDiscriminatorCSV_gen/F");

        ###### btag counters

        self.otree.Branch("nbjets_csvl_veto", self.nbjets_csvl_veto_ , "nbjets_csvl_veto/F");
        self.otree.Branch("nbjets_csvm_veto", self.nbjets_csvm_veto_ , "nbjets_csvm_veto/F");
        self.otree.Branch("nbjets_csvt_veto", self.nbjets_csvt_veto_ , "nbjets_csvt_veto/F");
        self.otree.Branch("nbjets_ssvhem_veto", self.nbjets_ssvhem_veto_ , "nbjets_ssvhem_veto/F");

        self.otree.Branch("nbjets_csvl_veto_cleaned", self.nbjets_csvl_veto_cleaned_ , "nbjets_csvl_veto_cleaned/F");
        self.otree.Branch("nbjets_csvm_veto_cleaned", self.nbjets_csvm_veto_cleaned_ , "nbjets_csvm_veto_cleaned/F");
        self.otree.Branch("nbjets_csvt_veto_cleaned", self.nbjets_csvt_veto_cleaned_ , "nbjets_csvt_veto_cleaned/F");
        self.otree.Branch("nbjets_ssvhem_veto_cleaned", self.nbjets_ssvhem_veto_cleaned_ , "nbjets_ssvhem_veto_cleaned/F");

        self.otree.Branch("njets", self.njets_ , "njets/F");

        self.otree.Branch("deltaR_lca8jet",     self.deltaR_lca8jet_ , "deltaR_lca8jet/F")
        self.otree.Branch("deltaphi_METca8jet", self.deltaphi_METca8jet_ , "deltaphi_METca8jet/F")
        self.otree.Branch("deltaphi_Vca8jet",   self.deltaphi_Vca8jet_ , "deltaphi_Vca8jet/F")
        self.otree.Branch("deltaphi_METca8jet_met", self.deltaphi_METca8jet_met_ , "deltaphi_METca8jet_met/F")
        self.otree.Branch("deltaphi_Vca8jet_met", self.deltaphi_Vca8jet_met_ , "deltaphi_Vca8jet_met/F")

        ### some generator level quantites
        self.otree.Branch("genHMass",     self.genHMass_ , "genHMass/F")
        self.otree.Branch("genHphi",     self.genHphi_ , "genHphi/F")
        self.otree.Branch("genHeta",     self.genHeta_ , "genHeta/F")
        self.otree.Branch("genHpt",     self.genHpt_ , "genHpt/F")
                
        self.otree.Branch("genTagQuark1W",       self.genTagQuark1E_ , "genTagQuark1E/F")
        self.otree.Branch("genTagQuark1phi",     self.genTagQuark1phi_ , "genTagQuark1phi/F")
        self.otree.Branch("genTagQuark1eta",     self.genTagQuark1eta_ , "genTagQuark1eta/F")
        self.otree.Branch("genTagQuark1pt",      self.genTagQuark1pt_ , "genTagQuark1pt/F")

        self.otree.Branch("genTagQuarkE",     self.genTagQuark2E_ , "genTagQuark2E/F")
        self.otree.Branch("genTagQuark2phi",     self.genTagQuark2phi_ , "genTagQuark2phi/F")
        self.otree.Branch("genTagQuark2eta",     self.genTagQuark2eta_ , "genTagQuark2eta/F")
        self.otree.Branch("genTagQuark2pt",     self.genTagQuark2pt_ , "genTagQuark2pt/F")
                                                                        
        # variables for ttbar control region
        self.otree.Branch("ttb_nak5_same",     self.ttb_nak5_same_ , "ttb_nak5_same/F")
        self.otree.Branch("ttb_nak5_same_csvl",     self.ttb_nak5_same_csvl_ , "ttb_nak5_same_csvl/F")
        self.otree.Branch("ttb_nak5_same_csvm",     self.ttb_nak5_same_csvm_ , "ttb_nak5_same_csvm/F")
        self.otree.Branch("ttb_nak5_same_csvt",     self.ttb_nak5_same_csvt_ , "ttb_nak5_same_csvt/F")

        self.otree.Branch("ttb_nak5_oppo",     self.ttb_nak5_oppo_ , "ttb_nak5_oppo/F")
        self.otree.Branch("ttb_nak5_oppo_csvl",     self.ttb_nak5_oppo_csvl_ , "ttb_nak5_oppo_csvl/F")
        self.otree.Branch("ttb_nak5_oppo_csvm",     self.ttb_nak5_oppo_csvm_ , "ttb_nak5_oppo_csvm/F")
        self.otree.Branch("ttb_nak5_oppo_csvt",     self.ttb_nak5_oppo_csvt_ , "ttb_nak5_oppo_csvt/F")

        self.otree.Branch("ttb_nak5_oppoveto",     self.ttb_nak5_oppoveto_ , "ttb_nak5_oppoveto/F")
        self.otree.Branch("ttb_nak5_oppoveto_csvl",     self.ttb_nak5_oppoveto_csvl_ , "ttb_nak5_oppoveto_csvl/F")
        self.otree.Branch("ttb_nak5_oppoveto_csvm",     self.ttb_nak5_oppoveto_csvm_ , "ttb_nak5_oppoveto_csvm/F")
        self.otree.Branch("ttb_nak5_oppoveto_csvt",     self.ttb_nak5_oppoveto_csvt_ , "ttb_nak5_oppoveto_csvt/F")

        self.otree.Branch("ttb_nak5_sameveto",     self.ttb_nak5_sameveto_ , "ttb_nak5_sameveto/F")
        self.otree.Branch("ttb_nak5_sameveto_csvl",     self.ttb_nak5_sameveto_csvl_ , "ttb_nak5_sameveto_csvl/F")
        self.otree.Branch("ttb_nak5_sameveto_csvm",     self.ttb_nak5_sameveto_csvm_ , "ttb_nak5_sameveto_csvm/F")
        self.otree.Branch("ttb_nak5_sameveto_csvt",     self.ttb_nak5_sameveto_csvt_ , "ttb_nak5_sameveto_csvt/F")

        self.otree.Branch("ttb_ht",     self.ttb_ht_ , "ttb_ht/F")
        self.otree.Branch("ttb_ca8_mass_pr",     self.ttb_ca8_mass_pr_ , "ttb_ca8_mass_pr/F")
        self.otree.Branch("ttb_ca8_charge",     self.ttb_ca8_charge_ , "ttb_ca8_charge/F")
        self.otree.Branch("ttb_ca8_charge_k05",     self.ttb_ca8_charge_k05_ , "ttb_ca8_charge_k05/F")
        self.otree.Branch("ttb_ca8_charge_k07",     self.ttb_ca8_charge_k07_ , "ttb_ca8_charge_k07/F")
        self.otree.Branch("ttb_ca8_charge_k10",     self.ttb_ca8_charge_k10_ , "ttb_ca8_charge_k10/F")

        self.otree.Branch("ttb_ca8_ungroomed_pt",     self.ttb_ca8_ungroomed_pt_ , "ttb_ca8_ungroomed_pt/F")
        self.otree.Branch("ttb_ca8_ungroomed_eta",     self.ttb_ca8_ungroomed_eta_ , "ttb_ca8_ungroomed_eta/F")
        self.otree.Branch("ttb_ca8_ungroomed_phi",     self.ttb_ca8_ungroomed_phi_ , "ttb_ca8_ungroomed_phi/F")
        self.otree.Branch("ttb_ca8_ungroomed_e",     self.ttb_ca8_ungroomed_e_ , "ttb_ca8_ungroomed_e/F")
        

        self.otree.Branch("ttb_ca8_ungroomed_gen_pt",     self.ttb_ca8_ungroomed_gen_pt_ , "ttb_ca8_ungroomed_gen_pt/F")
        self.otree.Branch("ttb_ca8_ungroomed_gen_eta",     self.ttb_ca8_ungroomed_gen_eta_ , "ttb_ca8_ungroomed_gen_eta/F")
        self.otree.Branch("ttb_ca8_ungroomed_gen_phi",     self.ttb_ca8_ungroomed_gen_phi_ , "ttb_ca8_ungroomed_gen_phi/F")
        self.otree.Branch("ttb_ca8_ungroomed_gen_e",     self.ttb_ca8_ungroomed_gen_e_ , "ttb_ca8_ungroomed_gen_e/F")
        
        self.otree.Branch("ttb_ca8_tau2tau1",     self.ttb_ca8_tau2tau1_ , "ttb_ca8_tau2tau1/F")
        self.otree.Branch("ttb_ca8_tau2tau1_exkT",     self.ttb_ca8_tau2tau1_exkT_ , "ttb_ca8_tau2tau1_exkT/F")
        self.otree.Branch("ttb_ca8_tau2tau1_pr",     self.ttb_ca8_tau2tau1_pr_ , "ttb_ca8_tau2tau1_pr/F")
        self.otree.Branch("ttb_ca8_GeneralizedECF",     self.ttb_ca8_GeneralizedECF_ , "ttb_ca8_GeneralizedECF/F")
        self.otree.Branch("ttb_ca8_mu",     self.ttb_ca8_mu_ , "ttb_ca8_mu/F")

        self.otree.Branch("ttb_ca8_mlvj_type0",     self.ttb_ca8_mlvj_type0_ , "ttb_ca8_mlvj_type0/F")
        self.otree.Branch("ttb_ca8_mlvj_type2",     self.ttb_ca8_mlvj_type2_ , "ttb_ca8_mlvj_type2/F")

        self.otree.Branch("ttb_ca8_mlvj_type0_met",     self.ttb_ca8_mlvj_type0_met_ , "ttb_ca8_mlvj_type0_met/F")
        self.otree.Branch("ttb_ca8_mlvj_type2_met",     self.ttb_ca8_mlvj_type2_met_ , "ttb_ca8_mlvj_type2_met/F")
    
        self.otree.Branch("isttbar",     self.isttbar_ , "isttbar/I")

        self.otree.Branch("gen_parton1_px_fromttbar",     self.gen_parton1_px_fromttbar_ , "gen_parton1_px_fromttbar/F")
        self.otree.Branch("gen_parton1_py_fromttbar",     self.gen_parton1_py_fromttbar_ , "gen_parton1_py_fromttbar/F")
        self.otree.Branch("gen_parton1_pz_fromttbar",     self.gen_parton1_pz_fromttbar_ , "gen_parton1_pz_fromttbar/F")
        self.otree.Branch("gen_parton1_e_fromttbar",     self.gen_parton1_e_fromttbar_ , "gen_parton1_e_fromttbar/F")
        self.otree.Branch("gen_parton1_id_fromttbar",     self.gen_parton1_id_fromttbar_ , "gen_parton1_id_fromttbar/F")

        self.otree.Branch("gen_parton2_px_fromttbar",     self.gen_parton2_px_fromttbar_ , "gen_parton2_px_fromttbar/F")
        self.otree.Branch("gen_parton2_py_fromttbar",     self.gen_parton2_py_fromttbar_ , "gen_parton2_py_fromttbar/F")
        self.otree.Branch("gen_parton2_pz_fromttbar",     self.gen_parton2_pz_fromttbar_ , "gen_parton2_pz_fromttbar/F")
        self.otree.Branch("gen_parton2_e_fromttbar",     self.gen_parton2_e_fromttbar_ , "gen_parton2_e_fromttbar/F")
        self.otree.Branch("gen_parton2_id_fromttbar",     self.gen_parton2_id_fromttbar_ , "gen_parton2_id_fromttbar/F")

        ### branches for bsm models 
        for i in range(len(self.cprimeVals)): 
            for j in range(len(self.brnewVals)): 
                brname = "bsmReweight_cPrime%02d_brNew%02d"%(self.cprimeVals[i],self.brnewVals[j]);
                self.otree.Branch(brname, self.bsmReweights[i][j] , brname);        

