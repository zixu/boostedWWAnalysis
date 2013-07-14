
# A class which takes histograms and plots them in a versatile way
# inputs are file names which can be "data" or "MC"

import ROOT
ROOT.gROOT.ProcessLine(".L ./tdrstyle.C")
from ROOT import setTDRStyle
from ROOT import TTree
from ROOT import gROOT
ROOT.setTDRStyle()
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(0)

ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);

import os
import sys

from array import array
import math
from optparse import OptionParser

from bsmReweighter import *

from ROOT import RooTrace

# FWLITE stuff
import sys
from DataFormats.FWLite import Events, Handle
ROOT.gSystem.Load('libCondFormatsJetMETObjects')
ROOT.gSystem.Load('libFWCoreFWLite');
ROOT.gSystem.Load('libFWCoreUtilities');  

### ------------ h e l p e r s --------------------

def getListRMS(list):
    mean = sum(list)/float(len(list));
    return math.sqrt(sum((n-mean)*(n-mean) for n in list)/len(list));
def getListMean(list):
    return sum(list)/float(len(list));

class sampleWrapperClass:

    ### ------------------------------------------------
    def __init__(self, label, file, channel, sampleEffLumi, lumi, treename, isData,outputfiledirectory):

        
        self.IsData_ = isData;
        self.FileName_ = file;
        self.File_ = ROOT.TFile(file);
        self.InputTree_ = self.File_.Get(treename);
        
        self.SampleWeight_ = lumi/sampleEffLumi;
        
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
            #print h_rwCPS.GetName();
            #prin"Integral of rWCPS histogram: ", h_rwCPS.Integral();
            self.x_rwCPS = self.h_rwCPS.GetXaxis();
        ######        
        
        
        # ---------- Set up jet corrections on the fly of R >= 0.7 jets
        fDir = "JECs/"      
        jecUncStr = ROOT.std.string(fDir + "GR_R_53_V10_Uncertainty_AK7PFchs.txt")
        self.jecUnc_ = ROOT.JetCorrectionUncertainty(jecUncStr)
    
    def getLabel(self):
        return self.Label_
            
    def createTrainingTree(self):
        
        print self.FileName_
        self.File_ = ROOT.TFile(self.FileName_);
        self.InputTree_ = self.File_.Get("WJet");
        self.NTree_ = self.InputTree_.GetEntries();
        
        print "Turning off branches...", self.FileName_

        # turn off unnecessary branches
        self.turnOffBranches();
        
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

        # read in original file and get lineshape 
        if not self.FitSMSignal: 
            fitSM = FitMassPoint(self.FileName_, massin, massmin[massin], massmax[massin]);
            self.FitSMSignal_mean = fitSM[0];
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

        otree = ROOT.TTree("otree","otree");        
        
        mass_lvj_         = array( 'f', [ 0. ] );
        v_pt_             = array( 'f', [ 0. ] );
        v_mt_             = array( 'f', [ 0. ] );

        jet_mass_pr_       = array( 'f', [ 0. ] );
        jet_pt_pr_         = array( 'f', [ 0. ] );
        ungroomed_jet_eta_ = array( 'f', [ 0. ] );
        ungroomed_jet_pt_  = array( 'f', [ 0. ] );
        ungroomed_jet_phi_ = array( 'f', [ 0. ] );

        jet_mass_pr_exo_      = array( 'f', [ 0. ] );
        jet_pt_pr_exo_        = array( 'f', [ 0. ] );
        ungroomed_jet_pt_exo_ = array( 'f', [ 0. ] );

        j_jecfactor_up_ = array( 'f', [ 0. ] );
        j_jecfactor_dn_ = array( 'f', [ 0. ] );
        jet_mass_pr_up_ = array( 'f', [ 0. ] );
        jet_mass_pr_dn_ = array( 'f', [ 0. ] );

        issignal_     = array( 'i', [ 0 ] );             
        issignal_exo_ = array( 'i', [ 0 ] );             
        
        l_pt_ = array( 'f', [ 0. ] );
        l_eta_ = array( 'f', [ 0. ] );
        l_phi_ = array( 'f', [ 0. ] );
        
        mvaMET_    = array( 'f', [ 0. ] );                
        pfMET_     = array( 'f', [ 0. ] );
        pfMET_Phi_ = array( 'f', [ 0. ] );

        nu_pz_type0_ = array( 'f', [ 0. ] );
        nu_pz_type1_ = array( 'f', [ 0. ] );
        nu_pz_type2_ = array( 'f', [ 0. ] );
        nu_pz_type3_ = array( 'f', [ 0. ] );

        W_pz_type0_ = array( 'f', [ 0. ] );
        W_pz_type1_ = array( 'f', [ 0. ] );
        W_pz_type2_ = array( 'f', [ 0. ] );
        W_pz_type3_ = array( 'f', [ 0. ] );

        nu_pz_gen_   = array( 'f', [ 0. ] );
        W_pz_gen_   = array( 'f', [ 0. ] );

        mass_lvj_type0_   = array( 'f', [ 0. ] );
        mass_lvj_type1_   = array( 'f', [ 0. ] );
        mass_lvj_type2_   = array( 'f', [ 0. ] );
        mass_lvj_type3_   = array( 'f', [ 0. ] );

        v_mass_type0_    = array( 'f', [ 0. ] );
        v_mass_type1_    = array( 'f', [ 0. ] );
        v_mass_type2_    = array( 'f', [ 0. ] );
        v_mass_type3_    = array( 'f', [ 0. ] );
                
        nPV_ = array( 'f', [ 0. ] );                        

        totalEventWeight_ = array( 'f', [ 0. ] );                                
        eff_and_pu_Weight_ = array( 'f', [ 0. ] );                                
        wSampleWeight_ = array( 'f', [ 0. ] ); #wSampleWeight*effwt*puwt                               

        interference_Weight_H600_ = array( 'f', [ 0. ] );                                
        interference_Weight_H700_ = array( 'f', [ 0. ] );                                
        interference_Weight_H800_ = array( 'f', [ 0. ] );                                
        interference_Weight_H900_ = array( 'f', [ 0. ] );                                
        interference_Weight_H1000_ = array( 'f', [ 0. ] );                                

        cps_Weight_H600_ = array( 'f', [ 0. ] );                                
        cps_Weight_H700_ = array( 'f', [ 0. ] );                                
        cps_Weight_H800_ = array( 'f', [ 0. ] );                                
        cps_Weight_H900_ = array( 'f', [ 0. ] );                                
        cps_Weight_H1000_ = array( 'f', [ 0. ] );                                
        
        jet_grsens_ft_ = array( 'f', [ 0. ] );
        jet_grsens_tr_ = array( 'f', [ 0. ] );
        jet_massdrop_pr_ = array( 'f', [ 0. ] );    
        jet_qjetvol_ = array( 'f', [ 0. ] ); 
        jet_tau2tau1_ = array( 'f', [ 0. ] );     
        jet_jetconstituents_ = array( 'f', [ 0. ] );     
        
        jet_rcore4_ = array( 'f', [ 0. ] );
        jet_rcore5_ = array( 'f', [ 0. ] );
        jet_rcore6_ = array( 'f', [ 0. ] );
        jet_rcore7_ = array( 'f', [ 0. ] );

        jet_grsens_ft_exo_ = array( 'f', [ 0. ] );
        jet_grsens_tr_exo_ = array( 'f', [ 0. ] );
        jet_massdrop_pr_exo_ = array( 'f', [ 0. ] );    
        jet_tau2tau1_exo_ = array( 'f', [ 0. ] );     
        jet_jetconstituents_exo_ = array( 'f', [ 0. ] );     

        jet_pt1frac_  = array ('f',[ 0. ]);
        jet_pt2frac_  = array ('f',[ 0. ]);
        jet_sjdr_     = array ('f',[ 0. ]);

        nbjets_ = array( 'f', [ 0. ] );
        nbjets_cvsl_ = array( 'f', [ 0. ] );
        nbjets_cvsm_ = array( 'f', [ 0. ] );
        nbjets_ssvhem_ = array( 'f', [ 0. ] );

        nbjets_csvl_veto_ = array( 'f', [ 0. ] );
        nbjets_csvm_veto_ = array( 'f', [ 0. ] );
        nbjets_csvt_veto_ = array( 'f', [ 0. ] );
        nbjets_ssvhem_veto_ = array( 'f', [ 0. ] );

        nbjets_csvl_veto_cleaned_ = array( 'f', [ 0. ] );
        nbjets_csvm_veto_cleaned_ = array( 'f', [ 0. ] );
        nbjets_csvt_veto_cleaned_ = array( 'f', [ 0. ] );
        nbjets_ssvhem_veto_cleaned_ = array( 'f', [ 0. ] );

        nbjetsCSV_ = array( 'f', [ 0. ] );
        nbjetsSSVHE_ = array( 'f', [ 0. ] );        
        njets_ = array( 'f', [ 0. ] );


        jet_planarlow04_ = array( 'f', [0.] );
        jet_planarlow05_ = array( 'f', [0.] );
        jet_planarlow06_ = array( 'f', [0.] );
        jet_planarlow07_ = array( 'f', [0.] );
        deltaR_lca8jet_ = array( 'f', [0.] );
        deltaphi_METca8jet_ = array( 'f', [0.] );
        deltaphi_Vca8jet_ = array( 'f', [0.] );
                                                                                
        event_runNo_      = array( 'i', [ 0 ] );
        event_lumi_       = array( 'i', [ 0 ] );
        event_ = array( 'i', [0] );

        otree.Branch("event", event_ , "event/I");
        otree.Branch("event_runNo", event_runNo_ , "event_runNo/I");
        otree.Branch("event_lumi", event_lumi_ , "event_lumi/I");

        otree.Branch("mass_lvj", mass_lvj_ , "mass_lvj/F");
        otree.Branch("v_pt", v_pt_ , "v_pt/F");
        otree.Branch("v_mt", v_mt_ , "v_mt/F");

        otree.Branch("nu_pz_type0", nu_pz_type0_ , "nu_pz_type0/F");
        otree.Branch("nu_pz_type1", nu_pz_type1_ , "nu_pz_type1/F");
        otree.Branch("nu_pz_type2", nu_pz_type2_ , "nu_pz_type2/F");
        otree.Branch("nu_pz_type3", nu_pz_type3_ , "nu_pz_type3/F");

        otree.Branch("W_pz_type0", W_pz_type0_ , "W_pz_type0/F");
        otree.Branch("W_pz_type1", W_pz_type1_ , "W_pz_type1/F");
        otree.Branch("W_pz_type2", W_pz_type2_ , "W_pz_type2/F");
        otree.Branch("W_pz_type3", W_pz_type3_ , "W_pz_type3/F");

        otree.Branch("nu_pz_gen", nu_pz_gen_ , "nu_pz_gen/F");
        otree.Branch("W_pz_gen", W_pz_gen_ , "W_pz_gen/F");

        otree.Branch("mass_lvj_type0", mass_lvj_type0_ , "mass_lvj_type0/F");
        otree.Branch("mass_lvj_type1", mass_lvj_type1_ , "mass_lvj_type1/F");
        otree.Branch("mass_lvj_type2", mass_lvj_type2_ , "mass_lvj_type2/F");
        otree.Branch("mass_lvj_type3", mass_lvj_type3_ , "mass_lvj_type3/F");

        otree.Branch("v_mass_type0", v_mass_type0_ , "v_mass_type0/F");
        otree.Branch("v_mass_type1", v_mass_type1_ , "v_mass_type1/F");
        otree.Branch("v_mass_type2", v_mass_type2_ , "v_mass_type2/F");
        otree.Branch("v_mass_type3", v_mass_type3_ , "v_mass_type3/F");

        otree.Branch("jet_pt_pr", jet_pt_pr_ , "jet_pt_pr/F");
        otree.Branch("ungroomed_jet_pt", ungroomed_jet_pt_, "ungroomed_jet_pt/F");
        otree.Branch("ungroomed_jet_phi", ungroomed_jet_phi_, "ungroomed_jet_phi/F");
        otree.Branch("ungroomed_jet_eta", ungroomed_jet_eta_, "ungroomed_jet_eta/F");

        otree.Branch("jet_mass_pr", jet_mass_pr_ , "jet_mass_pr/F");
        otree.Branch("j_jecfactor_up", j_jecfactor_up_ , "j_jecfactor_up/F");
        otree.Branch("j_jecfactor_dn", j_jecfactor_dn_ , "j_jecfactor_dn/F");        
        otree.Branch("jet_mass_pr_up", jet_mass_pr_up_ , "jet_mass_pr_up/F");
        otree.Branch("jet_mass_pr_dn", jet_mass_pr_dn_ , "jet_mass_pr_dn/F");

        otree.Branch("jet_pt_pr_exo", jet_pt_pr_exo_ , "jet_pt_pr_exo/F");
        otree.Branch("ungroomed_jet_pt_exo", ungroomed_jet_pt_exo_, "ungroomed_jet_pt_exo/F");
        otree.Branch("jet_mass_pr_exo", jet_mass_pr_exo_ , "jet_mass_pr_exo/F");

        otree.Branch("issignal", issignal_ , "issignal/I");        
        otree.Branch("issignal_exo", issignal_exo_ , "issignal_exo/I");        
        
        otree.Branch("l_pt", l_pt_ , "l_pt/F");
        otree.Branch("l_eta", l_eta_ , "l_eta/F");
        otree.Branch("l_phi", l_phi_ , "l_phi/F");

        otree.Branch("mvaMET", mvaMET_ , "mvaMET/F");
        otree.Branch("pfMET", pfMET_ , "pfMET/F");        
        otree.Branch("pfMET_Phi", pfMET_Phi_ , "pfMET_Phi/F");        
        otree.Branch("nPV", nPV_ , "nPV/F");

        otree.Branch("totalEventWeight", totalEventWeight_ , "totalEventWeight/F");
        otree.Branch("eff_and_pu_Weight", eff_and_pu_Weight_ , "eff_and_pu_Weight/F");
        otree.Branch("wSampleWeight", wSampleWeight_ , "wSampleWeight/F");
        otree.Branch("interference_Weight_H600", interference_Weight_H600_ , "interference_Weight_H600/F");
        otree.Branch("interference_Weight_H700", interference_Weight_H700_ , "interference_Weight_H700/F");
        otree.Branch("interference_Weight_H800", interference_Weight_H800_ , "interference_Weight_H800/F");
        otree.Branch("interference_Weight_H900", interference_Weight_H900_ , "interference_Weight_H900/F");
        otree.Branch("interference_Weight_H1000", interference_Weight_H1000_ , "interference_Weight_H1000/F");
        otree.Branch("cps_Weight_H600", cps_Weight_H600_ , "cps_Weight_H600/F");
        otree.Branch("cps_Weight_H700", cps_Weight_H700_ , "cps_Weight_H700/F");
        otree.Branch("cps_Weight_H800", cps_Weight_H800_ , "cps_Weight_H800/F");
        otree.Branch("cps_Weight_H900", cps_Weight_H900_ , "cps_Weight_H900/F");
        otree.Branch("cps_Weight_H1000", cps_Weight_H1000_ , "cps_Weight_H1000/F");

        otree.Branch("jet_grsens_ft", jet_grsens_ft_ , "jet_grsens_ft/F");
        otree.Branch("jet_grsens_tr", jet_grsens_tr_ , "jet_grsens_tr/F");
        otree.Branch("jet_massdrop_pr", jet_massdrop_pr_ , "jet_massdrop_pr/F");
        otree.Branch("jet_qjetvol", jet_qjetvol_ , "jet_qjetvol/F");
        otree.Branch("jet_tau2tau1", jet_tau2tau1_ , "jet_tau2tau1/F");
        otree.Branch("jet_jetconstituents", jet_jetconstituents_ , "jet_jetconstituents/F");
        otree.Branch("jet_rcore4", jet_rcore4_ , "jet_rcore4/F");
        otree.Branch("jet_rcore5", jet_rcore5_ , "jet_rcore5/F");
        otree.Branch("jet_rcore6", jet_rcore6_ , "jet_rcore6/F");
        otree.Branch("jet_rcore7", jet_rcore7_ , "jet_rcore7/F");


        otree.Branch("jet_grsens_ft_exo", jet_grsens_ft_exo_ , "jet_grsens_ft_exo/F");
        otree.Branch("jet_grsens_tr_exo", jet_grsens_tr_exo_ , "jet_grsens_tr_exo/F");
        otree.Branch("jet_massdrop_pr_exo", jet_massdrop_pr_exo_ , "jet_massdrop_pr_exo/F");
        otree.Branch("jet_tau2tau1_exo", jet_tau2tau1_exo_ , "jet_tau2tau1_exo/F");

        otree.Branch("nbjets", nbjets_ , "nbjets/F");
        otree.Branch("nbjets_cvsl", nbjets_cvsl_ , "nbjets_cvsl/F");
        otree.Branch("nbjets_cvsm", nbjets_cvsm_ , "nbjets_cvsm/F");
        otree.Branch("nbjets_ssvhem", nbjets_ssvhem_ , "nbjets_ssvhem/F");

        otree.Branch("nbjets_csvl_veto", nbjets_csvl_veto_ , "nbjets_csvl_veto/F");
        otree.Branch("nbjets_csvm_veto", nbjets_csvm_veto_ , "nbjets_csvm_veto/F");        
        otree.Branch("nbjets_csvt_veto", nbjets_csvt_veto_ , "nbjets_csvt_veto/F");        
        otree.Branch("nbjets_ssvhem_veto", nbjets_ssvhem_veto_ , "nbjets_ssvhem_veto/F");                

        otree.Branch("nbjets_csvl_veto_cleaned",nbjets_csvl_veto_cleaned_, "nbjets_csvl_veto_cleaned/F");
        otree.Branch("nbjets_csvm_veto_cleaned", nbjets_csvm_veto_cleaned_ , "nbjets_csvm_veto_cleaned/F");        
        otree.Branch("nbjets_csvt_veto_cleaned", nbjets_csvt_veto_cleaned_ , "nbjets_csvt_veto_cleaned/F");        
        otree.Branch("nbjets_ssvhem_veto_cleaned", nbjets_ssvhem_veto_cleaned_ , "nbjets_ssvhem_veto_cleaned/F");                

        otree.Branch("nbjetsCSV", nbjetsCSV_ , "nbjetsCSV_/F");
        otree.Branch("nbjetsSSVHE", nbjetsSSVHE_ , "nbjetsSSVHE_/F");        
        otree.Branch("njets", njets_ , "njets/F");
        otree.Branch("jet_pt1frac", jet_pt1frac_ , "jet_pt1frac/F");
        otree.Branch("jet_pt2frac", jet_pt2frac_ , "jet_pt2frac/F");
        otree.Branch("jet_sjdr", jet_sjdr_ , "jet_sjdr/F");


        otree.Branch("jet_planarflow04", jet_planarlow04_, "jet_planarflow04/F");
        otree.Branch("jet_planarflow05", jet_planarlow05_, "jet_planarflow05/F");
        otree.Branch("jet_planarflow06", jet_planarlow06_, "jet_planarflow06/F");
        otree.Branch("jet_planarflow07", jet_planarlow07_, "jet_planarflow07/F");
        
        otree.Branch("deltaR_lca8jet", deltaR_lca8jet_, "deltaR_lca8jet/F");
        otree.Branch("deltaphi_METca8jet", deltaphi_METca8jet_, "deltaphi_METca8jet/F");
        otree.Branch("deltaphi_Vca8jet", deltaphi_Vca8jet_, "deltaphi_Vca8jet/F");
        

        ##########################     
        # variables for ttbar control region
        ttb_nak5_same_       = array( 'f', [ 0. ] );        
        ttb_nak5_same_csvl_  = array( 'f', [ 0. ] );        
        ttb_nak5_same_csvm_  = array( 'f', [ 0. ] );        
        ttb_nak5_same_csvt_  = array( 'f', [ 0. ] );        
        ttb_nak5_oppo_       = array( 'f', [ 0. ] );        
        ttb_nak5_oppo_csvl_  = array( 'f', [ 0. ] );        
        ttb_nak5_oppo_csvm_  = array( 'f', [ 0. ] );        
        ttb_nak5_oppo_csvt_  = array( 'f', [ 0. ] );        
        ttb_nak5_oppoveto_       = array( 'f', [ 0. ] );        
        ttb_nak5_oppoveto_csvl_  = array( 'f', [ 0. ] );        
        ttb_nak5_oppoveto_csvm_  = array( 'f', [ 0. ] );        
        ttb_nak5_oppoveto_csvt_  = array( 'f', [ 0. ] );        
        ttb_ht_  = array( 'f', [ 0. ] );                
        ttb_ca8_mass_pr_  = array( 'f', [ 0. ] );                        

        ttb_ca8_px_  = array( 'f', [ 0. ] );                        
        ttb_ca8_py_  = array( 'f', [ 0. ] );                                
        ttb_ca8_pz_  = array( 'f', [ 0. ] );                                
        ttb_ca8_e_  = array( 'f', [ 0. ] );                                
        
        ttb_ca8_ungroomed_pt_  = array( 'f', [ 0. ] );                        
        ttb_ca8_tau2tau1_  = array( 'f', [ 0. ] );                        
        ttb_ca8_mu_  = array( 'f', [ 0. ] );                        
        ttb_mlvj_  = array( 'f', [ 0. ] );                        
        
        isttbar_ = array( 'i', [ 0 ] );            
        
        otree.Branch("ttb_nak5_same", ttb_nak5_same_ , "ttb_nak5_same/F");        
        otree.Branch("ttb_nak5_same_csvl", ttb_nak5_same_csvl_ , "ttb_nak5_same_csvl/F");        
        otree.Branch("ttb_nak5_same_csvm", ttb_nak5_same_csvm_ , "ttb_nak5_same_csvm/F");        
        otree.Branch("ttb_nak5_same_csvt", ttb_nak5_same_csvt_ , "ttb_nak5_same_csvt/F");        
        otree.Branch("ttb_nak5_oppo", ttb_nak5_oppo_ , "ttb_nak5_oppo/F");                
        otree.Branch("ttb_nak5_oppo_csvl", ttb_nak5_oppo_csvl_ , "ttb_nak5_oppo_csvl/F");                
        otree.Branch("ttb_nak5_oppo_csvm", ttb_nak5_oppo_csvm_ , "ttb_nak5_oppo_csvm/F");        
        otree.Branch("ttb_nak5_oppo_csvt", ttb_nak5_oppo_csvt_ , "ttb_nak5_oppo_csvt/F");        
        otree.Branch("ttb_nak5_oppoveto", ttb_nak5_oppoveto_ , "ttb_nak5_oppoveto/F");        
        otree.Branch("ttb_nak5_oppoveto_csvl", ttb_nak5_oppoveto_csvl_ , "ttb_nak5_oppoveto_csvl/F");        
        otree.Branch("ttb_nak5_oppoveto_csvm", ttb_nak5_oppoveto_csvm_ , "ttb_nak5_oppoveto_csvm/F");        
        otree.Branch("ttb_nak5_oppoveto_csvt", ttb_nak5_oppoveto_csvt_ , "ttb_nak5_oppoveto_csvt/F");        
        
        otree.Branch("ttb_ht", ttb_ht_ , "ttb_ht_/F");                
        otree.Branch("ttb_ca8_mass_pr", ttb_ca8_mass_pr_ , "ttb_ca8_mass_pr/F");                

        otree.Branch("ttb_ca8_px", ttb_ca8_px_ , "ttb_ca8_px/F");                
        otree.Branch("ttb_ca8_py", ttb_ca8_py_ , "ttb_ca8_py/F");                
        otree.Branch("ttb_ca8_pz", ttb_ca8_pz_ , "ttb_ca8_pz/F");                
        otree.Branch("ttb_ca8_e", ttb_ca8_e_ , "ttb_ca8_e/F");                

        otree.Branch("ttb_ca8_ungroomed_pt", ttb_ca8_ungroomed_pt_ , "ttb_ca8_ungroomed_pt/F");                
        otree.Branch("ttb_ca8_tau2tau1", ttb_ca8_tau2tau1_ , "ttb_ca8_tau2tau1/F");                
        otree.Branch("ttb_ca8_mu", ttb_ca8_mu_ , "ttb_ca8_mu/F");                
        otree.Branch("ttb_mlvj", ttb_mlvj_ , "ttb_mlvj/F");                
        
        otree.Branch("isttbar", isttbar_ , "isttbar/I"); 
        
        ##########################
        
        ##########################
        # variables for BSM reweighting
        genHMass_ = array( 'f', [ 0. ] );              
        otree.Branch("genHMass", genHMass_ , "genHMass/F");        
        
        cprimeVals = [1,2,3,4,5,6,7,8,9,10]
        brnewVals = [00, 01, 02, 03, 04, 05]
        bsmReweights = [];        
        for i in range(len(cprimeVals)): 
            col_bsmReweights = [];
            for j in range(len(brnewVals)): 
                col_bsmReweights.append( array( 'f', [ 0. ] ) );
            bsmReweights.append( col_bsmReweights );

        for i in range(len(cprimeVals)): 
            for j in range(len(brnewVals)): 
                brname = "bsmReweight_cPrime%02d_brNew%02d"%(cprimeVals[i],brnewVals[j]);
                #print brname
                otree.Branch(brname, bsmReweights[i][j] , brname);        
        ##########################

        ##########################
        # variables for GEN
        gen_parton1_px_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton1_py_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton1_pz_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton1_e_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton1_id_fromttbar_ = array( 'i', [ 0 ] );
        gen_parton2_px_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton2_py_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton2_pz_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton2_e_fromttbar_ = array( 'f', [ 0. ] );
        gen_parton2_id_fromttbar_ = array( 'i', [ 0 ] );

        otree.Branch("gen_parton1_px_fromttbar", gen_parton1_px_fromttbar_ , "gen_parton1_px_fromttbar/F");        
        otree.Branch("gen_parton1_py_fromttbar", gen_parton1_py_fromttbar_ , "gen_parton1_py_fromttbar/F");        
        otree.Branch("gen_parton1_pz_fromttbar", gen_parton1_pz_fromttbar_ , "gen_parton1_pz_fromttbar/F");        
        otree.Branch("gen_parton1_e_fromttbar", gen_parton1_e_fromttbar_ , "gen_parton1_e_fromttbar/F");        
        otree.Branch("gen_parton1_id_fromttbar", gen_parton1_id_fromttbar_ , "gen_parton1_id_fromttbar/F");        

        otree.Branch("gen_parton2_px_fromttbar", gen_parton2_px_fromttbar_ , "gen_parton2_px_fromttbar/F");        
        otree.Branch("gen_parton2_py_fromttbar", gen_parton2_py_fromttbar_ , "gen_parton2_py_fromttbar/F");        
        otree.Branch("gen_parton2_pz_fromttbar", gen_parton2_pz_fromttbar_ , "gen_parton2_pz_fromttbar/F");        
        otree.Branch("gen_parton2_e_fromttbar", gen_parton2_e_fromttbar_ , "gen_parton2_e_fromttbar/F");        
        otree.Branch("gen_parton2_id_fromttbar", gen_parton2_id_fromttbar_ , "gen_parton2_id_fromttbar/F");        

        ##########################
        

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
        for i in range(NLoop):
#        for i in range(2):            
            
            if i % 10000 == 0: 
                print "i = ", i

            self.InputTree_.GetEntry(i);


            # ---------- ttbar control region

            ttbarlike = 0;

            if self.Channel_ == 'mu': lepLabel = "muon";
            if self.Channel_ == 'el': lepLabel = "electron";
            index_ca8_in_oppoHemi = [];
            for i in range(6):
                if getattr( self.InputTree_, "GroomedJet_CA8_pt" )[i] > 200:
                    j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[i]
                    j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[i]
                    l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" )
                    l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" )
                    dR_lj = math.sqrt( (l_eta - j_ca8_eta)**2 + (l_phi - j_ca8_phi)**2 );
                    if dR_lj > ROOT.TMath.Pi()/2.: index_ca8_in_oppoHemi.append(i);

            minMass = -1;
            theca8Index = -1;
            for i in range(len(index_ca8_in_oppoHemi)):
                curmass = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[index_ca8_in_oppoHemi[i]]
                if curmass > minMass: 
                    minMass = curmass;
                    theca8Index = index_ca8_in_oppoHemi[i];

            index_ak5_in_sameHemi = [];
            index_ak5_in_oppoHemi = [];
            index_ak5_in_oppoHemi_vetoca8 = [];
            index_ak5_in_sameHemi_csvl = [];
            index_ak5_in_oppoHemi_csvl = [];
            index_ak5_in_oppoHemi_vetoca8_csvl = [];
            index_ak5_in_sameHemi_csvm = [];
            index_ak5_in_oppoHemi_csvm = [];
            index_ak5_in_oppoHemi_vetoca8_csvm = [];
            index_ak5_in_sameHemi_csvt = [];
            index_ak5_in_oppoHemi_csvt = [];
            index_ak5_in_oppoHemi_vetoca8_csvt = [];
            ttb_ht = getattr( self.InputTree_, "W_"+lepLabel+"_pt" );
            ttb_ht += getattr( self.InputTree_, "event_met_pfmet" ); 
        
            if theca8Index >= 0:
                for i in range(6):
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[i] > 30:
                        ttb_ht += getattr( self.InputTree_, "JetPFCor_Pt" )[i]
                        
                        j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index]
                        j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index]
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[i]
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[i]
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" )
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" )                
                        dR_jj = math.sqrt( (j_ak5_eta - j_ca8_eta)**2 + (j_ak5_phi - j_ca8_phi)**2 );
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + (l_phi - j_ak5_phi)**2 );
                        
                        if dR_lj < ROOT.TMath.Pi()/2.: 
                            index_ak5_in_sameHemi.append( i );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_sameHemi_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_sameHemi_csvm.append(i);                        
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_sameHemi_csvt.append(i);                        
                        else: 
                            index_ak5_in_oppoHemi.append( i );    
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_oppoHemi_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_oppoHemi_csvm.append(i);                                            
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_oppoHemi_csvt.append(i);                                            
                        if dR_lj > ROOT.TMath.Pi()/2. and dR_jj > 0.8: 
                            index_ak5_in_oppoHemi_vetoca8.append( i );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.244: index_ak5_in_oppoHemi_vetoca8_csvl.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.679: index_ak5_in_oppoHemi_vetoca8_csvm.append(i);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_in_oppoHemi_vetoca8_csvt.append(i);
                            

                ttb_nak5_same_[0] = len(index_ak5_in_sameHemi);
                ttb_nak5_same_csvl_[0] = len(index_ak5_in_sameHemi_csvl);
                ttb_nak5_same_csvm_[0] = len(index_ak5_in_sameHemi_csvm);
                ttb_nak5_same_csvt_[0] = len(index_ak5_in_sameHemi_csvt);
                ttb_nak5_oppo_[0] = len(index_ak5_in_oppoHemi);
                ttb_nak5_oppo_csvl_[0] = len(index_ak5_in_oppoHemi_csvl);
                ttb_nak5_oppo_csvm_[0] = len(index_ak5_in_oppoHemi_csvm);
                ttb_nak5_oppo_csvt_[0] = len(index_ak5_in_oppoHemi_csvt);
                ttb_nak5_oppoveto_[0] = len(index_ak5_in_oppoHemi_vetoca8);
                ttb_nak5_oppoveto_csvl_[0] = len(index_ak5_in_oppoHemi_vetoca8_csvl);
                ttb_nak5_oppoveto_csvm_[0] = len(index_ak5_in_oppoHemi_vetoca8_csvm);
                ttb_nak5_oppoveto_csvt_[0] = len(index_ak5_in_oppoHemi_vetoca8_csvt);
                ttb_ht_[0] = ttb_ht;
                ttb_ca8_mass_pr_[0] = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[theca8Index];
                ttb_ca8_ungroomed_pt_[0] = getattr( self.InputTree_, "GroomedJet_CA8_pt" )[theca8Index];
            
                ttb_ca8_tau2tau1_[0] = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1" )[theca8Index];
                ttb_ca8_mu_[0] = getattr( self.InputTree_, "GroomedJet_CA8_massdrop_pr" )[theca8Index];

                ttb_ca8J_p4 = ROOT.TLorentzVector();        
                ttb_ca8J_pt = getattr( self.InputTree_, "GroomedJet_CA8_pt_pr" )[theca8Index];    
                ttb_ca8J_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta_pr" )[theca8Index];    
                ttb_ca8J_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi_pr" )[theca8Index];    
                ttb_ca8J_e = getattr( self.InputTree_, "GroomedJet_CA8_e_pr" )[theca8Index];                        
                ttb_ca8J_p4.SetPtEtaPhiE(ttb_ca8J_pt, ttb_ca8J_eta, ttb_ca8J_phi, ttb_ca8J_e)
                ttb_V_p4 = ROOT.TLorentzVector(self.InputTree_.W_px,self.InputTree_.W_py,self.InputTree_.W_pz,self.InputTree_.W_e);                    
                ttb_mlvj_[0] = (ttb_V_p4+ttb_ca8J_p4).M();
            
                ttb_ca8_px_[0] = ttb_ca8J_p4.Px();
                ttb_ca8_py_[0] = ttb_ca8J_p4.Py();
                ttb_ca8_pz_[0] = ttb_ca8J_p4.Pz();
                ttb_ca8_e_[0] =  ttb_ca8J_p4.E();

                oppo1same1 = ttb_nak5_same_csvl_[0] > 0 or ttb_nak5_oppo_csvl_[0] > 0;
                oppo2same0 = ttb_nak5_same_csvl_[0] == 0 and ttb_nak5_oppo_csvl_[0] > 1;
                oppo2same0 = False
                if oppo1same1 or oppo2same0:
                    isttbar_[0] = 1
                    
            else:

                ttb_nak5_same_[0] = -1;
                ttb_nak5_same_csvl_[0] = -1;
                ttb_nak5_same_csvm_[0] = -1;
                ttb_nak5_same_csvt_[0] = -1;
                ttb_nak5_oppo_[0] = -1;
                ttb_nak5_oppo_csvl_[0] = -1;
                ttb_nak5_oppo_csvm_[0] = -1;
                ttb_nak5_oppo_csvt_[0] = -1;
                ttb_nak5_oppoveto_[0] = -1;
                ttb_nak5_oppoveto_csvl_[0] = -1;
                ttb_nak5_oppoveto_csvm_[0] = -1;
                ttb_nak5_oppoveto_csvt_[0] = -1;
                ttb_ht_[0] = -1;
                ttb_ca8_mass_pr_[0] = -1;
                ttb_ca8_ungroomed_pt_[0] = -1;
                ttb_ca8_tau2tau1_[0] = -1;
                ttb_ca8_mu_[0] = -1;
                ttb_mlvj_[0] = -1;
                isttbar_[0] = 0
                        
                ttb_ca8_px_[0] = 0;
                ttb_ca8_py_[0] = 0;
                ttb_ca8_pz_[0] = 0;
                ttb_ca8_e_[0] = 0;


            # ---------- end ttbar control region

            #### Exo CA8 Jet ID
            index_ca8jet_exo_                  = 0;
            mass_pr_ca8jet_exo_                = 0.;
            deltaphi_met_ca8jet_exo_           = 0. ;
            deltaphi_lep_ca8jet_exo_           = 0;
            GroomedJet_CA8_deltaR_lca8jet_exo_ = 0;
            W_mass = 80.38;
            
            for i in range(6):
                if( math.fabs( getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[i] - W_mass) < math.fabs( mass_pr_ca8jet_exo_ - W_mass)):
                    index_ca8jet_exo_=i ;
                    mass_pr_ca8jet_exo_ = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[index_ca8jet_exo_] ;
  
            if(math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_phi" )[index_ca8jet_exo_]-getattr( self.InputTree_, "event_met_pfmetPhi" ))<=ROOT.TMath.Pi()):
              deltaphi_met_ca8jet_exo_  = math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_phi" )[index_ca8jet_exo_]-getattr( self.InputTree_, "event_met_pfmetPhi" ));
            else : deltaphi_met_ca8jet_exo_ = 2*ROOT.TMath.Pi()-math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_phi" )[index_ca8jet_exo_]-getattr( self.InputTree_,"event_met_pfmetPhi" ));

            if(math.fabs(getattr( self.InputTree_,"GroomedJet_CA8_phi")[index_ca8jet_exo_]-getattr( self.InputTree_,"W_"+lepLabel+"_phi"))<=ROOT.TMath.Pi()):
              deltaphi_lep_ca8jet_exo_  = math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_phi" )[index_ca8jet_exo_]-getattr( self.InputTree_, "W_"+lepLabel+"_phi" ));
            else : deltaphi_lep_ca8jet_exo_ = 2*ROOT.TMath.Pi()-math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_phi" )[index_ca8jet_exo_]-getattr( self.InputTree_,"W_"+lepLabel+"_phi"));

            GroomedJet_CA8_deltaR_lca8jet_exo_ = math.sqrt((getattr(self.InputTree_,"GroomedJet_CA8_eta")[i] -getattr( self.InputTree_,"W_"+lepLabel+"_eta"))**2+(deltaphi_lep_ca8jet_exo_)**2 );

                   
            
            # make cuts
            leptonCut = 30;
            leptonCutString = "W_muon_pt";
            metCut = 40;
            if self.Channel_ == "el": 
                leptonCut = 35;    
                leptonCutString = "W_electron_pt";
                metCut = 40;

            signallike = 0;
            if getattr( self.InputTree_, "W_pt" ) > 200 and getattr( self.InputTree_, "GroomedJet_CA8_pt" )[0] > 200 and getattr( self.InputTree_, "event_met_pfmet" ) > metCut and getattr( self.InputTree_, leptonCutString ) > leptonCut and getattr( self.InputTree_, "GroomedJet_CA8_deltaphi_METca8jet_type2") > 2.0 and getattr( self.InputTree_, "GroomedJet_CA8_deltaR_lca8jet") > 1.57 :
                if ( self.Channel_ == "mu" and math.fabs(getattr( self.InputTree_, "W_muon_dz000")) < 0.02 and math.fabs(getattr( self.InputTree_, "W_muon_dzPV")) < 0.5 and math.fabs( getattr( self.InputTree_, "W_muon_eta" )) < 2.1 ) : signallike = 1 ; 
                elif self.Channel_ == "el" : signallike = 1 ; 
                              
            signallike_exo_ = 0;
            if getattr( self.InputTree_, "W_pt" ) > 200 and   getattr( self.InputTree_, "GroomedJet_CA8_pt" )[index_ca8jet_exo_] > 200 and getattr( self.InputTree_, "event_met_pfmet" ) > metCut and getattr( self.InputTree_, leptonCutString ) > leptonCut and deltaphi_met_ca8jet_exo_ > 2.0 and GroomedJet_CA8_deltaR_lca8jet_exo_> 1.57 :signallike_exo_ = 1;

            

            ttbarlike = 0;
            if isttbar_[0] == 1 and getattr( self.InputTree_, "event_met_pfmet" ) > metCut and getattr( self.InputTree_, leptonCutString ) > leptonCut:
                ttbarlike = 1;

            if ttbarlike == 1 or signallike == 1 or signallike_exo_ ==1:                
             
             event_[0] = self.InputTree_.event_evtNo;
             event_runNo_[0]      = getattr( self.InputTree_, "event_runNo" );
             event_lumi_[0]       = getattr( self.InputTree_, "event_lumi" );

             
             effwt = getattr( self.InputTree_, "effwt" );
             puwt = getattr( self.InputTree_, "puwt" );
             totSampleWeight = 1.;
             if self.IsData_: totSampleWeight = wSampleWeight;
             else: totSampleWeight = wSampleWeight*effwt*puwt;
                #print puwt;

             rwCPS = 1;
             if self.SignalMass_ > 0:
                 binVal = self.x_rwCPS.FindBin(self.InputTree_.W_H_mass_gen);
                 if binVal > self.h_rwCPS.GetNbinsX(): binVal = self.h_rwCPS.GetNbinsX();
                 if binVal < 1: binVal = 1;
                    #print "binVal: ",binVal
                 rwCPS = self.h_rwCPS.GetBinContent( binVal );

             #interference weight
             complexpolewtggH600    = getattr(self.InputTree_,"complexpolewtggH600")*rwCPS;
             interferencewtggH600   = getattr(self.InputTree_,"interferencewtggH600");
             avecomplexpolewtggH600 = getattr(self.InputTree_,"avecomplexpolewtggH600"); 
             infe_Weight_H600 = complexpolewtggH600*interferencewtggH600/avecomplexpolewtggH600;

             complexpolewtggH700    = getattr(self.InputTree_,"complexpolewtggH700")*rwCPS; 
             interferencewtggH700   = getattr(self.InputTree_,"interferencewtggH700");
             avecomplexpolewtggH700 = getattr(self.InputTree_,"avecomplexpolewtggH700"); 
             infe_Weight_H700 = complexpolewtggH700*interferencewtggH700/avecomplexpolewtggH700;

             complexpolewtggH800    = getattr(self.InputTree_,"complexpolewtggH800")*rwCPS; 
             interferencewtggH800   = getattr(self.InputTree_,"interferencewtggH800");
             avecomplexpolewtggH800 = getattr(self.InputTree_,"avecomplexpolewtggH800"); 
             infe_Weight_H800 = complexpolewtggH800*interferencewtggH800/avecomplexpolewtggH800;

             complexpolewtggH900    = getattr(self.InputTree_,"complexpolewtggH900")*rwCPS; 
             interferencewtggH900   = getattr(self.InputTree_,"interferencewtggH900");
             avecomplexpolewtggH900 = getattr(self.InputTree_,"avecomplexpolewtggH900"); 
             infe_Weight_H900 = complexpolewtggH900*interferencewtggH900/avecomplexpolewtggH900;

             complexpolewtggH1000    = getattr(self.InputTree_,"complexpolewtggH1000")*rwCPS; 
             interferencewtggH1000   = getattr(self.InputTree_,"interferencewtggH1000");
             avecomplexpolewtggH1000 = getattr(self.InputTree_,"avecomplexpolewtggH1000"); 
             infe_Weight_H1000 = complexpolewtggH1000*interferencewtggH1000/avecomplexpolewtggH1000;

              # produce weights for alternative models
             if self.SignalMass_ > 0:
#                if False:                    
                    curIntfRw = getattr(self.InputTree_,"interferencewtggH%03d"%(self.SignalMass_));
                    
                    
                    genHMass_[0] = getattr(self.InputTree_,"W_H_mass_gen");
                    for i in range(len(cprimeVals)): 
                        for j in range(len(brnewVals)): 
                            curCprime = float(cprimeVals[i])/10.;
                            curBRnew = float(brnewVals[j])/10.;
                            if self.isVBF_: 
                                bsmReweights[i][j][0] = self.GetInteferenceWeights( getattr(self.InputTree_,"W_H_mass_gen"), curCprime, curBRnew );  
                            else:
                                bsmReweights[i][j][0] = self.GetInteferenceWeights( getattr(self.InputTree_,"W_H_mass_gen"), curCprime, curBRnew )*IntfRescale(curIntfRw,curCprime,curBRnew);
                
             else:   
                    genHMass_[0] = -1;
                    for i in range(len(cprimeVals)): 
                        for j in range(len(brnewVals)): 
                            bsmReweights[i][j][0] = -1;
            
                
             ###################################
             # make training tree
             mass_lvj_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type2" );

             v_pt_[0] = getattr( self.InputTree_, "W_pt" );
             v_mt_[0] = getattr( self.InputTree_, "W_mt" );

             nu_pz_type0_[0] = getattr( self.InputTree_, "W_nu1_pz_type0" );
             nu_pz_type1_[0] = getattr( self.InputTree_, "W_nu1_pz_type1" );
             nu_pz_type2_[0] = getattr( self.InputTree_, "W_nu1_pz_type2" );
             nu_pz_type3_[0] = getattr( self.InputTree_, "W_nu1_pz_type3" );

             W_pz_type0_[0] = getattr( self.InputTree_, "W_pz_type0" );
             W_pz_type1_[0] = getattr( self.InputTree_, "W_pz_type1" );
             W_pz_type2_[0] = getattr( self.InputTree_, "W_pz_type2" );
             W_pz_type3_[0] = getattr( self.InputTree_, "W_pz_type3" );
             
             if self.IsData_ == False  : nu_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" );
             if self.Channel_ == "mu" and self.IsData_ == False: W_pz_gen_[0]   =  getattr( self.InputTree_, "W_neutrino_pz_gen" ) + getattr( self.InputTree_, "W_muon_pz_gen" );
             if self.Channel_ == "el" and self.IsData_ == False: W_pz_gen_[0]   =  getattr( self.InputTree_, "W_neutrino_pz_gen" ) + getattr( self.InputTree_, "W_electron_pz_gen" );

              
             mass_lvj_type0_[0]   = getattr( self.InputTree_, "boostedW_lvj_m_type0" );
             mass_lvj_type1_[0]   = getattr( self.InputTree_, "boostedW_lvj_m_type1" );
             mass_lvj_type2_[0]   = getattr( self.InputTree_, "boostedW_lvj_m_type2" );
             mass_lvj_type3_[0]   = getattr( self.InputTree_, "boostedW_lvj_m_type3" );

             v_mass_type0_[0]    = getattr( self.InputTree_, "W_mass_type0" );
             v_mass_type1_[0]    = getattr( self.InputTree_, "W_mass_type1" );
             v_mass_type2_[0]    = getattr( self.InputTree_, "W_mass_type2" );
             v_mass_type3_[0]    = getattr( self.InputTree_, "W_mass_type3" );
             
             jet_mass_pr_[0] = getattr( self.InputTree_, prefix + "_mass_pr" )[0];
               
             jet_mass_pr_exo_[0] = getattr( self.InputTree_, prefix + "_mass_pr" )[index_ca8jet_exo_];
                 
             issignal_[0] = signallike;
             issignal_exo_[0] = signallike_exo_;

             jet_pt_pr_[0] = getattr( self.InputTree_, prefix + "_pt_pr" )[0];
             ungroomed_jet_pt_[0] = getattr( self.InputTree_, prefix+"_pt" )[0];
             ungroomed_jet_eta_[0] = getattr( self.InputTree_, prefix+"_eta" )[0];
             ungroomed_jet_phi_[0] = getattr( self.InputTree_, prefix+"_phi" )[0];

             jet_pt_pr_exo_[0] = getattr( self.InputTree_, prefix + "_pt_pr" )[index_ca8jet_exo_];
             ungroomed_jet_pt_exo_[0] = getattr( self.InputTree_, prefix+"_pt" )[index_ca8jet_exo_];

             self.jecUnc_.setJetEta( getattr( self.InputTree_, prefix+"_eta" )[0] )
             self.jecUnc_.setJetPt( getattr( self.InputTree_, prefix+"_pt" )[0] );                        
             j_jecfactor_up_[0] = self.jecUnc_.getUncertainty( True )
             self.jecUnc_.setJetEta( getattr( self.InputTree_, prefix+"_eta" )[0] )
             self.jecUnc_.setJetPt( getattr( self.InputTree_, prefix+"_pt" )[0] );               
             j_jecfactor_dn_[0] = self.jecUnc_.getUncertainty( False ) 

             j_jecfactor_up_[0] = math.sqrt( j_jecfactor_up_[0]**2 + 0.02**2 );
             j_jecfactor_dn_[0] = math.sqrt( j_jecfactor_dn_[0]**2 + 0.02**2 );
             curjes_up = 1 + j_jecfactor_up_[0]
             curjes_dn = 1 - j_jecfactor_dn_[0]

             jorig_pt = getattr( self.InputTree_, prefix + "_pt_pr" )[0];    
             jorig_eta = getattr( self.InputTree_, prefix + "_eta_pr" )[0];    
             jorig_phi = getattr( self.InputTree_, prefix + "_phi_pr" )[0];    
             jorig_e = getattr( self.InputTree_, prefix + "_e_pr" )[0];                        
             jdef_ptetaphie = ROOT.TLorentzVector();
             jdef_ptetaphie.SetPtEtaPhiE(jorig_pt, jorig_eta, jorig_phi, jorig_e)
             jdef_up = ROOT.TLorentzVector(jdef_ptetaphie.Px() * curjes_up, jdef_ptetaphie.Py() * curjes_up, jdef_ptetaphie.Pz() * curjes_up, jdef_ptetaphie.E() * curjes_up)
             jdef_dn = ROOT.TLorentzVector(jdef_ptetaphie.Px() * curjes_dn, jdef_ptetaphie.Py() * curjes_dn, jdef_ptetaphie.Pz() * curjes_dn, jdef_ptetaphie.E() * curjes_dn)
             jet_mass_pr_up_[0] = jdef_up.M();  
             jet_mass_pr_dn_[0] = jdef_dn.M();  
                    
             if self.Channel_ == "mu" :
                    l_pt_[0] = getattr( self.InputTree_, "W_muon_pt" );
                    l_eta_[0] = getattr( self.InputTree_, "W_muon_eta" );
                    l_phi_[0] = getattr( self.InputTree_, "W_muon_phi" );
                    
             elif self.Channel_ == "el":
                    l_pt_[0] = getattr( self.InputTree_, "W_electron_pt" );
                    l_eta_[0] = getattr( self.InputTree_, "W_electron_eta" );
                    l_phi_[0] = getattr( self.InputTree_, "W_electron_phi" );

             mvaMET_[0] = getattr( self.InputTree_, "event_metMVA_met" );
             pfMET_[0] = getattr( self.InputTree_, "event_met_pfmet" );        
             pfMET_Phi_[0] = getattr( self.InputTree_, "event_met_pfmetPhi" );        

             nPV_[0] = getattr( self.InputTree_, "event_nPV" );
             totalEventWeight_[0] = totSampleWeight;
             eff_and_pu_Weight_[0] = effwt*puwt;
             wSampleWeight_[0]     = wSampleWeight;

             interference_Weight_H600_[0] = infe_Weight_H600;
             interference_Weight_H700_[0] = infe_Weight_H700;
             interference_Weight_H800_[0] = infe_Weight_H800;
             interference_Weight_H900_[0] = infe_Weight_H900;
             interference_Weight_H1000_[0] = infe_Weight_H1000;
             cps_Weight_H600_[0] = complexpolewtggH600/avecomplexpolewtggH600;
             cps_Weight_H700_[0] = complexpolewtggH700/avecomplexpolewtggH700;
             cps_Weight_H800_[0] = complexpolewtggH800/avecomplexpolewtggH800;
             cps_Weight_H900_[0] = complexpolewtggH900/avecomplexpolewtggH900;
             cps_Weight_H1000_[0] = complexpolewtggH1000/avecomplexpolewtggH1000;

             jet_grsens_ft_[0] = getattr( self.InputTree_, prefix + "_mass_ft" )[0] / getattr( self.InputTree_, prefix + "_mass" )[0];
             jet_grsens_tr_[0] = getattr( self.InputTree_, prefix + "_mass_tr" )[0] / getattr( self.InputTree_, prefix + "_mass" )[0];
             jet_massdrop_pr_[0] = getattr( self.InputTree_, prefix + "_massdrop_pr" )[0];    

             jet_grsens_ft_exo_[0] = getattr( self.InputTree_, prefix + "_mass_ft" )[index_ca8jet_exo_] / getattr( self.InputTree_, prefix + "_mass" )[index_ca8jet_exo_];
             jet_grsens_tr_exo_[0] = getattr( self.InputTree_, prefix + "_mass_tr" )[index_ca8jet_exo_] / getattr( self.InputTree_, prefix + "_mass" )[index_ca8jet_exo_];
             jet_massdrop_pr_exo_[0] = getattr( self.InputTree_, prefix + "_massdrop_pr" )[index_ca8jet_exo_];    

             qjetmassdistribution = getattr( self.InputTree_, prefix+"_qjetmass" );
             qjetvol = getListRMS(qjetmassdistribution)/getListMean(qjetmassdistribution);
             jet_qjetvol_[0] = qjetvol;

             jet_tau2tau1_exo_[0] = getattr( self.InputTree_, prefix + "_tau2tau1" )[index_ca8jet_exo_];     
             jet_jetconstituents_exo_[0] = getattr( self.InputTree_, prefix + "_jetconstituents" )[index_ca8jet_exo_];     

             jet_tau2tau1_[0] = getattr( self.InputTree_, prefix + "_tau2tau1" )[0];     
             jet_jetconstituents_[0] = getattr( self.InputTree_, prefix + "_jetconstituents" )[0];     

             jet_rcore4_[0] = getattr( self.InputTree_, prefix + "_rcores")[3*6 + 0];
             jet_rcore5_[0] = getattr( self.InputTree_, prefix + "_rcores")[4*6 + 0];
             jet_rcore6_[0] = getattr( self.InputTree_, prefix + "_rcores")[5*6 + 0];
             jet_rcore7_[0] = getattr( self.InputTree_, prefix + "_rcores")[6*6 + 0];

             jet_planarlow04_[0] = getattr( self.InputTree_, prefix + "_planarflow04");
             jet_planarlow05_[0] = getattr( self.InputTree_, prefix + "_planarflow05");
             jet_planarlow06_[0] = getattr( self.InputTree_, prefix + "_planarflow06");
             jet_planarlow07_[0] = getattr( self.InputTree_, prefix + "_planarflow07");
                
#                nbjets_[0] = getattr( self.InputTree_, "GroomedJet_numberbjets" );

             nbjets_cvsl_[0]   = getattr( self.InputTree_, "GroomedJet_numberbjets_csvl" );
             nbjets_cvsm_[0]   = getattr( self.InputTree_, "GroomedJet_numberbjets_csvm" );
             nbjets_ssvhem_[0] = getattr( self.InputTree_, "GroomedJet_numberbjets_ssvhem" );

             nbjets_csvl_veto_cleaned_[0]   = getattr( self.InputTree_, "GroomedJet_numberbjets_csvl_veto" );
             nbjets_csvm_veto_cleaned_[0]   = getattr( self.InputTree_, "GroomedJet_numberbjets_csvm_veto" );
             nbjets_ssvhem_veto_cleaned_[0] = getattr( self.InputTree_, "GroomedJet_numberbjets_ssvhem_veto" );

             index_ak5_cvst = array( 'f', [0.] );

             nbjets_csvt_veto_cleaned_[0] = 0. ;
             dR_lj = 0. ;
             j_ca8_eta = 0 ; j_ca8_phi = 0;

             if getattr( self.InputTree_, "GroomedJet_CA8_pt" )[0] > 200:
                    j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[0]
 
                    j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[0]
                    l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" )
                    l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" )
                    dR_lj = math.sqrt( (l_eta - j_ca8_eta)**2 + (l_phi - j_ca8_phi)**2 );

             l_phi = 0. ; dR_jj = 0. ; dR_lj = 0. ;

             for i in range(6):
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[i] > 30:

                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" )                
                        dR_jj = math.sqrt( (j_ak5_eta - j_ca8_eta)**2 + (j_ak5_phi - j_ca8_phi)**2 );
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + (l_phi - j_ak5_phi)**2 );
                        
                        if dR_jj > 0.8: 
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] > 0.898: index_ak5_cvst.append(i);

               
             nbjets_csvt_veto_cleaned_[0] = len(index_ak5_in_oppoHemi_csvt);


             nbjetsCSV_[0] =0 ;
             nbjets_csvl_veto_[0] = 0 ;
             nbjets_csvm_veto_[0] = 0 ;
             nbjets_csvt_veto_[0] = 0 ;

             for i in range(0,6):
                    if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.244: nbjetsCSV_[0]=nbjetsCSV_[0]+1;
                    if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.244: nbjets_csvl_veto_[0]=nbjets_csvl_veto_[0]+1;
                    if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.679: nbjets_csvm_veto_[0]=nbjets_csvm_veto_[0]+1;
                    if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[i] >=0.898: nbjets_csvt_veto_[0]=nbjets_csvt_veto_[0]+1;

                    
             nbjetsSSVHE_[0]      = getattr( self.InputTree_, "numPFCorJetBTags");
             nbjets_ssvhem_veto_  = getattr( self.InputTree_, "numPFCorJetBTags");

                
             njets_[0] = getattr( self.InputTree_, "GroomedJet_numberjets" );
             pt1FracVal = max( getattr( self.InputTree_, prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, prefix + "_prsubjet2ptoverjetpt" ) );
             pt2FracVal = min( getattr( self.InputTree_, prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, prefix + "_prsubjet2ptoverjetpt" ) );
             jet_pt1frac_[0] = pt1FracVal;
             jet_pt2frac_[0] = pt2FracVal;
             jet_sjdr_[0] = getattr( self.InputTree_, prefix + "_prsubjet1subjet2_deltaR" );       
                
             deltaR_lca8jet_[0] = getattr( self.InputTree_, prefix + "_deltaR_lca8jet" );       
             deltaphi_METca8jet_[0] = getattr( self.InputTree_, prefix + "_deltaphi_METca8jet_type2" );       
             deltaphi_Vca8jet_[0] = getattr( self.InputTree_, prefix + "_deltaphi_Vca8jet_type2" );       
                
             listOfVarArray1[0][0] = jet_qjetvol_[0]
             listOfVarArray1[1][0] = jet_tau2tau1_[0]

             listOfVarArray2[0][0] = jet_qjetvol_[0]
             listOfVarArray2[1][0] = jet_tau2tau1_[0]


             otree.Fill();
             ###################################

        self.OFile_.cd();
        otree.Write();
        self.OFile_.Close();

    ## training tree name
    def getTrainingTreeName(self):
        return self.OFileName_
    
    def getSampleLabel(self):
        return self.Label_

    ### ------------------------------------------

    def turnOffBranches(self):
        
        prefix = self.JetPrefix_;
        
        self.InputTree_.SetBranchStatus("*",0);

     
        self.InputTree_.SetBranchStatus("event_evtNo",1);
        self.InputTree_.SetBranchStatus("event_runNo",1);
        self.InputTree_.SetBranchStatus("event_lumi",1);
        self.InputTree_.SetBranchStatus("ggdboostedWevt",1);
        self.InputTree_.SetBranchStatus("effwt",1);
        self.InputTree_.SetBranchStatus("puwt",1);
        self.InputTree_.SetBranchStatus("puwt_up",1);
        self.InputTree_.SetBranchStatus("puwt_down",1);

        if self.Channel_ == "mu" :
            self.InputTree_.SetBranchStatus("W_muon_pt",1);
            self.InputTree_.SetBranchStatus("W_muon_eta",1);
            self.InputTree_.SetBranchStatus("W_muon_phi",1);
            self.InputTree_.SetBranchStatus("W_muon_dz000");
            self.InputTree_.SetBranchStatus("W_muon_dzPV");
         
        elif self.Channel_ == "el":
            self.InputTree_.SetBranchStatus("W_electron_pt",1);
            self.InputTree_.SetBranchStatus("W_electron_eta",1);
            self.InputTree_.SetBranchStatus("W_electron_phi",1);


        self.InputTree_.SetBranchStatus("event_metMVA_met",1);
        self.InputTree_.SetBranchStatus("event_met_pfmet",1);
        self.InputTree_.SetBranchStatus("event_met_pfmetPhi",1);
        self.InputTree_.SetBranchStatus("event_nPV",1);
        self.InputTree_.SetBranchStatus("boostedW_lvj_m_type0",1);
        self.InputTree_.SetBranchStatus("boostedW_lvj_m_type1",1);
        self.InputTree_.SetBranchStatus("boostedW_lvj_m_type2",1);
        self.InputTree_.SetBranchStatus("boostedW_lvj_m_type3",1);
        self.InputTree_.SetBranchStatus("W_pt",1);
        self.InputTree_.SetBranchStatus("W_px",1);
        self.InputTree_.SetBranchStatus("W_py",1);
        self.InputTree_.SetBranchStatus("W_pz",1);
        self.InputTree_.SetBranchStatus("W_e",1);
        self.InputTree_.SetBranchStatus("W_mt",1);

        self.InputTree_.SetBranchStatus("W_nu1_pz_type0",1);
        self.InputTree_.SetBranchStatus("W_nu1_pz_type1",1);
        self.InputTree_.SetBranchStatus("W_nu1_pz_type2",1);
        self.InputTree_.SetBranchStatus("W_nu1_pz_type3",1);

        self.InputTree_.SetBranchStatus("W_pz_type0",1);
        self.InputTree_.SetBranchStatus("W_pz_type1",1);
        self.InputTree_.SetBranchStatus("W_pz_type2",1);
        self.InputTree_.SetBranchStatus("W_pz_type3",1);

        self.InputTree_.SetBranchStatus("W_neutrino_pz_gen",1);

        if self.Channel_ == "mu" : self.InputTree_.SetBranchStatus("W_muon_pz_gen",1);
        if self.Channel_ == "el" : self.InputTree_.SetBranchStatus("W_electron_pz_gen",1);
 

        self.InputTree_.SetBranchStatus("W_mass_type0",1);
        self.InputTree_.SetBranchStatus("W_mass_type1",1);
        self.InputTree_.SetBranchStatus("W_mass_type2",1);
        self.InputTree_.SetBranchStatus("W_mass_type3",1);

        self.InputTree_.SetBranchStatus(prefix + "_pt_pr",1);
        self.InputTree_.SetBranchStatus(prefix + "_eta_pr",1);
        self.InputTree_.SetBranchStatus(prefix + "_phi_pr",1);
        self.InputTree_.SetBranchStatus(prefix + "_e_pr",1);
        self.InputTree_.SetBranchStatus(prefix + "_eta",1);
        self.InputTree_.SetBranchStatus(prefix + "_phi",1);
        self.InputTree_.SetBranchStatus(prefix + "_pt",1);
        self.InputTree_.SetBranchStatus(prefix + "_e",1);
        self.InputTree_.SetBranchStatus(prefix + "_massdrop_pr",1);
        self.InputTree_.SetBranchStatus(prefix + "_mass",1);
        self.InputTree_.SetBranchStatus(prefix + "_mass_pr",1);        
        self.InputTree_.SetBranchStatus(prefix + "_mass_tr",1);
        self.InputTree_.SetBranchStatus(prefix + "_mass_ft",1);
        self.InputTree_.SetBranchStatus(prefix + "_tau2tau1",1);
        self.InputTree_.SetBranchStatus(prefix + "_qjetmass",1);
        self.InputTree_.SetBranchStatus(prefix + "_rcores",1);
        self.InputTree_.SetBranchStatus(prefix + "_jetconstituents",1);

        self.InputTree_.SetBranchStatus("GroomedJet_numberbjets_csvl",1);
        self.InputTree_.SetBranchStatus("GroomedJet_numberbjets_csvm",1);
        self.InputTree_.SetBranchStatus("GroomedJet_numberbjets_ssvhem",1);
        self.InputTree_.SetBranchStatus("GroomedJet_numberbjets_csvl_veto",1);
        self.InputTree_.SetBranchStatus("GroomedJet_numberbjets_csvm_veto",1);
        self.InputTree_.SetBranchStatus("GroomedJet_numberbjets_ssvhem_veto",1);

        self.InputTree_.SetBranchStatus("JetPFCor_Pt",1);
        self.InputTree_.SetBranchStatus("JetPFCor_Eta",1);
        self.InputTree_.SetBranchStatus("JetPFCor_Phi",1);
        self.InputTree_.SetBranchStatus("JetPFCor_bDiscriminatorCSV",1);
        self.InputTree_.SetBranchStatus("numPFCorJetBTags",1);
        self.InputTree_.SetBranchStatus(prefix + "_prsubjet1ptoverjetpt",1);
        self.InputTree_.SetBranchStatus(prefix + "_prsubjet2ptoverjetpt",1);
        self.InputTree_.SetBranchStatus(prefix + "_prsubjet1subjet2_deltaR",1);
        self.InputTree_.SetBranchStatus(prefix + "_planarflow04",1);
        self.InputTree_.SetBranchStatus(prefix + "_planarflow05",1);
        self.InputTree_.SetBranchStatus(prefix + "_planarflow06",1);
        self.InputTree_.SetBranchStatus(prefix + "_planarflow07",1);

        self.InputTree_.SetBranchStatus(prefix + "_deltaR_lca8jet",1);
        self.InputTree_.SetBranchStatus(prefix + "_deltaphi_METca8jet_type2",1);
        self.InputTree_.SetBranchStatus(prefix + "_deltaphi_Vca8jet_type2",1);

        if self.SignalMass_ > 0: self.InputTree_.SetBranchStatus("W_H_mass_gen",1)
        self.InputTree_.SetBranchStatus("avecomplexpolewtggH600",1)
        self.InputTree_.SetBranchStatus("avecomplexpolewtggH700",1)
        self.InputTree_.SetBranchStatus("avecomplexpolewtggH800",1)
        self.InputTree_.SetBranchStatus("avecomplexpolewtggH900",1)
        self.InputTree_.SetBranchStatus("avecomplexpolewtggH1000",1)
        self.InputTree_.SetBranchStatus("complexpolewtggH600",1)
        self.InputTree_.SetBranchStatus("complexpolewtggH700",1)
        self.InputTree_.SetBranchStatus("complexpolewtggH800",1)
        self.InputTree_.SetBranchStatus("complexpolewtggH900",1)
        self.InputTree_.SetBranchStatus("complexpolewtggH1000",1)
        self.InputTree_.SetBranchStatus("interferencewtggH600",1)
        self.InputTree_.SetBranchStatus("interferencewtggH700",1)
        self.InputTree_.SetBranchStatus("interferencewtggH800",1)
        self.InputTree_.SetBranchStatus("interferencewtggH900",1)
        self.InputTree_.SetBranchStatus("interferencewtggH1000",1)

        # GEN INFO  

        self.InputTree_.SetBranchStatus("W_tb_px",1)
        self.InputTree_.SetBranchStatus("W_tb_py",1)
        self.InputTree_.SetBranchStatus("W_tb_pz",1)
        self.InputTree_.SetBranchStatus("W_tb_E",1)
        self.InputTree_.SetBranchStatus("W_tb_Id",1)
        self.InputTree_.SetBranchStatus("W_tbbar_px",1)
        self.InputTree_.SetBranchStatus("W_tbbar_py",1)
        self.InputTree_.SetBranchStatus("W_tbbar_pz",1)
        self.InputTree_.SetBranchStatus("W_tbbar_E",1)
        self.InputTree_.SetBranchStatus("W_tbbar_Id",1)
        self.InputTree_.SetBranchStatus("W_tMet_px",1)
        self.InputTree_.SetBranchStatus("W_tMet_py",1)
        self.InputTree_.SetBranchStatus("W_tMet_pz",1)
        self.InputTree_.SetBranchStatus("W_tMet_E",1)
        self.InputTree_.SetBranchStatus("W_tMet_Id",1)
        self.InputTree_.SetBranchStatus("W_tLepton_px",1)
        self.InputTree_.SetBranchStatus("W_tLepton_py",1)
        self.InputTree_.SetBranchStatus("W_tLepton_pz",1)
        self.InputTree_.SetBranchStatus("W_tLepton_E",1)
        self.InputTree_.SetBranchStatus("W_tLepton_Id",1)

        self.InputTree_.SetBranchStatus("W_tParton_*",1)

        self.InputTree_.SetBranchStatus("W_neutrino_px_gen",1)
        self.InputTree_.SetBranchStatus("W_neutrino_py_gen",1)
        self.InputTree_.SetBranchStatus("W_neutrino_pz_gen",1)
        self.InputTree_.SetBranchStatus("W_neutrino_e_gen",1)
        if self.Channel_ == "el":
            self.InputTree_.SetBranchStatus("W_electron_px_gen",1)
            self.InputTree_.SetBranchStatus("W_electron_py_gen",1)
            self.InputTree_.SetBranchStatus("W_electron_pz_gen",1)
            self.InputTree_.SetBranchStatus("W_electron_e_gen",1)
        if self.Channel_ == "mu":
            self.InputTree_.SetBranchStatus("W_muon_px_gen",1)
            self.InputTree_.SetBranchStatus("W_muon_py_gen",1)
            self.InputTree_.SetBranchStatus("W_muon_pz_gen",1)
            self.InputTree_.SetBranchStatus("W_muon_e_gen",1)
