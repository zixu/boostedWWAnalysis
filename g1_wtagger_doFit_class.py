#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser


from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG, RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite


############################################
# Job steering #
############################################

def foo_callback(option, opt, value, parser):
          setattr(parser.values, option.dest, value.split(','))


parser = OptionParser()

parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
parser.add_option('--cprime', action="store",type="int",dest="cprime",default=10)
parser.add_option('--BRnew', action="store",type="int",dest="BRnew",default=0)

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")


parser.add_option('--fitwtagger', action='store_true', dest='fitwtagger', default=True, help='fit wtagger jet in ttbar control sample')
parser.add_option('--fitwtaggersim', action='store_true', dest='fitwtaggersim', default=False, help='fit wtagger jet in ttbar control sample with mu and el samples simultaneously')

parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")

parser.add_option('--herwig', action="store",type="int",dest="herwig",default=0)
parser.add_option('--pTbin', action="callback",callback=foo_callback,type="string",dest="pTbin",default="")


(options, args) = parser.parse_args()


ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

###############################
## doFit Class Implemetation ##
###############################

class doFit_wj_and_wlvj:

    ## contructor taking channel, signal name, range in mj, label and a workspace
    def __init__(self, in_channel,in_signal_sample, in_mj_min=30, in_mj_max=140, label="", input_workspace=None):

        print " "
        print "#####################################################"
        print "## Constructor of the fit object doFit_wj_and_wlvj ##"
        print "#####################################################"
        print " "

        #set plots style
        self.setTDRStyle();

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        ### set the channel type --> electron or muon
        self.channel=in_channel;

        ### Set the mj binning for plots
        self.BinWidth_mj=5.;
        
        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.narrow_factor=1.;

        ## set the range max properly in order to have a integer number of bins
        self.BinWidth_mj = self.BinWidth_mj/self.narrow_factor;
        nbins_mj = int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max = in_mj_min+nbins_mj*self.BinWidth_mj;

        ## declare the RooRealVar + binning
        rrv_mass_j = RooRealVar("rrv_mass_j","pruned m_{J}",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV/c^{2}");
        rrv_mass_j.setBins(nbins_mj);

        ## create the workspace and import the variable
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit"+label+"_","Workspace4fit"+label+"_");
        else:
            self.workspace4fit_ = input_workspace;
        getattr(self.workspace4fit_,"import")(rrv_mass_j);

        ## Region definition --> signal region between 65 and 105 GeV
        self.mj_sideband_lo_min=in_mj_min;
        self.mj_sideband_lo_max=65;
        self.mj_signal_min=65;
        self.mj_signal_max=105;
        self.mj_sideband_hi_min=105;
        self.mj_sideband_hi_max=in_mj_max;

        ## define ranges on mj
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("controlsample_fitting_range",40,130);

        ## directory where are the trees to be run
        self.file_Directory="AnaSigTree_forTT_new/";
             

        ## taking the root file for data and mc
        self.file_data=("treeEDBR_data_xww.root");
        self.signal_sample=in_signal_sample;
        self.file_signal=("treeEDBR_%s_xww.root"%(self.signal_sample));
        
      
        self.file_pseudodata=("treeEDBR_allBkg_xww.root");

        #WJets0 is the default PS model, WJets1 is the alternative PS model
        self.file_WJets0_mc = ("treeEDBR_WJetsPt180_xww.root");
        self.file_VV_mc     = ("treeEDBR_VV_xww.root");# WW+WZ
        self.file_TTbar_mc  = ("treeEDBR_TTBARpowheg_xww.root");
        self.file_STop_mc   = ("treeEDBR_SingleTop_xww.root");#single Top

                                                                                                                                  
        ## Define the workin poit on the N-subjettines cut

        self.wtagger_label=options.category;

        if self.wtagger_label=="HP" :
            if self.channel=="el":self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
            if self.channel=="mu":self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;

        if self.wtagger_label=="HP" :
         if self.channel=="el":self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
         if self.channel=="mu":self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
         if self.wtagger_label=="LP": self.wtagger_cut=0.75 ; self.wtagger_cut_min=0.5 ;
         if self.wtagger_label=="nocut": self.wtagger_cut=10000;

        self.categoryID=-1;
        if self.wtagger_label=="LP" and self.channel=="el": self.categoryID=0;
        if self.wtagger_label=="HP" and self.channel=="el": self.categoryID=1;
        if self.wtagger_label=="LP" and self.channel=="mu": self.categoryID=2;
        if self.wtagger_label=="HP" and self.channel=="mu": self.categoryID=3;

               
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : 2,
            'VV' : 4,
            'STop' : 7,
            'TTbar' : 210,
            'TTbar_herwig' :53,
            'ggH' : 1,
            'vbfH' : 12,
            'Signal': 1,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        ## Some basic cut vaule

        self.vpt_cut = 200; ## hadronic and leptonic W cut
        self.mass_lvj_max = 3000.; ## invariant mass of 3 body max
        self.mass_lvj_min = 200.; ## invariant mass of 3 body min
        self.pfMET_cut = 50; ## missing transverse energy
        self.lpt_cut = 50; ## lepton pT
        self.deltaPhi_METj_cut = 2.0;

        ## binning in the W jet pT for differential SF study in bins of pT

        if options.pTbin != "" :
         self.ca8_ungroomed_pt_min = int(options.pTbin[0]);
         self.ca8_ungroomed_pt_max = int(options.pTbin[1]);
        else:
         self.ca8_ungroomed_pt_min = 200;
         self.ca8_ungroomed_pt_max = 10000;

        
        ## tighter cut for the electron channel
        if self.channel=="el":
            self.pfMET_cut= 80; self.lpt_cut = 90;

        ## out txt file with info about fit and event couting
        self.file_ttbar_control_txt = "ttbar_control_%s_%s_wtaggercut%s_g1.txt"%(self.signal_sample,self.channel,self.wtagger_label);
        self.file_out_ttbar_control=open(self.file_ttbar_control_txt,"w");



    ## Set basic TDR style for canvas, pad ..etc ..
    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");
        #For the canvas:
        self.tdrStyle.SetCanvasBorderMode(0);
        self.tdrStyle.SetCanvasColor(kWhite);
        self.tdrStyle.SetCanvasDefH(600); #Height of canvas
        self.tdrStyle.SetCanvasDefW(600); #Width of canvas
        self.tdrStyle.SetCanvasDefX(0); #POsition on screen
        self.tdrStyle.SetCanvasDefY(0);
      
        #For the Pad:
        self.tdrStyle.SetPadBorderMode(0);
        self.tdrStyle.SetPadColor(kWhite);
        self.tdrStyle.SetPadGridX(False);
        self.tdrStyle.SetPadGridY(False);
        self.tdrStyle.SetGridColor(0);
        self.tdrStyle.SetGridStyle(3);
        self.tdrStyle.SetGridWidth(1);
      
        #For the frame:
        self.tdrStyle.SetFrameBorderMode(0);
        self.tdrStyle.SetFrameBorderSize(1);
        self.tdrStyle.SetFrameFillColor(0);
        self.tdrStyle.SetFrameFillStyle(0);
        self.tdrStyle.SetFrameLineColor(1);
        self.tdrStyle.SetFrameLineStyle(1);
        self.tdrStyle.SetFrameLineWidth(1);
      
        #For the histo:
        self.tdrStyle.SetHistLineColor(1);
        self.tdrStyle.SetHistLineStyle(0);
        self.tdrStyle.SetHistLineWidth(1);
        self.tdrStyle.SetEndErrorSize(2);
        self.tdrStyle.SetErrorX(0.);
        self.tdrStyle.SetMarkerStyle(20);
      
        #For the fit/function:
        self.tdrStyle.SetOptFit(1);
        self.tdrStyle.SetFitFormat("5.4g");
        self.tdrStyle.SetFuncColor(2);
        self.tdrStyle.SetFuncStyle(1);
        self.tdrStyle.SetFuncWidth(1);
      
        #For the date:
        self.tdrStyle.SetOptDate(0);
      
        #For the statistics box:
        self.tdrStyle.SetOptFile(0);
        self.tdrStyle.SetOptStat(0); #To display the mean and RMS:
        self.tdrStyle.SetStatColor(kWhite);
        self.tdrStyle.SetStatFont(42);
        self.tdrStyle.SetStatFontSize(0.025);
        self.tdrStyle.SetStatTextColor(1);
        self.tdrStyle.SetStatFormat("6.4g");
        self.tdrStyle.SetStatBorderSize(1);
        self.tdrStyle.SetStatH(0.1);
        self.tdrStyle.SetStatW(0.15);
      
        #Margins:
        self.tdrStyle.SetPadTopMargin(0.05);
        self.tdrStyle.SetPadBottomMargin(0.13);
        self.tdrStyle.SetPadLeftMargin(0.18);
        self.tdrStyle.SetPadRightMargin(0.06);
      
        #For the Global title:
        self.tdrStyle.SetOptTitle(0);
        self.tdrStyle.SetTitleFont(42);
        self.tdrStyle.SetTitleColor(1);
        self.tdrStyle.SetTitleTextColor(1);
        self.tdrStyle.SetTitleFillColor(10);
        self.tdrStyle.SetTitleFontSize(0.05);
      
        #For the axis titles:
        self.tdrStyle.SetTitleColor(1, "XYZ");
        self.tdrStyle.SetTitleFont(42, "XYZ");
        self.tdrStyle.SetTitleSize(0.03, "XYZ");
        self.tdrStyle.SetTitleXOffset(0.9);
        self.tdrStyle.SetTitleYOffset(1.5);
      
        #For the axis labels:
        self.tdrStyle.SetLabelColor(1, "XYZ");
        self.tdrStyle.SetLabelFont(42, "XYZ");
        self.tdrStyle.SetLabelOffset(0.007, "XYZ");
        self.tdrStyle.SetLabelSize(0.03, "XYZ");
      
        #For the axis:
        self.tdrStyle.SetAxisColor(1, "XYZ");
        self.tdrStyle.SetStripDecimals(kTRUE);
        self.tdrStyle.SetTickLength(0.03, "XYZ");
        self.tdrStyle.SetNdivisions(510, "XYZ");
        self.tdrStyle.SetPadTickX(1); #To get tick marks on the opposite side of the frame
        self.tdrStyle.SetPadTickY(1);
      
        #Change for log plots:
        self.tdrStyle.SetOptLogx(0);
        self.tdrStyle.SetOptLogy(0);
        self.tdrStyle.SetOptLogz(0);
      
        #Postscript options:
        self.tdrStyle.SetPaperSize(20.,20.);
        self.tdrStyle.cd();

   #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters
        
    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[]):

        if TString(mass_spectrum).Contains("_mj"): rrv_x = self.workspace4fit_.var("rrv_mass_j"); ## inport the mass spectrum variable as rrv_x

        ## Voigtian Pdf
        if in_model_name == "Voig":
              
            rrv_mean_voig = RooRealVar("rrv_mean_voig"+label+"_"+self.channel,"rrv_mean_voig"+label+"_"+self.channel,84,78,88);# W mass: 80.385
            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.channel,"rrv_width_voig"+label+"_"+self.channel,7.,1,40);
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel,"rrv_sigma_voig"+label+"_"+self.channel,5,0.01,20);

            model_pdf = RooVoigtian("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

        ## Breit Wigner Pdf
        if in_model_name == "BW":

            rrv_mean_BW = RooRealVar("rrv_mean_BW"+label+"_"+self.channel,"rrv_mean_BW"+label+"_"+self.channel,84,78, 88);
            rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.channel,"rrv_width_BW"+label+"_"+self.channel,20,1,40);

            model_pdf = RooBreitWigner("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_BW,rrv_width_BW);

        ## Double Voigtian W + Z mass peak
        if in_model_name == "2Voig":
            rrv_mean_voig = RooRealVar("rrv_mean_voig"+label+"_"+self.channel,"rrv_mean_voig"+label+"_"+self.channel,84,78,88);#W mass 80.385
            rrv_shift_2Voig = RooRealVar("rrv_shift_2Voig"+label+"_"+self.channel,"rrv_shift_2Voig"+label+"_"+self.channel,10.8026) # Z mass: 91.1876; shift=91.1876-80.385=10.8026
            rrv_mean_shifted = RooFormulaVar("rrv_mean_voig2"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean_voig,rrv_shift_2Voig));

            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.channel,"rrv_width_voig"+label+"_"+self.channel,16.,6,26);
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel,"rrv_sigma_voig"+label+"_"+self.channel,0.);

            rrv_frac=RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,0.5,0.,1.);
            
            model_voig1 = RooVoigtian("model_voig1"+label+"_"+self.channel+mass_spectrum,"model_voig1"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);
            model_voig2 = RooVoigtian("model_voig2"+label+"_"+self.channel+mass_spectrum,"model_voig2"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_shifted,rrv_width_voig,rrv_sigma_voig);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, RooArgList(model_voig1,model_voig2), RooArgList(rrv_frac));

            ## Gaussian function
        if in_model_name == "Gaus":

           rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,84,78,92);
           rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,7,1,15);

           model_pdf = RooGaussian("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_mean_gaus,rrv_sigma_gaus);
                                                    

        ## Crystal Ball function
        if in_model_name == "CB":

            rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,84,78,88);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,7,4,15);
            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,-2,-4,-0.5);
            rrv_n_CB = RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel,2,0.,4);
            model_pdf = RooCBShape("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

        ## Convolution BW with a CB
        if in_model_name == "CBBW": # FFT: BreitWigner*CBShape

            rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,84.0,78,88);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,7,4,15);
            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,-2,-4,-1);
            rrv_n_CB = RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel,0.5,0.,2);

            cbshape = RooCBShape("cbshape"+label+"_"+self.channel,"cbshape"+label+"_"+self.channel, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

            rrv_mean_BW = RooRealVar("rrv_mean_BW"+label+"_"+self.channel,"rrv_mean_BW"+label+"_"+self.channel,0);
            rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.channel,"rrv_width_BW"+label+"_"+self.channel,10,5,20);

            bwshape = RooBreitWigner("bw"+label+"_"+self.channel,"bw"+label+"_"+self.channel, rrv_x,rrv_mean_BW,rrv_width_BW);
            
            model_pdf = RooFFTConvPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, cbshape, bwshape);

        ## Convolution between landau and gaussian
        if in_model_name == "LDGaus": # FFT: Landau*Gaus

            rrv_mean_landau = RooRealVar("rrv_mean_landau"+label+"_"+self.channel,"rrv_mean_landau"+label+"_"+self.channel,84.0,78,88);
            rrv_sigma_landau = RooRealVar("rrv_sigma_landau"+label+"_"+self.channel,"rrv_sigma_landau"+label+"_"+self.channel,7,4,15);

            landau = RooLandau("landau"+label+"_"+self.channel,"landau"+label+"_"+self.channel, rrv_x,rrv_mean_landau,rrv_sigma_landau);

            rrv_mean_gaus = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,0);
            rrv_sigma_gaus=RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,16,10,20);

            gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            model_pdf = RooFFTConvPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, landau, gaus);


        ## Product between Erf and Expnential
        if in_model_name == "ErfExp" :

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.1, -1.e-4);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,60.,30.,120);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10.,60.);
            model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        if in_model_name == "ErfExp_2" :

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2,-1e-4);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,85.,25.,98);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,45.,20, 70.);
            model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);


        ## ErfExp + Gaussian
        if in_model_name == "ErfExpGaus":
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.4,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,100.,10.,300.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,100.);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,82,78,87);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,7,4,10);
            gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.7,0.,1.);
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        if in_model_name == "ErfExpGaus_sp":#offset == mean
            rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2,0.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,200.);
            rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,78,88);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp);

            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,10);
            gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.5,0.,1.);
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ## Exp +Gaus
        if in_model_name == "ExpGaus":
            rrv_c_Exp = RooRealVar("rrv_c_ExpGaus"+label+"_"+self.channel,"rrv_c_ExpGaus"+label+"_"+self.channel,0.05,-0.5,0.5);
            exp = ROOT.RooExponential("ExpGaus"+label+"_"+self.channel,"ExpGaus"+label+"_"+self.channel,rrv_x,rrv_c_Exp);

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,78,92);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,16);
            gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.5,0.,1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))


        ## Erf*Exp + Gaus
        if in_model_name == "ErfExpGaus_v0":
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,100.,10.,140.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,100.);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean_gaus = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,84,78,88);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,7,4,10);
            gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.7,0.,1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))
    

        if in_model_name == "ErfExpGaus_v2":
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-10.,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,100.,10.,140.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,100.);
            rrv_mean_gaus = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,84,78,88);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,7,4,10);
            rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,200.,0.,1000.);

            model_pdf = ROOT.RooErfExp_Gaus_Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_mean_gaus,rrv_sigma_gaus,rrv_high );

        ## Erf*exp + 2Gaus
        if in_model_name == "ErfExp2Gaus":
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.2,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,100.,10.,140.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,100.);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean1_gaus=RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,78,88);
            rrv_mean2_gaus=RooRealVar("rrv_mean2_gaus"+label+"_"+self.channel,"rrv_mean2_gaus"+label+"_"+self.channel,180,170,190);
            rrv_sigma1_gaus=RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,10);
            rrv_sigma2_gaus=RooRealVar("rrv_sigma2_gaus"+label+"_"+self.channel,"rrv_sigma2_gaus"+label+"_"+self.channel,10,7,15);
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_high1 = RooRealVar("rrv_high1"+label+"_"+self.channel,"rrv_high1"+label+"_"+self.channel,0.6,0.,1.);
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.channel,"rrv_high2"+label+"_"+self.channel,0.4,0.,1.);

            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus1,gaus2),RooArgList(rrv_high1,rrv_high2))

        if in_model_name == "2Gaus": ## double gaussian for the resonant of mJ

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,83.5, 80., 86.5);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7.6,3.5,10.);
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel, 4.5, 0.,9.);
            rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,3.,0.,5.5);
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel, 0.5, 0., 1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

        if in_model_name == "GausChebychev":

            rrv_mean_gaus = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,83.5, 80., 86.5);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,7.6,3.5,30.);
            gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);
                
            rrv_p0_cheb = RooRealVar("rrv_p0_cheb"+label+"_"+self.channel,"rrv_p0_cheb"+label+"_"+self.channel, -0.3, -30., 30.)
            rrv_p1_cheb = RooRealVar("rrv_p1_cheb"+label+"_"+self.channel,"rrv_p1_cheb"+label+"_"+self.channel, -0.1, -30., 30.)
            cheb = RooChebychev("cheb"+label+"_"+self.channel,"cheb"+label+"_"+self.channel, rrv_x, RooArgList(rrv_p0_cheb, rrv_p1_cheb) );

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,0.5, 0., 1.)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus,cheb),RooArgList(rrv_frac),1)

        if in_model_name == "Exp":

          rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel, -0.05, -0.2, 0.);
          model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp);
                                                                  

        if in_model_name == "2_2Gaus":# 2Gaus + 2Gaus for the resonant W and Z part of the mJ spectra
              
            mean1_tmp =8.3145e+01; mean1_tmp_err =1.63e-01;
            deltamean_tmp =6.6321e+00; deltamean_tmp_err =1.21e+00;
            sigma1_tmp =7.5097e+00; sigma1_tmp_err =2.01e-01;
            scalesigma_tmp=3.8707e+00; scalesigma_tmp_err=2.20e-01;
            frac_tmp =6.4728e-01; frac_tmp_err =2.03e-02;


            if self.wtagger_cut==0.50:
                mean1_tmp =8.3141e+01; mean1_tmp_err =1.63e-01;
                deltamean_tmp =6.9129e+00; deltamean_tmp_err =1.24e+00;
                sigma1_tmp =7.5145e+00; sigma1_tmp_err =1.99e-01;
                scalesigma_tmp=3.6819e+00; scalesigma_tmp_err=2.11e-01;
                frac_tmp =6.7125e-01; frac_tmp_err =2.09e-02;

            rrv_shift = RooRealVar("rrv_shift"+label+"_"+self.channel,"rrv_shift"+label+"_"+self.channel,10.8026) # Z mass: 91.1876; shift=91.1876-80.385=10.8026

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,0.,-8,10);
            rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac1 = RooRealVar("rrv_frac1"+label+"_"+self.channel,"rrv_frac1"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            gausguas_1 = RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum,"gausguas_1"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)

            rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_shift));
            rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean2_gaus, rrv_shift));
            gaus3 = RooGaussian("gaus3"+label+"_"+self.channel,"gaus3"+label+"_"+self.channel, rrv_x,rrv_mean3_gaus,rrv_sigma1_gaus);
            gaus4 = RooGaussian("gaus4"+label+"_"+self.channel,"gaus4"+label+"_"+self.channel, rrv_x,rrv_mean4_gaus,rrv_sigma2_gaus);

            gausguas_2 = RooAddPdf("gausguas_2"+label+"_"+self.channel+mass_spectrum,"gausguas_2"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus3,gaus4),RooArgList(rrv_frac1),1)

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,0.74)#,0.5,1.0);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gausguas_1,gausguas_2),RooArgList(rrv_frac),1)

        ## ErfExp + 2Gaus
        if in_model_name == "2Gaus_ErfExp":

            mean1_tmp =8.3145e+01; mean1_tmp_err =1.63e-01;
            deltamean_tmp =6.6321e+00; deltamean_tmp_err =1.21e+00;
            sigma1_tmp =7.5097e+00; sigma1_tmp_err =2.01e-01;
            scalesigma_tmp=3.8707e+00; scalesigma_tmp_err=2.20e-01;
            frac_tmp =6.4728e-01; frac_tmp_err =2.03e-02;


            if self.wtagger_cut==0.50:
                mean1_tmp =8.3141e+01; mean1_tmp_err =1.63e-01;
                deltamean_tmp =6.9129e+00; deltamean_tmp_err =1.24e+00;
                sigma1_tmp =7.5145e+00; sigma1_tmp_err =1.99e-01;
                scalesigma_tmp=3.6819e+00; scalesigma_tmp_err=2.11e-01;
                frac_tmp =6.7125e-01; frac_tmp_err =2.09e-02;

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp)#, deltamean_tmp, deltamean_tmp);
            rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp);
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel,"rrv_frac_2gaus"+label+"_"+self.channel,frac_tmp);#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            c0_tmp = -2.8628e-02 ; c0_tmp_err = 6.08e-03;
            offset_tmp= 7.6259e+01 ; offset_tmp_err = 9.17e+00;
            width_tmp = 3.4207e+01 ; width_tmp_err = 3.18e+00;

            if self.wtagger_cut==0.50:
                c0_tmp = -2.9893e-02 ; c0_tmp_err = 6.83e-03;
                offset_tmp= 7.9350e+01 ; offset_tmp_err = 9.35e+00;
                width_tmp = 3.3083e+01 ; width_tmp_err = 2.97e+00;

            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
            rrv_offset_ErfExp= RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp)#, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-10, width_tmp+10);
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum,"erfexp"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel, 0.5,0.,1.);
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)

        ## 2Gaus per ttbar fit
        if in_model_name == "2Gaus_ttbar":

            mean1_tmp = 8.3141e+01; mean1_tmp_err = 1.63e-01;
            deltamean_tmp = 6.9129e+00; deltamean_tmp_err = 1.24e+00;
            sigma1_tmp = 7.5145e+00; sigma1_tmp_err = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp = 6.7125e-01; frac_tmp_err = 2.09e-02;

            if TString(label).Contains("herwig") and not TString(label).Contains("data"):

             mean1_tmp = 8.3141e+01; mean1_tmp_err = 1.63e-01;
             deltamean_tmp = 6.9129e+00; deltamean_tmp_err = 1.24e+00;
             sigma1_tmp = 7.5145e+00; sigma1_tmp_err = 1.99e-01;
             scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err= 2.11e-01;
             frac_tmp = 6.7125e-01; frac_tmp_err = 2.09e-02;

            ## constrain the peak and the width among electron and muon channel to be the same in case of simultaneous fit
            if self.channel=="el":
                if self.workspace4fit_.var("rrv_mean1_gaus%s_mu"%(label)) and self.workspace4fit_.var("rrv_sigma1_gaus%s_mu"%(label)):
                    rrv_mean1_gaus = self.workspace4fit_.var("rrv_mean1_gaus%s_mu"%(label));
                    rrv_sigma1_gaus = self.workspace4fit_.var("rrv_sigma1_gaus%s_mu"%(label));
                else:
                    rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
                    rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4);
            if self.channel=="mu":
                if self.workspace4fit_.var("rrv_mean1_gaus%s_el"%(label)) and self.workspace4fit_.var("rrv_sigma1_gaus%s_el"%(label)):
                    rrv_mean1_gaus = self.workspace4fit_.var("rrv_mean1_gaus%s_el"%(label));
                    rrv_sigma1_gaus = self.workspace4fit_.var("rrv_sigma1_gaus%s_el"%(label));
                else:
                    rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
                    rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );

            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp)#,deltamean_tmp-deltamean_tmp_err*5, deltamean_tmp+deltamean_tmp_err*5);
            rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp)#,scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp)#,frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

        ## Gaus + pol 2 background for resonant part in the fail sample
        if in_model_name == "GausChebychev_ttbar_failtau2tau1cut":

          p0_tmp = -3.5459e-01; p0_tmp_err = 5.04e-02;
          p1_tmp = -1.2790e-01; p1_tmp_err = 6.74e-02;
          frac_tmp = 2.7324e-01; frac_tmp_err = 2.48e-02;

          if TString(label).Contains("herwig") and not TString(label).Contains("data") and self.channel == "mu":

           p0_tmp = -6.6824e-01; p0_tmp_err = 5.04e-02;
           p1_tmp = -5.3365e-02; p1_tmp_err = 6.74e-02;
           frac_tmp = 2.0421e-01; frac_tmp_err = 2.48e-02;

          ## take the same gaussian used in the pass sample
          if TString(label).Contains("data") and not TString(label).Contains("herwig"):
                gaus = self.workspace4fit_.pdf("gaus1_ttbar_data_"+self.channel);
          if TString(label).Contains("TotalMC") and not TString(label).Contains("herwig"):
                gaus = self.workspace4fit_.pdf("gaus1_ttbar_TotalMC_"+self.channel);

          if TString(label).Contains("data") and TString(label).Contains("herwig"):
                gaus = self.workspace4fit_.pdf("gaus1_ttbar_data_herwig_"+self.channel);
          if TString(label).Contains("TotalMC") and TString(label).Contains("herwig"):
                gaus = self.workspace4fit_.pdf("gaus1_ttbar_TotalMC_herwig_"+self.channel);
                
          rrv_p0_cheb = RooRealVar("rrv_p0_cheb"+label+"_"+self.channel,"rrv_p0_cheb"+label+"_"+self.channel,p0_tmp)#,p0_tmp-p0_tmp_err*4,p0_tmp+p0_tmp_err*4);
          rrv_p1_cheb = RooRealVar("rrv_p1_cheb"+label+"_"+self.channel,"rrv_p1_cheb"+label+"_"+self.channel,p1_tmp)#,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
          cheb = RooChebychev("cheb"+label+"_"+self.channel,"cheb"+label+"_"+self.channel, rrv_x, RooArgList(rrv_p0_cheb, rrv_p1_cheb) );

          rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp)#,frac_tmp-frac_tmp_err*4,frac_tmp+frac_tmp_err*4);
          model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus,cheb),RooArgList(rrv_frac),1)

        ## double gauss for fail sample
        if in_model_name == "2Gaus_ttbar_failtau2tau1cut":

            mean1_tmp = 8.3209e+01; mean1_tmp_err = 1.17e-01;
            deltamean_tmp = 1.1427e+00; deltamean_tmp_err = 1.03e+00;
            sigma1_tmp = 7.4932e+00; sigma1_tmp_err = 1.44e-01;
            scalesigma_tmp = 4.5922e+00; scalesigma_tmp_err = 2.87e-01;
            frac_tmp = 5.7910e-01; frac_tmp_err = 1.59e-02;

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp, deltamean_tmp-deltamean_tmp_err*4, deltamean_tmp+deltamean_tmp_err*4);
            rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

        ## Erf*Exp in the ttbar sample
        if in_model_name == "ErfExp_ttbar":

            c0_tmp = -2.9893e-02 ; c0_tmp_err = 6.83e-03;
            offset_tmp = 7.9350e+01 ; offset_tmp_err = 9.35e+00;
            width_tmp = 3.3083e+01 ; width_tmp_err = 2.97e+00;

            if TString(label).Contains("herwig") and not TString(label).Contains("data"):

              c0_tmp = -2.9357e-02 ; c0_tmp_err = 6.83e-03;
              offset_tmp = 7.9350e+01 ; offset_tmp_err = 9.35e+00;
              width_tmp = 3.3216e+01 ; width_tmp_err = 2.97e+00;
                                                                      
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp,width_tmp-10, width_tmp+10);
                
            model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            ## Constrint the parameter within one sigma on the background part of the non resonant component in the pass sample
            self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);
            self.addConstraint(rrv_offset_ErfExp,offset_tmp, offset_tmp_err,ConstraintsList);
            self.addConstraint(rrv_width_ErfExp,width_tmp, width_tmp_err,ConstraintsList);

        ## Erf*Exp for failing events
        if in_model_name == "ErfExp_ttbar_failtau2tau1cut":

            c0_tmp = -1.0143e-01 ; c0_tmp_err = 1.46e-02;
            offset_tmp = 2.7718e+02 ; offset_tmp_err = 4.92e+01;
            width_tmp = 7.1891e+01 ; width_tmp_err = 4.69e+00;

            if TString(label).Contains("herwig") and not TString(label).Contains("data"):

              c0_tmp = -1.0141e-01 ; c0_tmp_err = 1.46e-02;
              offset_tmp = 2.6730e+02 ; offset_tmp_err = 4.92e+01;
              width_tmp = 7.1505e+01 ; width_tmp_err = 4.69e+00;

                                                                      
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,offset_tmp)#offset_tmp-offset_tmp_err*5,offset_tmp+offset_tmp_err*5);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp)#width_tmp-width_tmp_err*5, width_tmp+width_tmp_err*5);
            model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);
            ## constrain just the slope of the exponential
            self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);

        ## extreme fail cut
        if in_model_name == "Exp_ttbar_extremefailtau2tau1cut":

            c0_tmp = -3.0278e-02 ; c0_tmp_err = 5.16e-03;

            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp);
            ## constraint within one sigma
            self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);

        if in_model_name == "Exp_bkg_extremefailtau2tau1cut":
           c0_tmp = -4.2105e-02 ; c0_tmp_err = 2.61e-03;

           rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
           model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp);
           ## constraint within one sigma
           self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);

        ## Erf*Exp + 2Gaus
        if in_model_name == "ErfExp2Gaus_ttbar":
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel, -9.72533e-02, -9.72533e-02-2*6.05691e-03, -9.72533e-02+2*6.05691e-03 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,2.21236e+02, 2.21236e+02-2*9.53939e+00, 2.21236e+02+2*9.53939e+00 );
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,5.79251e+01,5.79251e+01-2*1.22221e+00,5.79251e+01+2*1.22221e+00);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean1_gaus=RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,8.36227e+01,8.36227e+01-2*1.60148e-01,8.36227e+01+2*1.60148e-01);
            rrv_mean2_gaus=RooRealVar("rrv_mean2_gaus"+label+"_"+self.channel,"rrv_mean2_gaus"+label+"_"+self.channel,8.45644e+01 ,8.45644e+01-2*9.21003e-01,8.45644e+01+2*9.21003e-01 );
            rrv_sigma1_gaus=RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7.45239e+00,7.45239e+00-2*2.05052e-01,7.45239e+00+2*2.05052e-01);
            rrv_sigma2_gaus=RooRealVar("rrv_sigma2_gaus"+label+"_"+self.channel,"rrv_sigma2_gaus"+label+"_"+self.channel,3.39773e+01,3.39773e+01-2*2.44179e+00,3.39773e+01+2*2.44179e+00);
            rrv_high1 = RooRealVar("rrv_high1"+label+"_"+self.channel,"rrv_high1"+label+"_"+self.channel,0.5,0.,1.);
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.channel,"rrv_high2"+label+"_"+self.channel,0.5622,0.51,0.61);
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus1,gaus2),RooArgList(rrv_high1,rrv_high2),1)

        ## Erf*Exp + 2Gaus fail sample
        if in_model_name == "ErfExp2Gaus_ttbar_failtau2tau1cut":
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-2.81053e-02 ,-2.81053e-02-2*5.82712e-03 ,-2.81053e-02+2*5.82712e-03 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,7.63158e+01,7.63158e+01-2*8.66322e+00,7.63158e+01+2*8.66322e+00);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,3.38156e+01,3.38156e+01-2*3.02392e+00,3.38156e+01+2*3.02392e+00);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean1_gaus=RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,8.34758e+01,8.34758e+01-2*1.62468e-01,8.34758e+01+2*1.62468e-01);
            rrv_mean2_gaus=RooRealVar("rrv_mean2_gaus"+label+"_"+self.channel,"rrv_mean2_gaus"+label+"_"+self.channel,8.99116e+01,8.99116e+01-2*7.47952e-01,8.99116e+01+2*7.47952e-01);
            rrv_sigma1_gaus=RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7.45827e+00,7.45827e+00-2*1.95598e-01,7.45827e+00+2*1.95598e-01);
            rrv_sigma2_gaus=RooRealVar("rrv_sigma2_gaus"+label+"_"+self.channel,"rrv_sigma2_gaus"+label+"_"+self.channel,2.90745e+01,2.90745e+01-2*1.77035e+00,2.90745e+01+2*1.77035e+00);
            rrv_high1 = RooRealVar("rrv_high1"+label+"_"+self.channel,"rrv_high1"+label+"_"+self.channel,0.5,0.,1.);
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.channel,"rrv_high2"+label+"_"+self.channel,0.6323,0.57,0.7);
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus1,gaus2),RooArgList(rrv_high1,rrv_high2),1)

        ## Erf*Exp +gaus+ Voig
        if in_model_name == "ErfExpGausGaus":
            rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.1,-10.,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,90.,10.,140.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,100.);
            rrv_mean_voig = RooRealVar("rrv_mean_voig"+label+"_"+self.channel,"rrv_mean_voig"+label+"_"+self.channel,84,78,88);
            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.channel,"rrv_width_voig"+label+"_"+self.channel,7,1,20);
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.channel,"rrv_sigma_voig"+label+"_"+self.channel,5,1,100);
            rrv_high1 = RooRealVar("rrv_high1"+label+"_"+self.channel,"rrv_high1"+label+"_"+self.channel,1,0.,200.);
            rrv_mean_gaus = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,92)#,160,187);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,20)#,0.1,100);
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.channel,"rrv_high2"+label+"_"+self.channel,0.)#,0.,0.);

            model_pdf = ROOT.RooErfExp_Voig_Gaus_Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig,rrv_high1,rrv_mean_gaus,rrv_sigma_gaus,rrv_high2 );


        ## import the pdf in the workspace and save it
        getattr(self.workspace4fit_,"import")(model_pdf)
        ## retun the pdf with the right name
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)


   ## Gaussian contraint of a parameter of a pdf
    def addConstraint(self, rrv_x, x_mean, x_sigma, ConstraintsList):

            rrv_x_mean = RooRealVar(rrv_x.GetName()+"_mean",rrv_x.GetName()+"_mean",x_mean );
            rrv_x_sigma = RooRealVar(rrv_x.GetName()+"_sigma",rrv_x.GetName()+"_sigma",x_sigma );
            constrainpdf_x = RooGaussian("constrainpdf_"+rrv_x.GetName(),"constrainpdf_"+rrv_x.GetName(),rrv_x, rrv_x_mean, rrv_x_sigma);
            ## import in the workspace and save the name of constriant pdf
            getattr(self.workspace4fit_,"import")(constrainpdf_x)
            ConstraintsList.append(constrainpdf_x.GetName());

  
    ### get an mj model from the workspace givin the label
    def get_mj_Model(self,label):
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_mj_Model(self, label ):
        rdataset_General_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_General = self.get_mj_Model(label);
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_mj);
        par=parameters_General.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")) and (not (options.fitwtaggersim 
or options.fitwtagger)):
                param.Print();
            param.setConstant(kTRUE);
            param=par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    ### fix only the ttbar component using the default label --> for extended pdf
    def get_TTbar_mj_Model(self,label="_TTbar"):
        return self.get_General_mj_Model(label);

    ### fix only the stop component using the default label --> for extended pdf
    def get_STop_mj_Model(self,label="_STop"):
        return self.get_General_mj_Model(label);

    ### fix only the VV component using the default label --> for extended pdf
    def get_VV_mj_Model(self,label="_VV"):
        return self.get_General_mj_Model(label);

    ### fix only the WJets model --> for extended pdf
    def get_WJets_mj_Model(self,label):
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_WJets = self.get_mj_Model(label);
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mj);
        par=parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            if (paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or 
paraName.Contains("rrv_p1_User1_WJets") or paraName.Contains("rrv_p1_User1_WJets")) :
             param.setConstant(kTRUE);
            else:
             param.setConstant(0);
            param=par.Next()
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))


    ### fix a given model taking the label, and the region --> for extended pdf
    def fix_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj"):
        rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
        model = self.get_mlvj_Model(label,mlvj_region);
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param=par.Next()

    ### fix a pdf in a different way --> for RooAbsPdf
    def fix_Pdf(self,model_pdf,argset_notparameter):
        parameters = model_pdf.getParameters(argset_notparameter);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()

    ### print the parameters of a given pdf --> only non constant ones
    def ShowParam_Pdf(self,model_pdf,argset_notparameter):
        model_pdf.Print()
        parameters = model_pdf.getParameters(argset_notparameter);
        par = parameters.createIterator(); par.Reset();
        param = par.Next()
        while (param):
            if not param.isConstant():
                param.Print();
                if (param.getVal()-param.getMin())< param.getError()*1 or (param.getMax()- param.getVal())< param.getError()*1:
                    param.Print();
            param=par.Next()

    ### change a dataset to a histpdf roofit object
    def change_dataset_to_histpdf(self,x,dataset):
        datahist = dataset.binnedClone(dataset.GetName()+"_binnedClone",dataset.GetName()+"_binnedClone")
        histpdf = RooHistPdf(dataset.GetName()+"_histpdf",dataset.GetName()+"_histpdf",RooArgSet(x),datahist)
        getattr(self.workspace4fit_,"import")(histpdf)

    ### change from a dataset to a histogramm of Roofit
    def change_dataset_to_histogram(self, x,dataset,label=""):
        datahist=dataset.binnedClone(dataset.GetName()+"_binnedClone",dataset.GetName()+"_binnedClone")
        nbin=int( (x.getMax()-x.getMin())/self.BinWidth_mj);
        if label=="":
            return datahist.createHistogram("histo_%s"%(dataset.GetName()),x, RooFit.Binning( nbin ,x.getMin(),x.getMax()));
        else:
            return datahist.createHistogram("histo_"+label,x, RooFit.Binning( nbin,x.getMin(),x.getMax()));


 ## Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[], area_init_value=500):

      ##### define an extended pdf from a standard Roofit One
      print " "
      print "###############################################"
      print "## Make model : ",label," ",in_model_name,"##";
      print "###############################################"
      print " "

      rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum,"rrv_number"+label+"_"+self.channel+mass_spectrum,area_init_value,0.,1e7);
      ## call the make RooAbsPdf method
      model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList)
      print " model pdf : "
      model_pdf.Print();
      
      ## create the extended pdf
      model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
      print " Extended pdf : "

      #### put all the parameters ant the shape in the workspace
      getattr(self.workspace4fit_,"import")(rrv_number)
      getattr(self.workspace4fit_,"import")(model)
      self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();
      ## return the total extended pdf
      return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum);


    ### Make a generic model for the ttbar control sample fit
    def make_Model_for_ttbar_controlsample(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],info=""):

        print " "
        print "###################################################################################"
        print " ## make_Model_for_ttbar_controlsample : ",label," ",in_model_name," ",info," "##";
        print "###################################################################################"
        print " "
        
        if TString(label).Contains("_ttbar_data") and not TString(label).Contains("failtau2tau1cut"):
            rrv_number_total = RooRealVar("rrv_number_total_ttbar_data"+info+"_"+self.channel,"rrv_number_total_ttbar_data"+info+"_"+self.channel,500,0.,1e7);
            eff_ttbar = RooRealVar("eff_ttbar_data"+info+"_"+self.channel,"eff_ttbar_data"+info+"_"+self.channel,0.7,0.3,0.9);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "@0*@1", RooArgList(rrv_number_total,eff_ttbar));

        ## fail fit has to be done after the pass one
        elif TString(label).Contains("_ttbar_data") and TString(label).Contains("failtau2tau1cut"):
            ## correlated paramters among the two cateogry of pass and fail sample
            rrv_number_total = self.workspace4fit_.var("rrv_number_total_ttbar_data"+info+"_"+self.channel);
            eff_ttbar = self.workspace4fit_.var("eff_ttbar_data"+info+"_"+self.channel);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total));

        elif TString(label).Contains("_ttbar_TotalMC") and not TString(label).Contains("failtau2tau1cut"):
            rrv_number_total = RooRealVar("rrv_number_total_ttbar_TotalMC"+info+"_"+self.channel,"rrv_number_total_ttbar_TotalMC"+info+"_"+self.channel,500,0.,1e7);
            eff_ttbar = RooRealVar("eff_ttbar_TotalMC"+info+"_"+self.channel,"eff_ttbar_TotalMC"+info+"_"+self.channel,0.7,0.3,0.9);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "@0*@1", RooArgList(eff_ttbar,rrv_number_total));

        elif TString(label).Contains("_ttbar_TotalMC") and TString(label).Contains("failtau2tau1cut"):
            rrv_number_total = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+info+"_"+self.channel);
            eff_ttbar = self.workspace4fit_.var("eff_ttbar_TotalMC"+info+"_"+self.channel);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total)) ;

        rrv_number_total.Print();
        eff_ttbar.Print();
        rrv_number.Print();
        ## make the pdf for the given model
        model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList)
        print " model Pdf :"
        model_pdf.Print();
        ## make the extended pdf via rrv_number
        print " extended PDF :"
        model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
        ## import in the worspace
        getattr(self.workspace4fit_,"import")(model)
        print "model"+label+"_"+self.channel+mass_spectrum
        self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();
        ## return the model
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum);



    ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

      ## import variable and dataset
      rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
      rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.channel+"_mj");
      rdataset_mj.Print();

      ## make the extended model
      model = self.make_Model(label,in_model_name);
      ## Fit Sequence setting the errors correctly
      rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
      rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
      rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));

      ## Plot the result
      mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins())));
      rdataset_mj.plotOn(mplot, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

      ## draw the error band for an extend pdf
      draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,6,"L");
      ## re-draw the dataset
      rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
      ## draw the function
      model.plotOn(mplot);

      ## Get the pull
      mplot_pull = self.get_pull(rrv_mass_j, mplot);
      mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
      
      parameters_list = model.getParameters(rdataset_mj);
      self.draw_canvas_with_pull( mplot,mplot_pull,parameters_list,"plots_%s_%s_%s_g1/m_j_fitting%s_wtaggercut%s/"%(options.additioninformation,self.channel,self.wtagger_label,additioninformation, self.wtagger_label),label+in_file_name+"_"+str(self.ca8_ungroomed_pt_min)+"_"+str(self.ca8_ungroomed_pt_max), in_model_name)
      

      self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal()  )
      self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal()  )
      if TString(label).Contains("ggH"):
            self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()  )
            self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError())
      self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();

      #correct the mean and sigma from the ttbar contral sample
      par=parameters_list.createIterator();
      par.Reset();
      param=par.Next()
      while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")) and (not (options.fitwtaggersim or options.fitwtagger)):
                #param.Print();
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                    #param.Print(); raw_input("mean"+label);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                    #param.Print(); raw_input("mean"+label);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                    #param.Print(); raw_input("sigma"+label);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
                    #param.Print(); raw_input("sigma"+label);
            param=par.Next()
                                       

    ############# ---------------------------------------------------
    def fit_mj_TTbar_controlsample(self,in_file_name):
        self.workspace4fit_.var("rrv_number_dataset_signal_region_data_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_WJets0_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_"+self.channel+"_mj").Print()

        number_dataset_signal_region_data_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_data_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_data_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_data_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region: %s +/- sqrt(%s)"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj);
        self.file_out_ttbar_control.write("%s channel SF: \n"%(self.channel));
        self.file_out_ttbar_control.write("event number of data in signal_region: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj));
        number_dataset_signal_region_TotalMC_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_TotalMC_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_TotalMC_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_TotalMC_"+self.channel+"_mj").getVal();
        print "event number of TotalMC in signal_region: %s +/- sqrt(%s) "%(number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj);
        self.file_out_ttbar_control.write("event number of TotalMC in signal_region: %s +/- sqrt(%s) \n"%(number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj));


        number_dataset_signal_region_before_mva_data_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_before_mva_data_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_mva_error2_data_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_before_mva_error2_data_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region before_mva: %s +/- sqrt(%s)"%(number_dataset_signal_region_before_mva_data_mj, number_dataset_signal_region_before_mva_error2_data_mj);
        self.file_out_ttbar_control.write("event number of data in signal_region before_mva: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_before_mva_data_mj, number_dataset_signal_region_before_mva_error2_data_mj));

        number_dataset_signal_region_before_mva_TotalMC_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_before_mva_TotalMC_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_mva_error2_TotalMC_mj=self.workspace4fit_.var("rrv_number_dataset_signal_region_before_mva_error2_TotalMC_"+self.channel+"_mj").getVal();
        print "event number of TotalMC in signal_region before_mva: %s +/- sqrt(%s) "%(number_dataset_signal_region_before_mva_TotalMC_mj, number_dataset_signal_region_before_mva_error2_TotalMC_mj);
        self.file_out_ttbar_control.write("event number of TotalMC in signal_region before_mva: %s +/- sqrt(%s) \n"%(number_dataset_signal_region_before_mva_TotalMC_mj, number_dataset_signal_region_before_mva_error2_TotalMC_mj));

        # wtagger_eff reweight: only reweight the efficiency difference between MC and data
        wtagger_eff_MC  = number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_before_mva_TotalMC_mj;
        wtagger_eff_data= number_dataset_signal_region_data_mj/number_dataset_signal_region_before_mva_data_mj;

        wtagger_eff_reweight=wtagger_eff_data/wtagger_eff_MC;
        wtagger_eff_reweight_err=wtagger_eff_reweight*TMath.Sqrt(
                number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj +  
                number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj +  
                number_dataset_signal_region_before_mva_error2_data_mj/number_dataset_signal_region_before_mva_data_mj/number_dataset_signal_region_data_mj +  
                number_dataset_signal_region_before_mva_error2_TotalMC_mj/number_dataset_signal_region_before_mva_TotalMC_mj/number_dataset_signal_region_before_mva_TotalMC_mj 
                );
        
        print "wtagger efficiency of %s channel"%(self.channel )
        print "wtagger_eff_MC   = %s "%(wtagger_eff_MC )
        print "wtagger_eff_data = %s "%(wtagger_eff_data )
        print "wtagger_eff_reweight = %s +/- %s"%(wtagger_eff_reweight, wtagger_eff_reweight_err)
        self.file_out_ttbar_control.write("wtagger_eff_MC   = %s \n"%(wtagger_eff_MC ));
        self.file_out_ttbar_control.write("wtagger_eff_data = %s \n"%(wtagger_eff_data ));
        self.file_out_ttbar_control.write("wtagger_eff_reweight = %s +/- %s\n"%(wtagger_eff_reweight, wtagger_eff_reweight_err));

        ## wtagger reweight: reweight the event number difference in signal region after mva cut between MC and data
        #wtagger_reweight=number_dataset_signal_region_data_mj/number_dataset_signal_region_TotalMC_mj;
        #wtagger_reweight_err=wtagger_reweight*TMath.Sqrt( 
        #        number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj +  
        #        number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj 
        #        );
        #print "wtagger_reweight = %s +/- %s"%(wtagger_reweight, wtagger_reweight_err)
        #self.file_out_ttbar_control.write("wtagger_reweight = %s +/- %s\n"%(wtagger_reweight, wtagger_reweight_err));

    ############# -----------
    def fit_mj_TTbar_controlsample(self,in_file_name,label=""):

        ##### Print the final result for the number of events passing the cut and before the cut + the efficiency for the W-tagging

        self.workspace4fit_.var("rrv_number_dataset_signal_region_data"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_VV"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_WJets0"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_STop"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar"+label+"_"+self.channel+"_mj").Print()

        number_dataset_signal_region_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_data"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_data"+label+"_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region: %s +/- sqrt(%s)"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj);

        self.file_out_ttbar_control.write("%s channel SF: \n"%(self.channel));
        self.file_out_ttbar_control.write("event number of data in signal_region: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj));

        number_dataset_signal_region_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_TotalMC"+label+"_"+self.channel+"_mj").getVal();

        print "event number of TotalMC %s in signal_region: %s +/- sqrt(%s) "%(label,number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj);
        self.file_out_ttbar_control.write("event number of TotalMC %s in signal_region: %s +/- sqrt(%s) \n"%(label,number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj));


        number_dataset_signal_region_before_cut_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_data"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_cut_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_data"+label+"_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region before_cut: %s +/- sqrt(%s)"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj);
        self.file_out_ttbar_control.write("event number of data in signal_region before_cut: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj));

        number_dataset_signal_region_before_cut_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_cut_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        print "event number of TotalMC %s in signal_region before_cut: %s +/- sqrt(%s) "%(label,number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj);
        self.file_out_ttbar_control.write("event number of TotalMC %s in signal_region before_cut: %s +/- sqrt(%s) \n"%(label,number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj));

                                                             
        # wtagger_eff reweight: only reweight the efficiency difference between MC and data
        wtagger_eff_MC = number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj;
        wtagger_eff_data= number_dataset_signal_region_data_mj/number_dataset_signal_region_before_cut_data_mj;

        wtagger_eff_reweight=wtagger_eff_data/wtagger_eff_MC;
        wtagger_eff_reweight_err=wtagger_eff_reweight*TMath.Sqrt(number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj +number_dataset_signal_region_before_cut_error2_data_mj/number_dataset_signal_region_before_cut_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_before_cut_error2_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj);
        
        print "wtagger efficiency of %s channel"%(self.channel )
        print "wtagger_eff_MC %s = %s "%(label,wtagger_eff_MC )
        print "wtagger_eff_data = %s "%(wtagger_eff_data )
        print "wtagger_eff_reweight %s = %s +/- %s"%(label,wtagger_eff_reweight, wtagger_eff_reweight_err)

        self.file_out_ttbar_control.write("wtagger_eff_MC %s = %s \n"%(label,wtagger_eff_MC ));
        self.file_out_ttbar_control.write("wtagger_eff_data = %s \n"%(wtagger_eff_data ));
        self.file_out_ttbar_control.write("wtagger_eff_reweight %s= %s +/- %s\n"%(label,wtagger_eff_reweight, wtagger_eff_reweight_err));


## method to do fail and pass fit for the scale factor evaluation
    def ScaleFactor_forPureWJet_TTbar_controlsample(self,in_file_name,label="",fit_or_not=1):

        print " ############################## Pass: dataset #################################### "
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_data_mj = self.workspace4fit_.data("rdataset_data"+label+"_"+self.channel+"_mj");
        rdataset_TTbar_mj = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+self.channel+"_mj");
        rdataset_STop_mj = self.workspace4fit_.data("rdataset_STop"+label+"_"+self.channel+"_mj");
        rdataset_VV_mj = self.workspace4fit_.data("rdataset_VV"+label+"_"+self.channel+"_mj");
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj)

        # inmport the pdf for each background
        print " ############################### Pass: model ################################### "
        model_histpdf_TTbar = self.workspace4fit_.pdf(rdataset_TTbar_mj.GetName()+"_histpdf")
        model_histpdf_STop = self.workspace4fit_.pdf(rdataset_STop_mj.GetName()+"_histpdf")
        model_histpdf_VV = self.workspace4fit_.pdf(rdataset_VV_mj.GetName()+"_histpdf")
        model_histpdf_WJets = self.workspace4fit_.pdf(rdataset_WJets_mj.GetName()+"_histpdf")
         
        number_TTbar = RooRealVar("rrv_number_TTbar"+label+"_"+self.channel ,"rrv_number_TTbar"+label+"_"+self.channel,rdataset_TTbar_mj.sumEntries());
        number_STop = RooRealVar("rrv_number_STop"+label+"_"+self.channel ,"rrv_number_STop"+label+"_"+self.channel,rdataset_STop_mj.sumEntries());
        number_VV = RooRealVar("rrv_number_VV"+label+"_"+self.channel ,"rrv_number_VV"+label+"_"+self.channel,rdataset_VV_mj.sumEntries());
        number_WJets = RooRealVar("rrv_number_WJets"+label+"_"+self.channel,"rrv_number_WJets"+label+"_"+self.channel,rdataset_WJets_mj.sumEntries());

        #### Total MC Model
        model_TTbar_STop_VV_WJets = RooAddPdf("model_TTbar_STop_VV_WJets"+label+"_"+self.channel,"model_TTbar_STop_VV_WJets"+label+"_"+self.channel,RooArgList(model_histpdf_TTbar, model_histpdf_STop, model_histpdf_VV, model_histpdf_WJets), RooArgList(number_TTbar, number_STop, number_VV, number_WJets) );
        getattr(self.workspace4fit_,"import")(model_TTbar_STop_VV_WJets);


        print " ############################## Fail: dataset #################################### "
        rdataset_data_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_TTbar_mj_fail = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_STop_mj_fail = self.workspace4fit_.data("rdataset_STop"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_VV_mj_fail = self.workspace4fit_.data("rdataset_VV"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_WJets_mj_fail = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj_fail);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj_fail)

        print " ############################### Fail: model ################################### "
        model_histpdf_TTbar_fail = self.workspace4fit_.pdf(rdataset_TTbar_mj_fail.GetName()+"_histpdf")
        model_histpdf_STop_fail = self.workspace4fit_.pdf(rdataset_STop_mj_fail.GetName()+"_histpdf")
        model_histpdf_VV_fail = self.workspace4fit_.pdf(rdataset_VV_mj_fail.GetName()+"_histpdf")
        model_histpdf_WJets_fail = self.workspace4fit_.pdf(rdataset_WJets_mj_fail.GetName()+"_histpdf")
 
        number_TTbar_fail = RooRealVar("rrv_number_TTbar_fail"+label+"_"+self.channel ,"rrv_number_TTbar_fail"+label+"_"+self.channel,rdataset_TTbar_mj_fail.sumEntries());
        number_STop_fail = RooRealVar("rrv_number_STop_fail"+label+"_"+self.channel ,"rrv_number_STop_fail"+label+"_"+self.channel,rdataset_STop_mj_fail.sumEntries());
        number_VV_fail = RooRealVar("rrv_number_VV_fail"+label+"_"+self.channel ,"rrv_number_VV_fail"+label+"_"+self.channel,rdataset_VV_mj_fail.sumEntries());
        number_WJets_fail = RooRealVar("rrv_number_WJets_fail"+label+"_"+self.channel,"rrv_number_WJets_fail"+label+"_"+self.channel,rdataset_WJets_mj_fail.sumEntries());

        model_TTbar_STop_VV_WJets_fail = RooAddPdf("model_TTbar_STop_VV_WJets_fail"+label+"_"+self.channel,"model_TTbar_STop_VV_WJets_fail"+label+"_"+self.channel, RooArgList(model_histpdf_TTbar_fail, model_histpdf_STop_fail, model_histpdf_VV_fail, model_histpdf_WJets_fail), RooArgList(number_TTbar_fail, number_STop_fail, number_VV_fail, number_WJets_fail) );
        getattr(self.workspace4fit_,"import")(model_TTbar_STop_VV_WJets_fail);

        ##### inclusive scale factor from total integral

        scale_number_TTbar_STop_VV_WJets=(rdataset_TTbar_mj.sumEntries()+rdataset_STop_mj.sumEntries()+rdataset_VV_mj.sumEntries() +rdataset_WJets_mj.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )
        scale_number_TTbar_STop_VV_WJets_fail=(rdataset_TTbar_mj_fail.sumEntries()+rdataset_STop_mj_fail.sumEntries()+rdataset_VV_mj_fail.sumEntries() +rdataset_WJets_mj_fail.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )

        rrv_scale_number_TTbar_STop_VV_WJets = RooRealVar("rrv_scale_number_TTbar_STop_VV_WJets"+label,"rrv_scale_number_TTbar_STop_VV_WJets"+label,scale_number_TTbar_STop_VV_WJets);
        rrv_scale_number_TTbar_STop_VV_WJets_fail = RooRealVar("rrv_scale_number_TTbar_STop_VV_WJets_fail"+label,"rrv_scale_number_TTbar_STop_VV_WJets_fail"+label,scale_number_TTbar_STop_VV_WJets_fail);
        getattr(self.workspace4fit_,"import")(rrv_scale_number_TTbar_STop_VV_WJets);
        getattr(self.workspace4fit_,"import")(rrv_scale_number_TTbar_STop_VV_WJets_fail);


        #### All the shape parameters and normalization are fixed

        print " ############################### Pass: Single MC model ################################### "
        
        model_STop = self.get_STop_mj_Model("_STop"+label);
        model_VV = self.get_VV_mj_Model("_VV"+label);
        model_WJets = self.get_General_mj_Model("_WJets0"+label);

        print " ############################### Fail: Single MC model ################################### "

        model_STop_fail = self.get_STop_mj_Model("_STop"+label+"_"+"failtau2tau1cut");
        model_VV_fail = self.get_VV_mj_Model("_VV"+label+"_"+"failtau2tau1cut");
        model_WJets_fail = self.get_General_mj_Model("_WJets0"+label+"_"+"failtau2tau1cut");

        self.constrainslist_data=[];
        ### Model for unmatched events passing and failing the cut --> ErfExp
        print " ############################### Pass: Data Bkg ################################### "
        model_bkg_data = self.make_Model("_bkg_data"+label,"ErfExp_ttbar","_mj",self.constrainslist_data); ## all the parameters contrained within the assigned error
        print " ############################### Fail: Data Bkg ################################### "
        model_bkg_data_fail = self.make_Model("_bkg_data"+label+"_"+"failtau2tau1cut","ErfExp_ttbar_failtau2tau1cut","_mj",self.constrainslist_data); 
## all the parameters contrained within the assigned error

        ### Model for matched events passing and failing the cut --> 2Gaus_ttbar and GausChebychev_ttbar_failtau2tau1cut
        print " ############################### Pass: Data ttbar resonant component ################################### "
        model_ttbar_data = self.make_Model_for_ttbar_controlsample("_ttbar_data"+label,"2Gaus_ttbar","_mj",self.constrainslist_data,label); ## No constraint
        print " ############################### Fail: Data ttbar resonant component ################################### "
        model_ttbar_data_fail = self.make_Model_for_ttbar_controlsample("_ttbar_data"+label+"_"+"failtau2tau1cut","GausChebychev_ttbar_failtau2tau1cut","_mj",self.constrainslist_data,label); 
## No Constraint


        print " ############################### Total Data Pdf Fail ################################### "
        model_data_fail=RooAddPdf("model_data"+label+"_"+"failtau2tau1cut"+"_"+self.channel,"model_data+"+label+"_"+"failtau2tau1cut"+"_"+self.channel, RooArgList(model_ttbar_data_fail, model_bkg_data_fail, model_STop_fail, model_VV_fail, model_WJets_fail));
        print " ############################### Total Data Pdf Pass ################################### "
        model_data=RooAddPdf("model_data"+label+"_"+self.channel,"model_data"+label+"_"+self.channel, RooArgList(model_ttbar_data, model_bkg_data, model_STop, model_VV, model_WJets));

        getattr(self.workspace4fit_,"import")(model_data);
        getattr(self.workspace4fit_,"import")(model_data_fail);

        ### take the event category from the workspace
        print " ############################### Data Event Category ################################### "
        category_p_f=self.workspace4fit_.cat("category_p_f"+label+"_"+self.channel);

        ### Define the simultaneous fit
        simPdf_data = RooSimultaneous("simPdf_data"+label+"_"+self.channel,"simPdf_data"+label+"_"+self.channel,category_p_f);
        simPdf_data.addPdf(model_data,"pass");
        simPdf_data.addPdf(model_data_fail,"fail");
        simPdf_data.Print();

        combData_p_f_data= self.workspace4fit_.data("combData_p_f_data"+label+"_"+self.channel);## combined dataset
        combData_p_f_data.Print();

        print " ############################### Data Total Fit ################################### "
        pdfconstrainslist_data=RooArgSet("pdfconstrainslist_data"+label+"_"+self.channel);
        for i in range(len(self.constrainslist_data)):
            self.workspace4fit_.pdf(self.constrainslist_data[i]).Print();
            pdfconstrainslist_data.add(self.workspace4fit_.pdf(self.constrainslist_data[i]) );
            
        if fit_or_not :
            rfresult_data=simPdf_data.fitTo(combData_p_f_data,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
            rfresult_data=simPdf_data.fitTo(combData_p_f_data,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
            rfresult_data.Print();


        #fit TotalMC
        self.constrainslist_TotalMC=[];

        print " ############################### Pass: MC Bkg ################################### "
        model_bkg_TotalMC = self.make_Model("_bkg_TotalMC"+label,"ErfExp_ttbar","_mj",self.constrainslist_TotalMC);
        print " ############################### Fail: MC Bkg ################################### "
        model_bkg_TotalMC_fail = self.make_Model("_bkg_TotalMC"+label+"_"+"failtau2tau1cut","ErfExp_ttbar_failtau2tau1cut","_mj",self.constrainslist_TotalMC);

        print " ############################### Pass: MC ttbar ################################### "
        model_ttbar_TotalMC = self.make_Model_for_ttbar_controlsample("_ttbar_TotalMC"+label,"2Gaus_ttbar","_mj",self.constrainslist_TotalMC,label);
        print " ############################### Fail: MC ttbar ################################### "
        model_ttbar_TotalMC_fail = self.make_Model_for_ttbar_controlsample("_ttbar_TotalMC"+label+"_"+"failtau2tau1cut","GausChebychev_ttbar_failtau2tau1cut","_mj",self.constrainslist_TotalMC,label);

        print " ############################### Fail: MC total PDF ################################### "
        model_TotalMC_fail = RooAddPdf("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+self.channel,"model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+self.channel, RooArgList(model_ttbar_TotalMC_fail, model_bkg_TotalMC_fail, model_STop_fail, model_VV_fail, model_WJets_fail));
        print " ############################### Pass: MC total PDF ################################### "
        model_TotalMC=RooAddPdf("model_TotalMC"+label+"_"+self.channel,"model_TotalMC"+label+"_"+self.channel, RooArgList(model_ttbar_TotalMC, model_bkg_TotalMC, model_STop, model_VV, model_WJets));
        getattr(self.workspace4fit_,"import")(model_TotalMC_fail);
        getattr(self.workspace4fit_,"import")(model_TotalMC);


        print " ############################### Data MC Category ################################### "
        simPdf_TotalMC=RooSimultaneous("simPdf_TotalMC"+label+"_"+self.channel,"simPdf_TotalMC"+label+"_"+self.channel,category_p_f);
        simPdf_TotalMC.addPdf(model_TotalMC,"pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail,"fail");
        
        combData_p_f_TotalMC=self.workspace4fit_.data("combData_p_f_TotalMC"+label+"_"+self.channel);


        print " ############################### MC Total Fit ################################### "
        pdfconstrainslist_TotalMC=RooArgSet("pdfconstrainslist_TotalMC"+label+"_"+self.channel);
        for i in range(len(self.constrainslist_TotalMC)):
            self.workspace4fit_.pdf(self.constrainslist_TotalMC[i]).Print();
            pdfconstrainslist_TotalMC.add(self.workspace4fit_.pdf(self.constrainslist_TotalMC[i]) );
        if fit_or_not :
            rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_p_f_TotalMC,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC), RooFit.Minimizer("Minuit2"));
            rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_p_f_TotalMC,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC), RooFit.Minimizer("Minuit2"));
            rfresult_TotalMC.Print();


        xframe_data = rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));
        xframe_data_fail = rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));
        xframe_data_extremefail=rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));

        ### plot data pass and fail
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        
        ### plot total pass mc normalizing it to data, passing + failing
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop.GetName(),model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("WJets"),RooFit.Components("%s"%(model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        ### plot total fail mc normalizing it to data, passing + failing
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("WJets"),RooFit.Components("%s"%( model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        ## plot data again
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## plot mc fit function
        combData_p_f_TotalMC.plotOn(xframe_data,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible() , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed),RooFit.Color(2))
        else:
         simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
              
        ## plot data fit function
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0));

        if TString(label).Contains("herwig"):
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
        else:
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(1))
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(1))

        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        #fail plots -> plot MC fit
        combData_p_f_TotalMC.plotOn(xframe_data_fail,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed),RooFit.Color(2))

        else:
         simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
              

        #fail plots -> plot data fit
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
        else:
         simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(1))
              
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        #signal window
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data.GetMaximum()); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kGray+2); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data.GetMaximum()); upperLine.SetLineWidth(2); upperLine.SetLineColor(kGray+2); upperLine.SetLineStyle(9);
        #xframe_data.addObject(lowerLine); xframe_data.addObject(upperLine);
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data_fail.GetMaximum()); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kGray+2); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data_fail.GetMaximum()); upperLine.SetLineWidth(2); upperLine.SetLineColor(kGray+2); upperLine.SetLineStyle(9);
        #xframe_data_fail.addObject(lowerLine); xframe_data_fail.addObject(upperLine);

        # make the legend
        leg_data = self.legend4Plot(xframe_data,0,1, 0.15)
        xframe_data.addObject(leg_data)
        leg_data_fail=self.legend4Plot(xframe_data_fail,0,1, 0.22)
        xframe_data_fail.addObject(leg_data_fail)

        #add mean and width
        rrv_mean_gaus_data = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_"+self.channel);
        rrv_sigma_gaus_data = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_"+self.channel);
        rrv_mean_gaus_TotalMC = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+self.channel);
        rrv_sigma_gaus_TotalMC = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+self.channel);

        if rrv_mean_gaus_TotalMC:
            tl_MC_mean =TLatex(0.25 ,0.62, ("Mean_{MC } = %3.1f #pm %2.1f")%(rrv_mean_gaus_TotalMC.getVal(), rrv_mean_gaus_TotalMC.getError()) );
            tl_MC_sigma =TLatex(0.25 ,0.57, ("Sigma_{MC }= %2.1f #pm %2.1f")%(rrv_sigma_gaus_TotalMC.getVal(), rrv_sigma_gaus_TotalMC.getError()) );
            tl_MC_mean.SetNDC(); tl_MC_sigma.SetNDC();
            tl_MC_mean.SetTextSize(0.03)
            tl_MC_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_MC_mean);
            # xframe_data.addObject(tl_MC_sigma);

        if rrv_mean_gaus_data:
            tl_data_mean =TLatex(0.25 ,0.52, ("Mean_{data} = %3.1f #pm %2.1f")%(rrv_mean_gaus_data.getVal(), rrv_mean_gaus_data.getError()) );
            tl_data_sigma=TLatex(0.25 ,0.47, ("Sigma_{data}= %2.1f #pm %2.1f")%(rrv_sigma_gaus_data.getVal(), rrv_sigma_gaus_data.getError()) );
            tl_data_mean.SetNDC(); tl_data_sigma.SetNDC();
            tl_data_mean.SetTextSize(0.03)
            tl_data_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_data_mean);
            # xframe_data.addObject(tl_data_sigma);
         
        xframe_data.GetYaxis().SetRangeUser(1e-2,xframe_data.GetMaximum()*1.2);
        xframe_data_fail.GetYaxis().SetRangeUser(1e-2,xframe_data_fail.GetMaximum()*1.2);


        self.draw_canvas(xframe_data,"plots_%s_%s_%s_g1/m_j_fitting_TTbar_controlsample_wtaggercut%s/"%(options.additioninformation, self.channel, self.wtagger_label, self.wtagger_label),"control%s_%s_%s_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));
        self.draw_canvas(xframe_data_fail,"plots_%s_%s_%s_g1/m_j_fitting_TTbar_controlsample_wtaggercut%s/"%(options.additioninformation, self.channel, self.wtagger_label,self.wtagger_label),"control%s_%s_%s_fail_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));


        self.ShowParam_Pdf(simPdf_data,RooArgSet(rrv_mass_j,category_p_f));
        self.ShowParam_Pdf(simPdf_TotalMC,RooArgSet(rrv_mass_j,category_p_f));
        
        if fit_or_not:
            rfresult_TotalMC.covarianceMatrix().Print();
            rfresult_data.covarianceMatrix().Print();
            rfresult_TotalMC.Print();
            rfresult_data.Print();

        self.ShowParam_Pdf(simPdf_data,RooArgSet(rrv_mass_j,category_p_f));


        # SF for LP
        print "########################## Extreme Fail Analysis ############################## "

        rdataset_data_mj_extremefail = self.workspace4fit_.data("rdataset_data"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_TotalMC_mj_extremefail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_TTbar_mj_extremefail = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_STop_mj_extremefail = self.workspace4fit_.data("rdataset_STop"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_VV_mj_extremefail = self.workspace4fit_.data("rdataset_VV"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_WJets_mj_extremefail = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj_extremefail);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj_extremefail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj_extremefail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj_extremefail)

        model_histpdf_TTbar_extremefail = self.workspace4fit_.pdf(rdataset_TTbar_mj_extremefail.GetName()+"_histpdf")
        model_histpdf_STop_extremefail = self.workspace4fit_.pdf(rdataset_STop_mj_extremefail.GetName()+"_histpdf")
        model_histpdf_VV_extremefail = self.workspace4fit_.pdf(rdataset_VV_mj_extremefail.GetName()+"_histpdf")
        model_histpdf_WJets_extremefail = self.workspace4fit_.pdf(rdataset_WJets_mj_extremefail.GetName()+"_histpdf")

        number_TTbar_extremefail = RooRealVar("rrv_number_TTbar"+label+"_"+"extremefail"+"_"+self.channel ,"rrv_number_TTbar"+label+"_"+"extremefail"+"_"+self.channel,rdataset_TTbar_mj_extremefail.sumEntries());
        number_STop_extremefail = RooRealVar("rrv_number_STop"+label+"_"+"extremefail"+"_"+self.channel ,"rrv_number_STop"+label+"_"+"extremefail"+"_"+self.channel,rdataset_STop_mj_extremefail.sumEntries());
        number_VV_extremefail = RooRealVar("rrv_number_VV"+label+"_"+"extremefail"+"_"+self.channel ,"rrv_number_VV"+label+"_"+"extremefail"+"_"+self.channel,rdataset_VV_mj_extremefail.sumEntries());
        number_WJets_extremefail = RooRealVar("rrv_number_WJets"+label+"_"+"extremefail"+"_"+self.channel,"rrv_number_WJets"+label+"_"+"extremefail"+"_"+self.channel,rdataset_WJets_mj_extremefail.sumEntries());
        
        model_TTbar_STop_VV_WJets_extremefail = RooAddPdf("model_TTbar_STop_VV_WJets"+label+"_"+"extremefail"+"_"+self.channel,"model_TTbar_STop_VV_WJets"+label+"_"+"extremefail""_"+self.channel, RooArgList(model_histpdf_TTbar_extremefail, model_histpdf_STop_extremefail, model_histpdf_VV_extremefail, model_histpdf_WJets_extremefail), RooArgList(number_TTbar_extremefail, number_STop_extremefail, number_VV_extremefail, number_WJets_extremefail) );
        getattr(self.workspace4fit_,"import")(model_TTbar_STop_VV_WJets_extremefail);

        scale_number_TTbar_STop_VV_WJets_extremefail = (rdataset_TTbar_mj_extremefail.sumEntries()+rdataset_STop_mj_extremefail.sumEntries()+rdataset_VV_mj_extremefail.sumEntries() +rdataset_WJets_mj_extremefail.sumEntries() )/( rdataset_data_mj_extremefail.sumEntries() )
        rrv_scale_number_TTbar_STop_VV_WJets_extremefail = RooRealVar("rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail","rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail",scale_number_TTbar_STop_VV_WJets_extremefail);
        getattr(self.workspace4fit_,"import")(rrv_scale_number_TTbar_STop_VV_WJets_extremefail);

        model_STop_extremefail = self.get_STop_mj_Model("_STop"+label+"_"+"extremefailtau2tau1cut");
        model_VV_extremefail = self.get_VV_mj_Model("_VV"+label+"_"+"extremefailtau2tau1cut");
        model_WJets_extremefail = self.get_General_mj_Model("_WJets0"+label+"_"+"extremefailtau2tau1cut");
                                                   

        tmp_constrainslist_data_LP=[];
        model_ttbar_data_extremefail = self.make_Model("_ttbar_data"+label+"_"+"extremefailtau2tau1cut","Exp_ttbar_extremefailtau2tau1cut","_mj",tmp_constrainslist_data_LP,50.);
        model_bkg_data_extremefail = self.make_Model("_bkg_data"+label+"_"+"extremefailtau2tau1cut","Exp_bkg_extremefailtau2tau1cut","_mj",tmp_constrainslist_data_LP,50);
        ## make the model for the extreme fail fit on data
        model_data_extremefail = RooAddPdf("model_data"+label+"_"+"extremefailtau2tau1cut_"+self.channel,"model_data"+label+"_"+"extremefailtau2tau1cut_"+self.channel, RooArgList(model_ttbar_data_extremefail, model_bkg_data_extremefail, model_STop_extremefail, model_VV_extremefail, model_WJets_extremefail));
        getattr(self.workspace4fit_,"import")(model_data_extremefail);

        ## constraint the parameter in the extreme fail fit
        pdfconstrainslist_data_LP = RooArgSet("pdfconstrainslist_data_LP"+label+"_"+self.channel);
        for i in range(len(tmp_constrainslist_data_LP)):
         self.workspace4fit_.pdf(tmp_constrainslist_data_LP[i]).Print();
         pdfconstrainslist_data_LP.add(self.workspace4fit_.pdf(tmp_constrainslist_data_LP[i]) );

        if fit_or_not :
          rfresult_data_LP = model_data_extremefail.fitTo(rdataset_data_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data_LP),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
          rfresult_data_LP=model_data_extremefail.fitTo(rdataset_data_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data_LP),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
          rfresult_data_LP.Print();

          rrv_number_bkg_data_extremefailtau2tau1cut = self.workspace4fit_.var("rrv_number_bkg_data%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel));
          rrv_number_ttbar_data_extremefailtau2tau1cut = self.workspace4fit_.var("rrv_number_ttbar_data%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel))
          rrv_number_ttbar_data_extremefailtau2tau1cut.setError( (rrv_number_ttbar_data_extremefailtau2tau1cut.getVal()+rrv_number_bkg_data_extremefailtau2tau1cut.getVal())/2. );
          rrv_number_bkg_data_extremefailtau2tau1cut.setError( (rrv_number_ttbar_data_extremefailtau2tau1cut.getVal()+rrv_number_bkg_data_extremefailtau2tau1cut.getVal())/2. );

          model_STop_extremefail.Print();
          model_VV_extremefail.Print();
          model_WJets_extremefail.Print();
          rdataset_data_mj_extremefail.Print();
          rdataset_TTbar_mj_extremefail.Print();
          rdataset_STop_mj_extremefail.Print();
          rdataset_VV_mj_extremefail.Print();
          rdataset_WJets_mj_extremefail.Print();

        ## extremefail fit of the mc
        tmp_constrainslist_TotalMC_LP=[];
        model_ttbar_TotalMC_extremefail = self.make_Model("_ttbar_TotalMC"+label+"_"+"extremefailtau2tau1cut","Exp_ttbar_extremefailtau2tau1cut","_mj",tmp_constrainslist_TotalMC_LP, 100.);
        model_bkg_TotalMC_extremefail = self.make_Model("_bkg_TotalMC"+label+"_"+"extremefailtau2tau1cut","Exp_bkg_extremefailtau2tau1cut","_mj",tmp_constrainslist_TotalMC_LP, 0.);

        model_TotalMC_extremefail=RooAddPdf("model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+self.channel,"model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+self.channel, RooArgList(model_ttbar_TotalMC_extremefail, model_bkg_TotalMC_extremefail, model_STop_extremefail, model_VV_extremefail, model_WJets_extremefail));
        getattr(self.workspace4fit_,"import")(model_TotalMC_extremefail);

        pdfconstrainslist_TotalMC_LP=RooArgSet("pdfconstrainslist_TotalMC_LP_"+self.channel);
        for i in range(len(tmp_constrainslist_TotalMC_LP)):
            self.workspace4fit_.pdf(tmp_constrainslist_TotalMC_LP[i]).Print();
            pdfconstrainslist_TotalMC_LP.add(self.workspace4fit_.pdf(tmp_constrainslist_TotalMC_LP[i]) );

        if fit_or_not :
                                 
           rfresult_TotalMC_LP=model_TotalMC_extremefail.fitTo(rdataset_TotalMC_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_LP),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
           rfresult_TotalMC_LP=model_TotalMC_extremefail.fitTo(rdataset_TotalMC_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_LP),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
           rfresult_TotalMC_LP.Print();
           rrv_number_bkg_TotalMC_extremefailtau2tau1cut=self.workspace4fit_.var("rrv_number_bkg_TotalMC%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel));
           rrv_number_ttbar_TotalMC_extremefailtau2tau1cut=self.workspace4fit_.var("rrv_number_ttbar_TotalMC%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel));
           rrv_number_ttbar_TotalMC_extremefailtau2tau1cut.setError( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut.getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut.getVal())/2. );
           rrv_number_bkg_TotalMC_extremefailtau2tau1cut.setError( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut.getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut.getVal())/2. );
           model_STop_extremefail.Print();
           model_VV_extremefail.Print();
           model_WJets_extremefail.Print();
           rdataset_TotalMC_mj_extremefail.Print();
                                                        
        ## plot the extreme fail fit result + distribution

        rdataset_data_mj_extremefail.plotOn(xframe_data_extremefail, RooFit.Name("data invisi"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_extremefail.GetName(), model_histpdf_VV_extremefail.GetName(), model_histpdf_WJets_extremefail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV_extremefail.GetName(), model_histpdf_WJets_extremefail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("WJets"),RooFit.Components("%s"%( model_histpdf_WJets_extremefail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())
        rdataset_data_mj_extremefail.plotOn(xframe_data_extremefail, RooFit.Name("data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        rdataset_TotalMC_mj_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("TotalMC invisi"), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig") :
         model_TotalMC_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("MC fit"),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed),RooFit.LineColor(2))
        else:
         model_TotalMC_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("MC fit"),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))

        rdataset_data_mj_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data invisi"), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig") :
         model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data fit"),RooFit.NormRange("controlsample_fitting_range"),RooFit.LineColor(2))
         model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data fit invisi"),RooFit.NormRange("controlsample_fitting_range"),RooFit.LineColor(2))
        else:
         model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data fit"),RooFit.NormRange("controlsample_fitting_range"))
         model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data fit invisi"),RooFit.NormRange("controlsample_fitting_range"))
              
        xframe_data_extremefail.GetYaxis().SetRangeUser(1e-2,xframe_data_extremefail.GetMaximum()*1.3);

        leg_data_extremefail = self.legend4Plot(xframe_data_extremefail,0,1, 0.15)
        xframe_data_extremefail.addObject(leg_data_extremefail)

        self.draw_canvas(xframe_data_extremefail,"plots_%s_%s_%s_g1/m_j_fitting_TTbar_controlsample_wtaggercut%s/"%(options.additioninformation, self.channel, self.wtagger_label,self.wtagger_label),"control%s_%s_%s_extremefail_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));
                                                                  

  
    ### function to be used for the simultaneous fit of electron and muon channel --> just drawing the results
    def draw_ScaleFactor_forPureWJet_TTbar_controlsample(self,in_file_name, label="",fit_or_not=1):

        #dataset after tau2tau1 cut
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_data_mj = self.workspace4fit_.data("rdataset_data"+label+"_"+self.channel+"_mj");
        rdataset_TTbar_mj = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+self.channel+"_mj");
        rdataset_STop_mj = self.workspace4fit_.data("rdataset_STop"+label+"_"+self.channel+"_mj");
        rdataset_VV_mj = self.workspace4fit_.data("rdataset_VV"+label+"_"+self.channel+"_mj");
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj)

        model_histpdf_TTbar = self.workspace4fit_.pdf(rdataset_TTbar_mj.GetName()+"_histpdf")
        model_histpdf_STop = self.workspace4fit_.pdf(rdataset_STop_mj.GetName()+"_histpdf")
        model_histpdf_VV = self.workspace4fit_.pdf(rdataset_VV_mj.GetName()+"_histpdf")
        model_histpdf_WJets = self.workspace4fit_.pdf(rdataset_WJets_mj.GetName()+"_histpdf")

        number_TTbar = RooRealVar("rrv_number_TTbar"+label+"_"+self.channel ,"rrv_number_TTbar"+label+"_"+self.channel,rdataset_TTbar_mj.sumEntries());
        number_STop = RooRealVar("rrv_number_STop"+label+"_"+self.channel ,"rrv_number_STop"+label+"_"+self.channel,rdataset_STop_mj.sumEntries());
        number_VV = RooRealVar("rrv_number_VV"+label+"_"+self.channel ,"rrv_number_VV"+label+"_"+self.channel,rdataset_VV_mj.sumEntries());
        number_WJets = RooRealVar("rrv_number_WJets"+label+"_"+self.channel,"rrv_number_WJets"+label+"_"+self.channel,rdataset_WJets_mj.sumEntries());

        model_TTbar_STop_VV_WJets = self.workspace4fit_.pdf("model_TTbar_STop_VV_WJets"+label+"_"+self.channel);

        #dataset fail tau2tau1 cut
        rdataset_data_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_TTbar_mj_fail = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_STop_mj_fail = self.workspace4fit_.data("rdataset_STop"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_VV_mj_fail = self.workspace4fit_.data("rdataset_VV"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_WJets_mj_fail = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
 
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj_fail);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj_fail)

        model_histpdf_TTbar_fail = self.workspace4fit_.pdf(rdataset_TTbar_mj_fail.GetName()+"_histpdf")
        model_histpdf_STop_fail = self.workspace4fit_.pdf(rdataset_STop_mj_fail.GetName()+"_histpdf")
        model_histpdf_VV_fail = self.workspace4fit_.pdf(rdataset_VV_mj_fail.GetName()+"_histpdf")
        model_histpdf_WJets_fail = self.workspace4fit_.pdf(rdataset_WJets_mj_fail.GetName()+"_histpdf")

        number_TTbar_fail = RooRealVar("rrv_number_TTbar_fail"+label+"_"+self.channel ,"rrv_number_TTbar_fail"+label+"_"+self.channel,rdataset_TTbar_mj_fail.sumEntries());
        number_STop_fail = RooRealVar("rrv_number_STop_fail"+label+"_"+self.channel ,"rrv_number_STop_fail"+label+"_"+self.channel,rdataset_STop_mj_fail.sumEntries());
        number_VV_fail = RooRealVar("rrv_number_VV_fail"+label+"_"+self.channel ,"rrv_number_VV_fail"+label+"_"+self.channel,rdataset_VV_mj_fail.sumEntries());
        number_WJets_fail = RooRealVar("rrv_number_WJets_fail"+label+"_"+self.channel,"rrv_number_WJets_fail"+label+"_"+self.channel,rdataset_WJets_mj_fail.sumEntries());

        model_TTbar_STop_VV_WJets_fail = self.workspace4fit_.pdf("model_TTbar_STop_VV_WJets_fail"+label+"_"+self.channel);

        ## rescale factor for plots
        scale_number_TTbar_STop_VV_WJets = (rdataset_TTbar_mj.sumEntries()+rdataset_STop_mj.sumEntries()+rdataset_VV_mj.sumEntries() +rdataset_WJets_mj.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )
        scale_number_TTbar_STop_VV_WJets_fail = (rdataset_TTbar_mj_fail.sumEntries()+rdataset_STop_mj_fail.sumEntries()+rdataset_VV_mj_fail.sumEntries() +rdataset_WJets_mj_fail.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )


        model_STop = self.get_STop_mj_Model("_STop"+label);
        model_VV = self.get_VV_mj_Model("_VV"+label);
        model_WJets = self.get_General_mj_Model("_WJets0"+label);

        model_STop_fail = self.get_STop_mj_Model("_STop"+label+"_"+"failtau2tau1cut");
        model_VV_fail = self.get_VV_mj_Model("_VV"+label+"_"+"failtau2tau1cut");
        model_WJets_fail = self.get_General_mj_Model("_WJets0"+label+"_"+"failtau2tau1cut");

        model_data_fail = self.workspace4fit_.pdf("model_data"+label+"_"+"failtau2tau1cut"+"_"+self.channel)
        model_data = self.workspace4fit_.pdf("model_data"+label+"_"+self.channel);

        category_p_f=self.workspace4fit_.cat("category_p_f"+label+"_"+self.channel);

        simPdf_data = RooSimultaneous("simPdf_data"+label+"_"+self.channel,"simPdf_data"+label+"_"+self.channel,category_p_f);
        simPdf_data.addPdf(model_data,"pass");
        simPdf_data.addPdf(model_data_fail,"fail");
        combData_p_f_data = self.workspace4fit_.data("combData_p_f_data"+label+"_"+self.channel);

        simPdf_data.Print();
        combData_p_f_data.Print();

        model_TotalMC_fail = self.workspace4fit_.pdf("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+self.channel);
        model_TotalMC = self.workspace4fit_.pdf("model_TotalMC"+label+"_"+self.channel);

        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC"+label+"_"+self.channel,"simPdf_TotalMC"+label+"_"+self.channel,category_p_f);
        simPdf_TotalMC.addPdf(model_TotalMC,"pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail,"fail");
        combData_p_f_TotalMC = self.workspace4fit_.data("combData_p_f_TotalMC"+label+"_"+self.channel);


        xframe_data=rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));
        xframe_data_fail=rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));

        ## plot data and mc on the frame
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        
        #pass plot
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("WJets"),RooFit.Components("%s"%(model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #pass plots
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("STop_line_invisible"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("VV_line_invisible"),RooFit.Components("%s,%s"%(model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("WJets_line_invisible"),RooFit.Components("%s"%(model_histpdf_WJets.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        #fail plot
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("WJets"),RooFit.Components("%s"%( model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #solid line
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("STop_line_invisible"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("VV_line_invisible"),RooFit.Components("%s,%s"%(model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("WJets_line_invisible"),RooFit.Components("%s"%( model_histpdf_WJets_fail.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())


        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## plot mc fit
        combData_p_f_TotalMC.plotOn(xframe_data,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible() , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed),RooFit.Color(2))
        else:
         simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))


        ## plot data fit
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
        else:
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
         simPdf_data.plotOn(xframe_data,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))

        combData_p_f_TotalMC.plotOn(xframe_data_fail,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed),RooFit.Color(2))
        else:
         simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("MC fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
              

        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
         simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"),RooFit.Color(2))
        else:
         simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
         simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))

        #signal window
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data.GetMaximum()); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kGray+2); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data.GetMaximum()); upperLine.SetLineWidth(2); upperLine.SetLineColor(kGray+2); upperLine.SetLineStyle(9);
        xframe_data.addObject(lowerLine); xframe_data.addObject(upperLine);
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data_fail.GetMaximum()); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kGray+2); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data_fail.GetMaximum()); upperLine.SetLineWidth(2); upperLine.SetLineColor(kGray+2); upperLine.SetLineStyle(9);
        xframe_data_fail.addObject(lowerLine); xframe_data_fail.addObject(upperLine);

        #legend
        leg_data=self.legend4Plot(xframe_data,0,1, 0.15)
        xframe_data.addObject(leg_data)
        leg_data_fail=self.legend4Plot(xframe_data_fail,0,1,0.18)
        xframe_data_fail.addObject(leg_data_fail)

        #add mean and width
        tmp_channel="el";
        if self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_"+"el"): tmp_channel="el";
        else: tmp_channel="mu";

        rrv_mean_gaus_data =self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_"+tmp_channel);
        rrv_sigma_gaus_data =self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_"+tmp_channel);
        rrv_mean_gaus_TotalMC =self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel);
        rrv_sigma_gaus_TotalMC=self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel);

        if rrv_mean_gaus_TotalMC:
            tl_MC_mean =TLatex(0.25 ,0.62, ("Mean_{MC } = %3.1f #pm %2.1f")%(rrv_mean_gaus_TotalMC.getVal(), rrv_mean_gaus_TotalMC.getError()) );
            tl_MC_sigma =TLatex(0.25 ,0.57, ("Sigma_{MC }= %2.1f #pm %2.1f")%(rrv_sigma_gaus_TotalMC.getVal(), rrv_sigma_gaus_TotalMC.getError()) );
            tl_MC_mean.SetNDC(); tl_MC_sigma.SetNDC();
            tl_MC_mean.SetTextSize(0.03)
            tl_MC_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_MC_mean);
            # xframe_data.addObject(tl_MC_sigma);

        if rrv_mean_gaus_data:
            tl_data_mean =TLatex(0.25 ,0.52, ("Mean_{data} = %3.1f #pm %2.1f")%(rrv_mean_gaus_data.getVal(), rrv_mean_gaus_data.getError()) );
            tl_data_sigma=TLatex(0.25 ,0.47, ("Sigma_{data}= %2.1f #pm %2.1f")%(rrv_sigma_gaus_data.getVal(), rrv_sigma_gaus_data.getError()) );
            tl_data_mean.SetNDC(); tl_data_sigma.SetNDC();
            tl_data_mean.SetTextSize(0.03)
            tl_data_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_data_mean); xframe_data.addObject(tl_data_sigma);
         
        xframe_data.GetYaxis().SetRangeUser(1e-2,xframe_data.GetMaximum()*1.2);
        xframe_data_fail.GetYaxis().SetRangeUser(1e-2,xframe_data_fail.GetMaximum()*1.2);

        self.draw_canvas(xframe_data,"plots_%s_%s_%s_g1/m_j_fitting_TTbar_controlsample_wtaggercut%s/"%(options.additioninformation, self.channel,self.wtagger_label, self.wtagger_label),"control%s_%s_%s_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));
        self.draw_canvas(xframe_data_fail,"plots_%s_%s_%s_g1/m_j_fitting_TTbar_controlsample_wtaggercut%s/"%(options.additioninformation, self.channel, self.wtagger_label,self.wtagger_label),"control%s_%s_%s_fail_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));



    ########## ---------------------------------------------------
    def get_mj_and_mlvj_dataset_TTbar_controlsample(self,in_file_name, label, jet_mass ="mJJNoKinFit"):# to get the shape of m_lvj jet_mass="jet_mass_pr" , "ttb_ca8_mass_pr"

        # read tree
        fileIn_name=TString(options.inPath+"/"+self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("SelectedCandidatesPlain");

        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ### dataset of m_j before tau2tau1 cut : Passed
        rdataset_mj = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j before tau2tau1 cut : Total
        rdataset_beforetau2tau1cut_mj = RooDataSet("rdataset"+label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_beforetau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j failed tau2tau1 cut :
        rdataset_failtau2tau1cut_mj = RooDataSet("rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_failtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j extreme failed tau2tau1 cut: >0.75
        rdataset_extremefailtau2tau1cut_mj = RooDataSet("rdataset"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_extremefailtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
                                                                        
          
        #### combine of dataset before and after tau2tau1 cut
        if TString(label).Contains("herwig"):
         category_cut = RooCategory("category_cut"+"_herwig_"+self.channel,"category_cut"+"_herwig_"+self.channel);
         category_cut.defineType("cut",1);
         category_cut.defineType("beforecut",2);
         combData4cut = RooDataSet("combData4cut"+"_herwig_"+self.channel,"combData4cut"+"_"+self.channel,RooArgSet(rrv_mass_j, category_cut, rrv_weight),RooFit.WeightVar(rrv_weight) );
        else:
         category_cut = RooCategory("category_cut"+"_"+self.channel,"category_cut"+"_"+self.channel);
         category_cut.defineType("cut",1);
         category_cut.defineType("beforecut",2);
         combData4cut = RooDataSet("combData4cut"+"_"+self.channel,"combData4cut"+"_"+self.channel,RooArgSet(rrv_mass_j, category_cut, rrv_weight),RooFit.WeightVar(rrv_weight) );
            
        ### create the category for the fit
        if TString(label).Contains("herwig"):
         if self.workspace4fit_.cat("category_p_f"+"_herwig_"+self.channel):
             category_p_f = self.workspace4fit_.cat("category_p_f"+"_herwig_"+self.channel);
         else:
             category_p_f = RooCategory("category_p_f"+"_herwig_"+self.channel,"category_p_f"+"_herwig_"+self.channel);
             category_p_f.defineType("pass");
             category_p_f.defineType("fail");
             getattr(self.workspace4fit_,"import")(category_p_f);

        else:
         if self.workspace4fit_.cat("category_p_f"+"_"+self.channel):
             category_p_f = self.workspace4fit_.cat("category_p_f"+"_"+self.channel);
         else:
             category_p_f = RooCategory("category_p_f"+"_"+self.channel,"category_p_f"+"_"+self.channel);
             category_p_f.defineType("pass");
             category_p_f.defineType("fail");
             getattr(self.workspace4fit_,"import")(category_p_f);
        
        combData_p_f=RooDataSet("combData_p_f"+label+"_"+self.channel,"combData_p_f"+label+"_"+self.channel,RooArgSet(rrv_mass_j, category_p_f, rrv_weight),RooFit.WeightVar(rrv_weight));
         
        print "N entries: ", treeIn.GetEntries()

        hnum_4region=TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_4region_error2=TH1D("hnum_4region_error2"+label+"_"+self.channel,"hnum_4region_error2"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_4region_before_cut=TH1D("hnum_4region_before_cut"+label+"_"+self.channel,"hnum_4region_before_cut"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_4region_before_cut_error2=TH1D("hnum_4region_before_cut_error2"+label+"_"+self.channel,"hnum_4region_before_cut_error2"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_2region=TH1D("hnum_2region"+label+"_"+self.channel,"hnum_2region"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total
        hnum_2region_error2=TH1D("hnum_2region_error2"+label+"_"+self.channel,"hnum_2region_error2"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total

        if self.channel=="el": tmp_lumi=19531.85;
        else: tmp_lumi=19538.85;
        for i in range(treeIn.GetEntries()):
            if i % 100000 == 0: print "i: ",i
            treeIn.GetEntry(i);
    
            discriminantCut = 0;# -1 -> fail wtagger; 1: pass wtagger; other: not belong to our selection
            lep_type =-1;
            if self.channel=="el": lep_type=0;
            if self.channel=="mu": lep_type=1;

            if treeIn.lep== lep_type and treeIn.nXjets == 1:
                discriminantCut = -1;

                if (self.categoryID == 0 or self.categoryID == 2):# LP
                    if treeIn.vTagPurity == 0: discriminantCut= 1;
                    else: discriminantCut = -1;

                if (self.categoryID == 1 or self.categoryID == 3):# HP
                    if treeIn.vTagPurity == 1: discriminantCut= 1;
                    else: discriminantCut = -1;

            tmp_jet_mass=getattr(treeIn, jet_mass);

            if i==0: tmp_scale_to_lumi=treeIn.LumiWeight*tmp_lumi;
            if not label=="_data":
                tmp_event_weight = treeIn.weight*tmp_lumi;
                tmp_event_weight4fit = treeIn.HLTweight*treeIn.PUweight*treeIn.GenWeight*treeIn.BTagWeight*treeIn.VTagWeight;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.LumiWeight*tmp_lumi/tmp_scale_to_lumi # for multi-sample, like STop and VV. There are two sample, and two wSampleWeight_value.Use the least wSampleWeight as scale.
            else:
                tmp_event_weight = 1.;
                tmp_event_weight4fit = 1.;


            if discriminantCut >0 and treeIn.mZZ> self.mass_lvj_min and treeIn.mZZ< self.mass_lvj_max and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() :

                rrv_mass_j.setVal( tmp_jet_mass );
                rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max:
                    hnum_4region.Fill(-1,tmp_event_weight );
                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                    hnum_4region.Fill(0,tmp_event_weight);
                    hnum_4region_error2.Fill(0,tmp_event_weight*tmp_event_weight);
                if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max:
                    hnum_4region.Fill(1,tmp_event_weight);
                hnum_4region.Fill(2,tmp_event_weight);

                category_cut.setLabel("cut");
                combData4cut.add(RooArgSet(rrv_mass_j,category_cut),tmp_event_weight4fit);
                category_p_f.setLabel("pass");
                combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight);

            if (discriminantCut >0 or discriminantCut <0) and treeIn.mZZ> self.mass_lvj_min and treeIn.mZZ< self.mass_lvj_max and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() :

                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                    hnum_4region_before_cut.Fill(0,tmp_event_weight);
                    hnum_4region_before_cut_error2.Fill(0,tmp_event_weight*tmp_event_weight);
                rrv_mass_j.setVal( tmp_jet_mass );
                rdataset_beforetau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_beforetau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                category_cut.setLabel("beforecut");
                combData4cut.add(RooArgSet(rrv_mass_j,category_cut),tmp_event_weight4fit);

            if discriminantCut <0 and treeIn.mZZ> self.mass_lvj_min and treeIn.mZZ< self.mass_lvj_max and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() :

                rrv_mass_j.setVal( tmp_jet_mass );
                rdataset_failtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_failtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                category_p_f.setLabel("fail");
                combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight);

                if treeIn.vTagPurity == -1:
                    rdataset_extremefailtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                    rdataset4fit_extremefailtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi);

        rrv_scale_to_lumi_failtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,tmp_scale_to_lumi);
        rrv_scale_to_lumi_extremefailtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,tmp_scale_to_lumi);
                
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_failtau2tau1cut)
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_extremefailtau2tau1cut)

        #prepare m_j dataset

        rrv_number_dataset_sb_lo_mj = RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));
        rrv_number_dataset_signal_region_error2_mj = RooRealVar("rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj",hnum_4region_error2.GetBinContent(2));

        rrv_number_dataset_signal_region_before_cut_mj = RooRealVar("rrv_number_dataset_signal_region_before_cut"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut"+label+"_"+self.channel+"_mj",hnum_4region_before_cut.GetBinContent(2));
        rrv_number_dataset_signal_region_before_cut_error2_mj = RooRealVar("rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj",hnum_4region_before_cut_error2.GetBinContent(2));
        rrv_number_dataset_sb_hi_mj = RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));

        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_error2_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_error2_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj)
        getattr(self.workspace4fit_,"import")(combData_p_f);
        
        print "N_rdataset_mj: "
        getattr(self.workspace4fit_,"import")(rdataset_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj)
        getattr(self.workspace4fit_,"import")(rdataset_beforetau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_beforetau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset_failtau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_failtau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset_extremefailtau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_extremefailtau2tau1cut_mj)

                
        rdataset_mj.Print();
        rdataset4fit_mj.Print();
        rrv_number_dataset_sb_lo_mj.Print()
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_signal_region_error2_mj.Print()
        rrv_number_dataset_signal_region_before_cut_mj.Print()
        rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
        rrv_number_dataset_sb_hi_mj.Print()

        rdataset_mj.Print();
        rdataset_beforetau2tau1cut_mj.Print();
        rdataset_failtau2tau1cut_mj.Print();
        rdataset_extremefailtau2tau1cut_mj.Print();
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_signal_region_error2_mj.Print()
        rrv_number_dataset_signal_region_before_cut_mj.Print()
        rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
        combData_p_f.Print("v");




    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

       hpull = mplot_orig.pullHist();
       x = ROOT.Double(0.); y = ROOT.Double(0) ;
       for ipoint in range(0,hpull.GetN()):
        hpull.GetPoint(ipoint,x,y);
        if(y == 0):
          hpull.SetPoint(ipoint,x,10)
       
       mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
       medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
       mplot_pull.addObject(medianLine);
       mplot_pull.addPlotable(hpull,"P");
       mplot_pull.SetTitle("PULL");
       mplot_pull.GetXaxis().SetTitle("");
       mplot_pull.GetYaxis().SetRangeUser(-5,5);
       mplot_pull.GetYaxis().SetTitleSize(0.10);
       mplot_pull.GetYaxis().SetLabelSize(0.10);
       mplot_pull.GetXaxis().SetTitleSize(0.10);
       mplot_pull.GetXaxis().SetLabelSize(0.10);
       mplot_pull.GetYaxis().SetTitleOffset(0.40);
       mplot_pull.GetYaxis().SetTitle("#frac{data-fit}{#sigma_{data}}");
       mplot_pull.GetYaxis().CenterTitle();

       return mplot_pull;

    #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull=0):

      if iswithpull:
       if self.channel=="el":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
       elif self.channel=="mu":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
       banner.SetNDC(); banner.SetTextSize(0.04);
      else:
       if self.channel=="el":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
       if self.channel=="mu":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
       banner.SetNDC(); banner.SetTextSize(0.033);
                                                                                                         
      return banner;

    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, xoffset=0., yoffset=0., x_right_offset=0., y_upper_offset=0.):

        if left==-1:
            theLeg = TLegend(0.65+xoffset, 0.65+yoffset, 0.95+xoffset, 0.87+yoffset, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetBorderSize(0);
            theLeg.SetLineColor(0);
            theLeg.SetFillColor(0);
            if isFill:
                theLeg.SetFillStyle(1001);
            else:
                theLeg.SetFillStyle(0);
            theLeg.SetLineWidth(0);
            theLeg.SetLineStyle(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.1+xoffset, 0.65+yoffset, 0.8+xoffset+x_right_offset, 0.92+yoffset+y_upper_offset, "", "NDC");
            theLeg.SetFillColor(0);
            theLeg.SetFillStyle(0);
            theLeg.SetNColumns(2);
            theLeg.SetTextSize(.033);

        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        
        entryCnt = 0;
        objName_before = "";

        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            print objName;
            if not (((plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data() =="STop" : theLeg.AddEntry(theObj, "single Top","F");
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ/ZZ","F");
                    elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F");
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;

        return theLeg;

    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, mplot, mplot_pull,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0):# mplot + pull

        mplot.GetXaxis().SetTitleOffset(1.1);
        mplot.GetYaxis().SetTitleOffset(1.3);
        mplot.GetXaxis().SetTitleSize(0.05);
        mplot.GetYaxis().SetTitleSize(0.05);
        mplot.GetXaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetLabelSize(0.045);
        mplot_pull.GetXaxis().SetLabelSize(0.15);
        mplot_pull.GetYaxis().SetLabelSize(0.15);
        mplot_pull.GetYaxis().SetTitleSize(0.15);
        mplot_pull.GetYaxis().SetNdivisions(205);

                                                                          
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
         pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
         pad3=TPad("pad3","pad3",0.8,0.,1,1);
         pad1.Draw();
         pad2.Draw();
         pad3.Draw();
        else:
         pad1=TPad("pad1","pad1",0.,0. ,0.99,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.99,1. );
         pad1.Draw();
         pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
        mplot.Draw();
        banner = self.banner4Plot(1);
        banner.Draw();

        pad1.cd();
        mplot_pull.Draw();

        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

        ## create the directory where store the plots
        Directory = TString(in_directory+self.signal_sample+"_%02d_%02d/"%(options.cprime,options.BRnew));        
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
        else:
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","");
            rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());
        
        rlt_file.ReplaceAll(".pdf",".root");
        cMassFit.SaveAs(rlt_file.Data());

        string_file_name = TString(in_file_name);
        if string_file_name.EndsWith(".root"):
            string_file_name.ReplaceAll(".root","_"+in_model_name);
        else:
            string_file_name.ReplaceAll(".root","");
            string_file_name.Append("_"+in_model_name);

        if logy:
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy);

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0):

        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);

        if is_range:
            h2=TH2D("h2","",100,400,1400,4,0.00001,4);
            h2.Draw();
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045);
        in_obj.GetXaxis().SetTitleOffset(1.15);
        in_obj.GetXaxis().SetLabelSize(0.04);

        in_obj.GetYaxis().SetTitleSize(0.045);
        in_obj.GetYaxis().SetTitleOffset(1.40);
        in_obj.GetYaxis().SetLabelSize(0.04);

        banner = self.banner4Plot();
        banner.Draw();
        
        Directory=TString(in_directory+self.signal_sample+"_%02d_%02d/"%(options.cprime,options.BRnew));
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
        else:
            rlt_file.ReplaceAll(".root","");
            rlt_file = rlt_file.Append(".png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".pdf",".root");
        cMassFit.SaveAs(rlt_file.Data());

        if logy:
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());



    ### Get the liminosity value
    def GetLumi(self):

        if self.channel=="el": return 19.5
        if self.channel=="mu": return 19.5


    #### defines two different way to fit depending on pythia or herwig analysis
    def fit_TTbar_controlsample(self, isherwig=0):

      if isherwig ==0:

        print " Herwig 0 --> start fit on Pythia "
        print "fit_TTbar_controlsample --> Pythia samples"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");

        print "##################################################"
        print "############### Single Top DataSet ###############"
        print "##################################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop");

        if self.channel=="mu":
         self.fit_mj_single_MC(self.file_STop_mc,"_STop","ErfExpGaus_sp","_TTbar_controlsample"); ## Erf*exp + Gaus pass sample
         self.fit_mj_single_MC(self.file_STop_mc,"_STop_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp+Gaus fail sample
        else :
         self.fit_mj_single_MC(self.file_STop_mc,"_STop","ExpGaus","_TTbar_controlsample"); ## Exp + Gaus pass sample
         self.fit_mj_single_MC(self.file_STop_mc,"_STop_failtau2tau1cut","ExpGaus","_TTbar_controlsample"); ## Exp+Gaus fail sample
            
        self.fit_mj_single_MC(self.file_STop_mc,"_STop_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp for extreme fail

        ### Build WJet fit pass and fail distributions
        print "###########################################"
        print "############### WJets Pythia ##############"
        print "###########################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets0_mc,"_WJets0");

        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0","ErfExp","_TTbar_controlsample"); ## Erf*Exp pass sample
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp fail sample
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp extreme fail


        ### Build VV fit pass and fail distributions
        print "#########################################"
        print "############### VV Pythia ###############"
        print "#########################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV");

        if self.channel == "mu" :
              self.fit_mj_single_MC(self.file_VV_mc,"_VV","ErfExpGaus_sp","_TTbar_controlsample"); ## Exp+Gaus pass
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_failtau2tau1cut","ExpGaus","_TTbar_controlsample"); ## Exp +Gaus fail
        else:
              self.fit_mj_single_MC(self.file_VV_mc,"_VV","Gaus","_TTbar_controlsample"); ## Exp+Gaus pass
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp +Gaus fail

        self.fit_mj_single_MC(self.file_VV_mc,"_VV_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp


        print "#########################################"
        print "############# TTbar Powheg ##############"
        print "#########################################"
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_TTbar_mc,"_TTbar");

        print "################################################"
        print "############## Pseudo Data Powheg ##############"
        print "################################################"
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_pseudodata,"_TotalMC");


        print "#################################"
        print "############# Data ##############"
        print "#################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_data,"_data");
        self.fit_mj_TTbar_controlsample(self.file_data);
        self.ScaleFactor_forPureWJet_TTbar_controlsample(self.file_data);

      else:

          print " Herwig 1 --> start fit on Herwig "
          print "###############################################"
          print "############# Single Top DataSet ##############"
          print "###############################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop_herwig");
          if self.channel=="mu":
           self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig","ErfExpGaus_sp","_TTbar_controlsample");
          else :
           self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig","ExpGaus","_TTbar_controlsample");

          self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig_failtau2tau1cut","ExpGaus","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig_extremefailtau2tau1cut","Exp","_TTbar_controlsample");

          print "##########################################"
          print "############## WJets Herwig ##############"
          print "##########################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets0_mc,"_WJets0_herwig");

          self.fit_mj_single_MC(self.file_WJets1_mc,"_WJets0_herwig","ErfExp","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_WJets1_mc,"_WJets0_herwig_failtau2tau1cut","Exp","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_WJets1_mc,"_WJets0_herwig_extremefailtau2tau1cut","Exp","_TTbar_controlsample");


          print "##########################################"
          print "################ VV Pythia ###############"
          print "##########################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV_herwig");
          self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig","ExpGaus","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig_failtau2tau1cut","ExpGaus","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig_extremefailtau2tau1cut","Exp","_TTbar_controlsample");

          print "##########################################"
          print "################ TTbar Herwig ############"
          print "##########################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_TTbar_herwig,"_TTbar_herwig");
 
          print "##############################################"
          print "############# Pseudo Data Herwig #############"
          print "##############################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_pseudodata_herwig,"_TotalMC_herwig");

          print "##################################"
          print "############## Data ##############"
          print "##################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_data,"_data_herwig");
          self.fit_mj_TTbar_controlsample(self.file_data,"_herwig");
          self.ScaleFactor_forPureWJet_TTbar_controlsample(self.file_data,"_herwig");

class doFit_wj_and_wlvj_simultaneous:
    def __init__(self,isherwig=0):

        label = "";
        if isherwig==1 : label = "_herwig" ;

        self.workspace4fit_ = RooWorkspace("workspace4fit"+label+"_","workspace4fit"+label+"_");

        print " Simuiltaneous class -> define electron fitter "
        self.boostedW_fitter_el=doFit_wj_and_wlvj("el","ggH600",40,130,label, self.workspace4fit_)
        print " Simuiltaneous class -> define muon fitter "
        self.boostedW_fitter_mu=doFit_wj_and_wlvj("mu","ggH600",40,130,label, self.workspace4fit_)

        print " Fit Electron ----> "
        self.boostedW_fitter_el.fit_TTbar_controlsample(isherwig);
        print " Fit Muon ----> "
        self.boostedW_fitter_mu.fit_TTbar_controlsample(isherwig);

        self.workspace4fit_.data("rdataset_data"+label+"_mu_mj").Print();
        self.workspace4fit_.data("rdataset_data"+label+"_el_mj").Print();

        #### Define simultaneusly 4 category
        sample_type=RooCategory("sample_type"+label,"sample_type"+label);
        sample_type.defineType("mu_pass");
        sample_type.defineType("mu_fail");
        sample_type.defineType("el_pass");
        sample_type.defineType("el_fail");
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 


        rdataset_data_mu_mj      = self.workspace4fit_.data("rdataset_data"+label+"_mu_mj");
        rdataset_data_el_mj      = self.workspace4fit_.data("rdataset_data"+label+"_el_mj");
        rdataset_data_mu_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_mu_mj"); 
        rdataset_data_el_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_el_mj"); 
 
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j");

        combData_data=RooDataSet("combData_data"+label,"combData_data"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_data_mu_mj),RooFit.Import("el_pass",rdataset_data_el_mj),RooFit.Import("mu_fail",rdataset_data_mu_mj_fail),RooFit.Import("el_fail",rdataset_data_el_mj_fail) );
        combData_data.Print();

        rdataset_TotalMC_mu_mj = self.workspace4fit_.data("rdataset_TotalMC"+label+"_mu_mj");
        rdataset_TotalMC_el_mj = self.workspace4fit_.data("rdataset_TotalMC"+label+"_el_mj");
        rdataset_TotalMC_mu_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_mu_mj"); 
        rdataset_TotalMC_el_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_el_mj"); 

        combData_TotalMC = RooDataSet("combData_TotalMC"+label,"combData_TotalMC"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_TotalMC_mu_mj),RooFit.Import("el_pass",rdataset_TotalMC_el_mj),RooFit.Import("mu_fail",rdataset_TotalMC_mu_mj_fail),RooFit.Import("el_fail",rdataset_TotalMC_el_mj_fail) );
        combData_TotalMC.Print();

        # fit data
        model_data_mu      = self.workspace4fit_.pdf("model_data"+label+"_mu");
        model_data_fail_mu = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_mu");
        model_data_el      = self.workspace4fit_.pdf("model_data"+label+"_el");
        model_data_fail_el = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_el");

        simPdf_data = RooSimultaneous("simPdf_data_em"+label,"simPdf_data_em"+label,sample_type);
        simPdf_data.addPdf(model_data_mu,"mu_pass");
        simPdf_data.addPdf(model_data_el,"el_pass");
        simPdf_data.addPdf(model_data_fail_mu,"mu_fail");
        simPdf_data.addPdf(model_data_fail_el,"el_fail");

        constrainslist_data_em=self.boostedW_fitter_el.constrainslist_data +self.boostedW_fitter_mu.constrainslist_data
        pdfconstrainslist_data_em=RooArgSet("pdfconstrainslist_data_em"+label);
        for i in range(len(constrainslist_data_em)):
            self.workspace4fit_.pdf(constrainslist_data_em[i]).Print();
            pdfconstrainslist_data_em.add(self.workspace4fit_.pdf(constrainslist_data_em[i]) );
        pdfconstrainslist_data_em.Print();

        rfresult_data=simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em) )
        rfresult_data=simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em) )

        # fit TotalMC
        model_TotalMC_mu      = self.workspace4fit_.pdf("model_TotalMC"+label+"_mu");
        model_TotalMC_fail_mu = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_mu");
        model_TotalMC_el      = self.workspace4fit_.pdf("model_TotalMC"+label+"_el");
        model_TotalMC_fail_el = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_el");

        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC_em"+label,"simPdf_TotalMC_em"+label,sample_type);
        simPdf_TotalMC.addPdf(model_TotalMC_mu,"mu_pass");
        simPdf_TotalMC.addPdf(model_TotalMC_el,"el_pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail_mu,"mu_fail");
        simPdf_TotalMC.addPdf(model_TotalMC_fail_el,"el_fail");

        constrainslist_TotalMC_em=self.boostedW_fitter_el.constrainslist_TotalMC +self.boostedW_fitter_mu.constrainslist_TotalMC
        pdfconstrainslist_TotalMC_em=RooArgSet("pdfconstrainslist_TotalMC_em"+label);

        for i in range(len(constrainslist_TotalMC_em)):
            self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]).Print();
            pdfconstrainslist_TotalMC_em.add(self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]) );
        pdfconstrainslist_TotalMC_em.Print();

        rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em) )
        rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em) )

        self.boostedW_fitter_el.draw_ScaleFactor_forPureWJet_TTbar_controlsample(self.boostedW_fitter_el.file_data,label,0);
        self.boostedW_fitter_mu.draw_ScaleFactor_forPureWJet_TTbar_controlsample(self.boostedW_fitter_mu.file_data,label,0);


        rfresult_TotalMC.Print();
        rfresult_data.Print();
         
        rrv_eff_MC_el   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_el");
        rrv_eff_MC_mu   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_mu");
        rrv_mean_MC_el  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_el");
        rrv_sigma_MC_el = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_el");

        rrv_eff_data_el   = self.workspace4fit_.var("eff_ttbar_data"+label+"_el");
        rrv_eff_data_mu   = self.workspace4fit_.var("eff_ttbar_data"+label+"_mu");
        rrv_mean_data_el  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_el");
        rrv_sigma_data_el = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_el");

        rrv_eff_MC_el.Print()   ; rrv_eff_MC_mu.Print();
        rrv_eff_data_el.Print() ; rrv_eff_data_mu.Print();
        rrv_mean_MC_el.Print()  ; rrv_mean_data_el.Print();  
        rrv_sigma_MC_el.Print() ; rrv_sigma_data_el.Print();  

        pure_wtagger_sf_el=rrv_eff_data_el.getVal()/ rrv_eff_MC_el.getVal(); 
        pure_wtagger_sf_mu=rrv_eff_data_mu.getVal()/ rrv_eff_MC_mu.getVal(); 
        pure_wtagger_mean_shift_el= rrv_mean_data_el.getVal()-rrv_mean_MC_el.getVal();
        pure_wtagger_sigma_enlarge_el= rrv_sigma_data_el.getVal()/rrv_sigma_MC_el.getVal();

        pure_wtagger_sf_el_err= ( (rrv_eff_data_el.getError()/rrv_eff_data_el.getVal())**2 + (rrv_eff_MC_el.getError()/rrv_eff_MC_el.getVal())**2 )**0.5* pure_wtagger_sf_el

        pure_wtagger_sf_mu_err= ( (rrv_eff_data_mu.getError()/rrv_eff_data_mu.getVal())**2 + (rrv_eff_MC_mu.getError()/rrv_eff_MC_mu.getVal())**2 )**0.5* pure_wtagger_sf_mu
        
        pure_wtagger_mean_shift_err_el= ( rrv_mean_data_el.getError()**2 + rrv_mean_MC_el.getError()**2 )**0.5

        pure_wtagger_sigma_enlarge_err_el= ( (rrv_sigma_data_el.getError()/rrv_sigma_data_el.getVal())**2 + (rrv_sigma_MC_el.getError()/rrv_sigma_MC_el.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_el

        print "Pure W-tagger SF of el %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el, pure_wtagger_sf_el_err)
        print "Pure W-tagger SF of mu %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu, pure_wtagger_sf_mu_err)
        print "Pure W-tagger mean shift el %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_el, pure_wtagger_mean_shift_err_el)
        print "Pure W-tagger sigma enlarge el %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_el, pure_wtagger_sigma_enlarge_err_el)

        self.boostedW_fitter_el.file_out_ttbar_control.write( "\n***************************************************" )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of el %s     : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el, pure_wtagger_sf_el_err) )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of mu %s        : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu, pure_wtagger_sf_mu_err) )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger mean shift el %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_el, pure_wtagger_mean_shift_err_el) )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger sigma enlarge el %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_el, pure_wtagger_sigma_enlarge_err_el) )

        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_el_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_mu_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_el_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_mu_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.Print();
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.Print();

        rrv_number_total_ttbar_TotalMC_el = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_el");
        rrv_number_total_ttbar_TotalMC_mu = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_mu");
        rrv_number_total_ttbar_data_el = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_el");
        rrv_number_total_ttbar_data_mu = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_mu");

        rrv_number_total_ttbar_TotalMC_el.Print();
        rrv_number_total_ttbar_TotalMC_mu.Print();
        rrv_number_total_ttbar_data_el.Print();
        rrv_number_total_ttbar_data_mu.Print();

        print "el TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal());
        print "mu TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal());
        print "el data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal());
        print "mu data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal());
        
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal()));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal()));

        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nmu TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal()));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nmu data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal()));

        tmp_eff_MC_el_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal();
        tmp_eff_MC_mu_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal();
        tmp_eff_data_el_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal();
        tmp_eff_data_mu_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal();

        tmp_eff_MC_el_extremefail_error =tmp_eff_MC_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_el.getError()/rrv_number_total_ttbar_TotalMC_el.getVal() )**2 );
        tmp_eff_MC_mu_extremefail_error =tmp_eff_MC_mu_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_mu.getError()/rrv_number_total_ttbar_TotalMC_mu.getVal() )**2 );
        tmp_eff_data_el_extremefail_error =tmp_eff_data_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_data_el.getError()/rrv_number_total_ttbar_data_el.getVal() )**2 );
        tmp_eff_data_mu_extremefail_error =tmp_eff_data_mu_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() )**2+ (rrv_number_total_ttbar_data_mu.getError()/rrv_number_total_ttbar_data_mu.getVal() )**2 );

        print "eff_MC_el_extremefail_error %s: %f"%(label,tmp_eff_MC_el_extremefail_error)
        print "eff_MC_mu_extremefail_error %s: %f"%(label,tmp_eff_MC_mu_extremefail_error)
        print "eff_data_el_extremefail_error %s: %f"%(label,tmp_eff_data_el_extremefail_error)
        print "eff_data_mu_extremefail_error %s: %f"%(label,tmp_eff_data_mu_extremefail_error)

        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_MC_el_extremefail_error %s: 0.9%f"%(label,tmp_eff_MC_el_extremefail_error) );
        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_data_el_extremefail_error %s: 0.9%f"%(label,tmp_eff_data_el_extremefail_error) );
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\neff_MC_mu_extremefail_error %s: 0.9%f"%(label,tmp_eff_MC_mu_extremefail_error) );
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\neff_data_mu_extremefail_error %s: 0.9%f"%(label,tmp_eff_data_mu_extremefail_error) );
        
        rrv_eff_MC_el.Print()
        rrv_eff_MC_mu.Print()
        rrv_eff_data_el.Print()
        rrv_eff_data_mu.Print()

        tmp_eff_MC_el_LP =1. - rrv_eff_MC_el.getVal() - tmp_eff_MC_el_extremefail;
        tmp_eff_MC_mu_LP =1. - rrv_eff_MC_mu.getVal() - tmp_eff_MC_mu_extremefail;
        tmp_eff_data_el_LP =1. - rrv_eff_data_el.getVal() - tmp_eff_data_el_extremefail;
        tmp_eff_data_mu_LP =1. - rrv_eff_data_mu.getVal() - tmp_eff_data_mu_extremefail;

        tmp_eff_MC_el_LP_err = TMath.Sqrt( rrv_eff_MC_el.getError()**2 + tmp_eff_MC_el_extremefail_error**2 );
        tmp_eff_MC_mu_LP_err = TMath.Sqrt( rrv_eff_MC_mu.getError()**2 + tmp_eff_MC_mu_extremefail_error**2 );
        tmp_eff_data_el_LP_err = TMath.Sqrt( rrv_eff_data_el.getError()**2 + tmp_eff_data_el_extremefail_error**2 );
        tmp_eff_data_mu_LP_err = TMath.Sqrt( rrv_eff_data_mu.getError()**2 + tmp_eff_data_mu_extremefail_error**2 );

        print "LP Eff of el data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err);
        print "LP Eff of el MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err);
        print "LP Eff of mu data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_mu_LP, tmp_eff_data_mu_LP_err);
        print "LP Eff of mu MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_mu_LP, tmp_eff_MC_mu_LP_err);

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of el data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of el MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nLP Eff of mu data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_mu_LP, tmp_eff_data_mu_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nLP Eff of mu MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_mu_LP, tmp_eff_MC_mu_LP_err));
        
        pure_wtagger_sf_el_LP = tmp_eff_data_el_LP / tmp_eff_MC_el_LP;
        pure_wtagger_sf_mu_LP = tmp_eff_data_mu_LP / tmp_eff_MC_mu_LP;
        pure_wtagger_sf_el_LP_err = pure_wtagger_sf_el_LP*TMath.Sqrt( (tmp_eff_data_el_LP_err/tmp_eff_data_el_LP)**2 + (tmp_eff_MC_el_LP_err/tmp_eff_MC_el_LP)**2 );
        pure_wtagger_sf_mu_LP_err = pure_wtagger_sf_mu_LP*TMath.Sqrt( (tmp_eff_data_mu_LP_err/tmp_eff_data_mu_LP)**2 + (tmp_eff_MC_mu_LP_err/tmp_eff_MC_mu_LP)**2 );

        print "Pure W-tagger LP SF of el %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el_LP, pure_wtagger_sf_el_LP_err)
        print "Pure W-tagger LP SF of mu %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu_LP, pure_wtagger_sf_mu_LP_err)

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nPure W-tagger LP SF of el %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el_LP, pure_wtagger_sf_el_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nPure W-tagger LP SF of mu %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu_LP, pure_wtagger_sf_mu_LP_err));



### function to call single channel fits
def control_sample(channel="mu",isherwig=0):

    print "control sample "+channel;

    if isherwig ==0:
        print " Option herwig 0 --> enter in the single fit of pythia"
        ## create the object and do the single for pythia
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,130)
        boostedW_fitter.fit_TTbar_controlsample(isherwig);

    elif isherwig ==1:
        print " Option herwig 0 --> enter in the single fit of herwig"
        ## create the object and do the single for herwig
        boostedW_fitter_herwig=doFit_wj_and_wlvj(channel,"ggH600",40,130,"_herwig")
        boostedW_fitter_herwig.fit_TTbar_controlsample(isherwig);

    elif isherwig==2:
        print " Option herwig 2 --> enter in the single fit of pythia and herwig"
        ## do Pythia and herwig analysis at the same time
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,130)
        boostedW_fitter.fit_TTbar_controlsample(0,ttbarMC);
        boostedW_fitter_herwig = doFit_wj_and_wlvj(channel,"ggH600",40,130,"_herwig")
        boostedW_fitter_herwig.fit_TTbar_controlsample(1);


### function to call simultaneous channel fits
def control_sample_simultaneous():

    print "control_sample_simultaneous";

    if options.herwig == 0 :
        print " Option herwig 0 --> enter in the single simultaneous fit "
        boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous(options.herwig)
    elif options.herwig ==1 :
        print " Option herwig 1 --> enter in the single simultaneous fit for herwig sample"
        boostedW_fitter_sim_herwig = doFit_wj_and_wlvj_simultaneous(options.herwig)
    elif isherwig==2:
        print " Option herwig 2 --> enter in the single simultaneous fit of pythia and herwig"
        boostedW_fitter_sim=doFit_wj_and_wlvj_simultaneous()
        boostedW_fitter_sim=doFit_wj_and_wlvj_simultaneous(1)

### main code
if __name__ == '__main__':

    channel=options.channel;

    if options.fitwtaggersim:
        print 'fitwtagger for el+mu sample'
        control_sample_simultaneous();

    elif options.fitwtagger:
        print 'fitwtagger for %s sample'%(channel)
        control_sample(channel,options.herwig);
 
        


