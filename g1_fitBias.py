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

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit,RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet,RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite, RooMCStudy, RooGlobalFunc,RooChi2MCSModule, RooCurve

############################################
#              Job steering                #
############################################

parser = OptionParser()
parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-j','--njets',    help='number of jets: 1 or 2 ,default:single' , type=int, default = 1) 
parser.add_option('-p','--category', help='purity category: LP (low purity) or HP (high purity) or NP (no purity selection), default:HP',type="string", default = "HP")
parser.add_option('-l','--channel',  help='lepton flavor: el or mu  or both , default:both' ,type="string", default ="mu" ) ## lepton flavour
parser.add_option('-f','--inPath',   help='directory with workspace' , default = "./" )
parser.add_option('-m','--mass',     help='test signal yield for this mass', type=int, default=-1)
parser.add_option('-n','--nexp',     help='number of toys', type=int, default=1000)
parser.add_option('-g','--fgen',     help='function to generate toys Exp,ExpTail,Pow2,ExpN)', type="string", default="ExpN")
parser.add_option('-r','--fres',     help='function to fit toys (Exp,ExpTail,Pow2,ExpN)',     type="string", default="ExpN")
parser.add_option('-s','--storeplot',help='in case of more than 10 toys just 1/3 stored, more than 100 1/10',     type=int, default=0)
parser.add_option('-z','--skipMC',     help='options to skip pure mc w+jets toys',     type=int, default=0)

(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")


from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf


class doBiasStudy_mlvj:

    def __init__(self,in_channel,in_signal_sample,in_mlvj_min=700., in_mlvj_max=3000., in_mj_min=40, in_mj_max=140, generation_model="ExpN", fit_model="ExpN", input_workspace=None):

        self.setTDRStyle();

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        print "###################### construnctor ############################# ";
        ### set the channel type --> electron or muon
        self.channel=in_channel;

        ## event categorization as a function of the purity and the applied selection
        self.wtagger_label = options.category;

        self.BinWidth_mj=5.;
        nbins_mj=int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max=in_mj_min+nbins_mj*self.BinWidth_mj;

        self.BinWidth_mlvj=100.;
        self.in_mlvj_min = in_mlvj_min;
        self.nbins_mlvj=int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        self.in_mlvj_max=in_mlvj_min+self.nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
        rrv_mass_j = RooRealVar("rrv_mass_j","pruned m_{J}",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV/c^{2}");
        rrv_mass_j.setBins(nbins_mj);

        self.mj_sideband_lo_min = in_mj_min;
        self.mj_sideband_lo_max = 65;
        self.mj_signal_min = 65;
        self.mj_signal_max = 105;
        self.mj_sideband_hi_min = 105;
        self.mj_sideband_hi_max = in_mj_max;
                                                                                            
        ### define invariant mass WW variable
        rrv_mass_lvj= RooRealVar("rrv_mass_lvj","m_{WW}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV/c^{2}");
        rrv_mass_lvj.setBins(self.nbins_mlvj);

        ### define shapes to be used to generate and fit
        self.generation_model =generation_model ;
        self.fit_model = fit_model ;

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4bias_ = RooWorkspace("workspace4bias_%s_%s_%s_%s"%(self.channel,self.wtagger_label,self.generation_model,self.fit_model),"workspace4bias_%s_%s_%s_%s"%(self.channel,self.wtagger_label,self.generation_model,self.fit_model));
        else:
            self.workspace4bias_ = input_workspace;
            
        getattr(self.workspace4bias_,"import")(rrv_mass_lvj);
        getattr(self.workspace4bias_,"import")(rrv_mass_j);

        ## zone definition in the jet mass
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

                                                                           
        ## set the signal sample
        self.file_Directory = "AnaSigTree_new/";
        
        self.signal_sample = in_signal_sample;

        self.file_data = ("treeEDBR_data_xww.root");#keep blind!!!!

        self.file_pseudodata = ("treeEDBR_allBkg_xww.root");#fake data

        self.file_signal = ("treeEDBR_%s_xww.root"%(self.signal_sample));

        self.file_WJets0_mc = ("treeEDBR_WJetsPt180_xww.root");

        self.file_VV_mc = ("treeEDBR_VV_xww.root");# WW+WZ

        self.file_TTbar_mc = ("treeEDBR_TTBARpowheg_xww.root");

        self.file_STop_mc = ("treeEDBR_SingleTop_xww.root");

        ## event categorization as a function of the label        
        if self.wtagger_label=="HP" :
            if self.channel=="el":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
            if self.channel=="mu":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;

        if self.wtagger_label=="LP":
            self.wtagger_cut=0.75 ;
            self.wtagger_cut_min=0.5 ;

        if self.wtagger_label=="NP":
            self.wtagger_cut=10000;

        self.categoryID=-1;
        if self.wtagger_label=="LP" and self.channel=="el": self.categoryID=0;
        if self.wtagger_label=="HP" and self.channel=="el": self.categoryID=1;
        if self.wtagger_label=="LP" and self.channel=="mu": self.categoryID=2;
        if self.wtagger_label=="HP" and self.channel=="mu": self.categoryID=3;

        ## color palette
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : 2,
            'VV' : 4,
            'STop' : 7,
            'TTbar' : 210,
            'ggH' : 1,
            'vbfH' : 12,
            'Signal': 1,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        ## for basic selection
        self.vpt_cut   = 200;
        self.pfMET_cut = 50;
        self.lpt_cut   = 50;
        if self.channel=="el":
            self.pfMET_cut= 80; self.lpt_cut = 90;#very tight
        self.deltaPhi_METj_cut =2.0;

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
        self.tdrStyle.SetOptStat(1111); #To display the mean and RMS:
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
    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):

        if TString(mass_spectrum).Contains("_mj"):   rrv_x = self.workspace4bias_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4bias_.var("rrv_mass_lvj");

        if in_model_name == "CB_v1":
            label_tstring=TString(label);
            if label_tstring.Contains("H600"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,600,580,620);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,67,40,80);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,-1,-2,-0.5);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,10,80 );
            elif label_tstring.Contains("H700"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,700,650,750);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,100,40,150);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,-1,-3,-0.1);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,10,40);
            elif label_tstring.Contains("ggH800"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,780,700,850);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,140,120,160);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,-1,-4,0);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,5 , 2, 7);
            elif label_tstring.Contains("vbfH800"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,800,750,850);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,140,120,160);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,-1,-4,0);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel,5 , 2, 7);
            elif label_tstring.Contains("ggH900"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,880,820,950);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,170,140,200);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,1,0,4);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel, 2., 0.5,5);
            elif label_tstring.Contains("vbfH900"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,900,880,920);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,170,140,200);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,1);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel, 2., 0.5,5);
            elif label_tstring.Contains("ggH1000"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,920,800,1150);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,200,100,300);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,1,0.1,3);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel,2.,0.5,4);
            elif label_tstring.Contains("vbfH1000"):
                rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel,"rrv_mean_CB"+label+"_"+self.channel,1000,980,1150);
                rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel,"rrv_sigma_CB"+label+"_"+self.channel,200,100,300);
                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel,"rrv_alpha_CB"+label+"_"+self.channel,0.72);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel,"rrv_n_CB"+label+"_"+self.channel,2.,0.5,4);
            else:
                if label_tstring.Contains("600") and (not label_tstring.Contains("1600") ):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 600, 550, 650);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);

                elif label_tstring.Contains("700") and (not label_tstring.Contains("1700") ):
                     rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 700, 600, 800);
                     rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);
                     
                elif label_tstring.Contains("800") and (not label_tstring.Contains("1800") ):
                     rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 800, 600, 800);
                     rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 40,10 ,90);
                     
                elif label_tstring.Contains("900") and (not label_tstring.Contains("1900") ):
                    rrv_mean_CB=RooRealVaDr("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 900, 600, 800);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 40,10 ,90);

                elif label_tstring.Contains("1000"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1000, 900,1100);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1100"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1100,1000,1200);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1200"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1200,1100,1300);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1300"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1300,1200,1400);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1400"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1400,1300,1500);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1500"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1500,1400,1600);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                    
                elif label_tstring.Contains("1600"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1600,1500,1700);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1700"):
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1700,1500,1800);

                elif label_tstring.Contains("1800"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1800,1500,1900);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1900"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1900,1500,2000);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2000"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2000,1800,2200);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2100"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2100,1800,2300);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2200"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2200,1800,2400);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2300"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2300,1800,2500);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2400"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2400,1800,2600);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2500"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2500,2000,2700);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                else :
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,700,550,2500);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,4,1,5);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,10,40);

            model_pdf = RooCBShape("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);



        ## Crystal ball shape for Bulk GR samples and higgs
        if in_model_name == "DoubleCB_v1":
            label_tstring=TString(label);
            print "########### Double CB for Bulk graviton mlvj ############"

            if label_tstring.Contains("M600") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M600_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 600, 550, 650);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);
                         
            elif label_tstring.Contains("M700") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M700_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 700, 600, 800);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);
                                          
            elif label_tstring.Contains("M800") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M800_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,820,790,880);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,50,40,70);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 15.,5.,25.);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.64,1.,1.9);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,15.,5.,25.);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,1.,1.9);

                                          
            elif label_tstring.Contains("M900") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M900_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,920,850,950);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,59,45,70);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 25.,2,45);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,25.,0.1,45);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.25,0.5,3.);
                      
            elif label_tstring.Contains("M1000") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1000_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1020,970,1070);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,55,40,65);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,45);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.4,0.5,3.5);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.,0.5,3.5);
                        
            elif label_tstring.Contains("M1100") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1100_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1120,1080,1150);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,65,55,75);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,25);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,25);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);

            elif label_tstring.Contains("M1200") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1200_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1220,1200,1250);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,65,55,75);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,30);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,5.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,30);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,5.);
         
            elif label_tstring.Contains("M1300") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1300_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1320,1300,1350);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,70,60,75);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.3,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.3,0.5,3.);
                    
                         
            elif label_tstring.Contains("M1400") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1400_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1420,1400,1440);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,77,65,85);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.5);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.5);

            elif label_tstring.Contains("M1500") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1500_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1515,1500,1530);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,81,71,91);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 15.,0.01,25);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.5);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,15.,0.01,25);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.5);

            elif label_tstring.Contains("M1600") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1600_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1620,1600,1640);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,81,70,90);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                         
            elif label_tstring.Contains("M1700") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1700_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1720,1700,1740);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,90,75,96);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                                            
            elif label_tstring.Contains("M1800") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1800_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1820,1800,1840);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,90,75,100);
                    rrv_n1_CB    = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                      
            elif label_tstring.Contains("M1900") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1900_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1920,1900,1940);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,95,80,115);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                       
            elif label_tstring.Contains("M2000") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2000_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2020,2000,2040);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,100,80,115);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                       
            elif label_tstring.Contains("M2100") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2100_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2120,2100,2140);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,105,85,115);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                      
            elif label_tstring.Contains("M2200") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2200_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2220,2200,2250);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,115,75,140);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                      
            elif label_tstring.Contains("M2300") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2300_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2320,2300,2340);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,115,95,120);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 15.,0.2,30);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,15.,0.2,20);                        
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                        
            elif label_tstring.Contains("M2400") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2400_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2420,2400,2440);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,115,100,125);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35)                        
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                       
            elif label_tstring.Contains("M2500") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2500_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2520,2500,2540);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,125,90,145);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);

            elif label_tstring.Contains("M1000_W150") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1000,500,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 150,50 ,500);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M1000_W300") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1000,500,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 300,50 ,800);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M1000_W50") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1000,500,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,10 ,200);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M1500_W175") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1500,1000,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 75,50 ,250);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M1500_W225") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1500,1000,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 225,150 ,450);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M1500_W450") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1500,1000,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 450,400 ,700);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M2100_W105") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2100,1500,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 105,90 ,300);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M2100_W315") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2100,1500,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 315,250 ,600);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("M2100_W630") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2100,1500,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 630,550 ,900);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            else :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,700,550,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                    rrv_n1_CB = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);


            model_pdf = ROOT.RooDoubleCB("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB);

        ## ExpN pdf for W+jets bkg fit
        if in_model_name == "ExpN":
            
            print "########### ExpN funtion for W+jets mlvj ############"
            rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-3e-3,-1e-1,-1e-5);
            if(ismc==1):
              rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 1e3, -1e2, 1e4);
            else :
              if self.channel == "el" :
                 rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 1e3, -1e2, 1e4);
              elif self.wtagger_label == "LP" :
                 rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 1e3, -1e2, 1e4);
              else:
                 rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 5e2, 0, 1e3);

            model_pdf = ROOT.RooExpNPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ExpN, rrv_n_ExpN);

        ## levelled exp for W+jets bkg fit
        if in_model_name == "ExpTail":

            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if self.wtagger_label == "LP":
                rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
                rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1e-1,-1.e2,1e6);
            else:
             if self.channel == "el" :
               if ismc == 1 and label_tstring.Contains("sb_lo"):
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 139,0.,355);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2e-2,-1.e-2,5.5e-2);
               elif ismc == 1 and label_tstring.Contains("signal_region"):
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 162,18,395);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1.6e-2,-1.e-2,5.5e-2);
               elif ismc == 0 :
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,70,240);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1.e-2,1.3e-1);
             if self.channel == "mu" :
               if ismc == 1 and label_tstring.Contains("sb_lo"):
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 99,10,255);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-1e-2,7.5e-2);
               elif ismc == 1 and label_tstring.Contains("signal_region"):
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 110,20,242);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1e-2,7.5e-2);
               elif ismc == 0 :
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,40,280);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);

            model_pdf = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

   
        ## sum of two exponential
        if in_model_name == "2Exp":

            print "########### 2Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c0_2Exp   = RooRealVar("rrv_c0_2Exp"+label+"_"+self.channel,"rrv_c0_2Exp"+label+"_"+self.channel, -5e-3, -8e-3,-4e-3);
            rrv_c1_2Exp   = RooRealVar("rrv_c1_2Exp"+label+"_"+self.channel,"rrv_c1_2Exp"+label+"_"+self.channel, -1e-3, -4e-3,-1e-4);
            rrv_frac_2Exp = RooRealVar("rrv_frac_2Exp"+label+"_"+self.channel,"rrv_frac_2Exp"+label+"_"+self.channel, 0., 0., 1e-2);
            model_pdf = ROOT.Roo2ExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0_2Exp,rrv_c1_2Exp,rrv_frac_2Exp);

        ## sum of two exponential
        if in_model_name == "Exp" or in_model_name == "Exp_sr":
            print "########### Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,-0.05,-0.1,0.);
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp);


 
        ## For mlvj fit -> Pow function can replace exp
        if in_model_name == "Pow2":
            print "########### Pow2 Pdf for mlvj fit ############"
            rrv_c0 = RooRealVar("rrv_c0_Pow2"+label+"_"+self.channel,"rrv_c0_Pow2"+label+"_"+self.channel, 0, -20., 20);
            rrv_c1 = RooRealVar("rrv_c1_Pow2"+label+"_"+self.channel,"rrv_c1_Pow2"+label+"_"+self.channel, 0, -5, 5);
            model_pdf = ROOT.RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"TMath::Power(@0,@1 + @0*@2)",RooArgList(rrv_x,rrv_c0,rrv_c1));
                                           
        ## return the pdf
        getattr(self.workspace4bias_,"import")(model_pdf)
        return self.workspace4bias_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)

    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500):

      ##### define an extended pdf from a standard Roofit One
      print " "
      print "###############################################"
      print "## Make model : ",label," ",in_model_name,"##";
      print "###############################################"
      print " "

      rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum,"rrv_number"+label+"_"+self.channel+mass_spectrum,500.,0.,1e5);
      ## call the make RooAbsPdf method
      model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList,ismc_wjet)
      print "######## Model Pdf ########"
      model_pdf.Print();
      
      ## create the extended pdf
      model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
      print "######## Model Extended Pdf ########"

      #### put all the parameters ant the shape in the workspace
      getattr(self.workspace4bias_,"import")(rrv_number)
      getattr(self.workspace4bias_,"import")(model)
      self.workspace4bias_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();

      ## return the total extended pdf
      return self.workspace4bias_.pdf("model"+label+"_"+self.channel+mass_spectrum);
  
    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self,label, mlvj_region):
        print "model"+label+mlvj_region+"_"+self.channel+"_mlvj"
        return self.workspace4bias_.pdf("model"+label+mlvj_region+"_"+self.channel+"_mlvj");

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region="_signal_region"):
        print "########### Fixing amd return a general mlvj model ############"
        rdataset_General_mlvj = self.workspace4bias_.data("rdataset4bias%s%s_%s_mlvj"%(label, mlvj_region,self.channel))
        model_General = self.get_mlvj_Model(label,mlvj_region);
        rdataset_General_mlvj.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_mlvj);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
          param.setConstant(kTRUE);
          param.Print();
          param=par.Next()
        return self.get_mlvj_Model(label,mlvj_region);


    ###### get TTbar model mlvj in a region
    def get_TTbar_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing TTbar mlvj model for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_TTbar",mlvj_region);

    ###### get Single Top model mlvj in a region
    def get_STop_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing Stop mlvj model for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_STop",mlvj_region);

    ###### get Signal model mlvj in a region
    def get_signal_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing signal mlvj model for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_%s"%(self.signal_sample),mlvj_region);

    ###### get VV mlvj in a region
    def get_VV_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing VV mlvj for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_VV",mlvj_region);
                                                                                                                    

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj",additional="",notExtended=0):
        print "########### Fixing an Extended Pdf for mlvj ############"
        rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
        if notExtended:
           label = "_pdf"+label;    
        model = self.get_mlvj_Model(label,mlvj_region+additional);
        rdataset.Print();
        model.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); 
        par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()

    ### clone Model Extend in a not extend skipping the normalization
    def clone_Model(self, inputPdf, label, mlvj_region="_signal_region",mass_spectrum="_mlvj"):
        print "########### Cloning an Extended Pdf for mlvj ############"
        rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum));
        model = self.get_mlvj_Model(label,mlvj_region);                
        rdataset.Print();
        inputPdf.Print();
        model.Print()
        parameters  = inputPdf.getParameters(rdataset);
        parameters2 = model.getParameters(rdataset);
        par         = parameters.createIterator();
        par2        = parameters2.createIterator();
        par.Reset();
        par2.Reset();
        param  = par.Next();
        param2 = par2.Next();

        if parameters.getSize() < parameters2.getSize() :
         while (param):
             if(not TString(param2.GetName()).Contains("number")):
                 param.setVal(param2.getVal());
                 param.setError(param2.getError());
                 param.Print(); param2.Print();
                 param  = par.Next();
                 param2 = par2.Next();
             else:
                 param2 = par2.Next();
        elif  parameters.getSize() > parameters2.getSize() :
         while (param2):
             if(not TString(param.GetName()).Contains("number")):
                 param.Print(); param2.Print();
                 param.Print(); param2.Print();
                 param.setVal(param2.getVal());
                 param.setError(param2.getError());
                 param  = par.Next();
                 param2 = par2.Next();
             else:
                param = par.Next();

        else:
         while (param):
                 param.Print(); param2.Print(); 
                 param.setVal(param2.getVal());
                 param.setError(param2.getError());
                 param  = par.Next();
                 param2 = par2.Next();
            
    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Deco_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj ############"
        rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum));
        model = self.workspace4bias_.pdf("model_pdf"+label+mlvj_region+"_"+self.channel+"_mlvj"+"_Deco"+label+mlvj_region+"_"+self.channel+"_"+self.wtagger_label+"_mlvj");
        rdataset.Print();
        model.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); 
        par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_mlvj_model_single_MC(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",in_file_name," ",label," ",mlvj_model," ",in_range," ##################"
        ## imporparam_generatedt variable and dataset
        rrv_mass_lvj = self.workspace4bias_.var("rrv_mass_lvj")
        rdataset     = self.workspace4bias_.data("rdataset4bias"+label+in_range+"_"+self.channel+"_mlvj");
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,mlvj_model,"_mlvj",constrainslist,ismc);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4bias_,"import")(rfresult)

        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(int(rrv_mass_lvj.getBins())));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,2,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull

        ## get the pull
        mplot_pull = self.get_pull(rrv_mass_lvj,mplot);
        parameters_list = model.getParameters(rdataset);
        
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres),in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, logy);

         
        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4bias_.pdf("model_pdf%s%s_%s_mlvj"%(label,in_range,self.channel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.channel+"_mlvj");
            Deco = PdfDiagonalizer("Deco"+label+in_range+"_"+self.channel+"_"+self.wtagger_label+"_mlvj",wsfit_tmp,rfresult_pdf); ## in order to have a good name
            model_pdf_deco = Deco.diagonalize(model_pdf); ## diagonalize
            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4bias_,"import")(model_pdf_deco);
            wsfit_tmp.Print("v");
            model_pdf_deco.getParameters(rdataset).Print("v");
            model_pdf.getParameters(rdataset).Print("v");
            model_pdf.Print();
            model_pdf_deco.Print();
            
            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame( RooFit.Bins(int(rrv_mass_lvj.getBins())));
            
            rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

            mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

            rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
            rrv_number_dataset.setError(0.)
            draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            leg = self.legend4Plot(mplot_deco,0); ## add the legend
            mplot_deco.addObject(leg);

            self.draw_canvas( mplot_deco, "plots_%s_%s_%s_g1/other_%s_%s/"%(options.additioninformation, self.channel, self.wtagger_label, options.fgen,options.fres), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco",0,logy);

     ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0):

        print "############### Fit mlvj in mj sideband: ",label," ",mlvj_region," ",mlvj_model," ##################"

        rrv_mass_lvj = self.workspace4bias_.var("rrv_mass_lvj");
        rdataset_data_mlvj = self.workspace4bias_.data("rdataset4bias_data%s_%s_mlvj"%(mlvj_region,self.channel))

        ## get and fix the minor component shapes in the sb low
        model_VV_backgrounds    = self.get_VV_mlvj_Model("_sb_lo");
        number_VV_sb_lo_mlvj    = self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)); ## get the normalization
        model_TTbar_backgrounds = self.get_TTbar_mlvj_Model("_sb_lo");
        number_TTbar_sb_lo_mlvj = self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)); ## get the normalization
        model_STop_backgrounds  = self.get_STop_mlvj_Model("_sb_lo");
        number_STop_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)); ## get the normalization

        self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)).Print();

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_sb_lo_from_fitting"%(label), mlvj_model,"_mlvj");
        model_pdf_WJets.Print();
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_sb_lo = self.workspace4bias_.var("rrv_number%s_sb_lo_%s_mlvj"%(label,self.channel)).clone("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)); 
        model_WJets= RooExtendPdf("model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),"model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),model_pdf_WJets,number_WJets_sb_lo);
        number_WJets_sb_lo.Print()

        ## Add the other bkg component fixed to the total model --> in the extended way
        model_data = RooAddPdf("model_data%s%s_%s_mlvj"%(label,mlvj_region,self.channel),"model_data%s%s_%s_mlvj"%(label,mlvj_region,self.channel),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds));
        
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        getattr(self.workspace4bias_,"import")(model_data)

        model_WJets.getParameters(rdataset_data_mlvj).Print("v");
        self.workspace4bias_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label,self.channel)).getParameters(rdataset_data_mlvj).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_sb_lo_%s_mlvj"%(self.channel),"rrv_number_data_sb_lo_%s_mlvj"%(self.channel),
                                                 self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getVal() );

        rrv_number_data_sb_lo_mlvj.setError( TMath.Sqrt(self.workspace4bias_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)).getError()));

        getattr(self.workspace4bias_,"import")(rrv_number_data_sb_lo_mlvj)


        ### plot for WJets default + default shape
        if label=="_WJets0":

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins())));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_sb_lo_%s_mlvj,model_STop_sb_lo_%s_mlvj,model_VV_sb_lo_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_sb_lo_%s_mlvj,model_STop_sb_lo_%s_mlvj,model_VV_sb_lo_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_sb_lo_%s_mlvj,model_STop_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_STop_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());

            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_sb_lo_%s_mlvj,model_STop_sb_lo_%s_mlvj,model_VV_sb_lo_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_sb_lo_%s_mlvj,model_STop_sb_lo_%s_mlvj,model_VV_sb_lo_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_sb_lo_%s_mlvj,model_STop_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_STop_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());
 

            rdataset_data_mlvj.plotOn(mplot,RooFit.Name("data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            ### draw the error band
            draw_error_band(rdataset_data_mlvj, model_data,self.workspace4bias_.var("rrv_number_data_sb_lo_%s_mlvj"%(self.channel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible());
            rdataset_data_mlvj.plotOn(mplot,RooFit.Name("data_invisible1"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

            ### Add the legend to the plot
            leg=self.legend4Plot(mplot,0,1,0., 0.06, 0.16, 0.);
            mplot.addObject(leg)

            ### calculate the chi2
            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof = nBinX-self.nPar_float_in_fitTo;
            print mplot.chiSquare();
            print "#################### nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_mass_lvj,mplot);
            parameters_list = model_data.getParameters(rdataset_data_mlvj);
                
            self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres), "m_lvj_sb_lo%s"%(label),"",1,1)

        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sb_lo_from_fitting_mlvj"%(label));
        Deco = PdfDiagonalizer("Deco%s_sb_lo_from_fitting_%s_%s_mlvj"%(label,self.channel,self.wtagger_label),wsfit_tmp,rfresult);
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets);
        model_pdf_WJets_deco.Print("v");
        model_pdf_WJets_deco.getParameters(rdataset_data_mlvj).Print("");
        wsfit_tmp.allVars().Print("v");
        getattr(self.workspace4bias_,"import")(model_pdf_WJets_deco);



    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass ):

        print "################### get_mj_and_mlvj_dataset : ",in_file_name," ",label," ################## ";

        fileIn_name = TString(options.inPath+"/"+self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("SelectedCandidatesPlain");
        
        rrv_mass_j = self.workspace4bias_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4bias_.var("rrv_mass_lvj")
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ##### dataset of m_lvj -> scaled and not scaled to lumi in different region
        rdataset4bias_sb_lo_mlvj = RooDataSet("rdataset4bias"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4bias_signal_region_mlvj = RooDataSet("rdataset4bias"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4bias_sb_hi_mlvj = RooDataSet("rdataset4bias"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### categorize the event in sideband and signal region --> combined dataset
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");
        combData4bias  = RooDataSet("combData4bias"+label+"_"+self.channel,"combData4bias"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        print "###### N entries: ", treeIn.GetEntries()
        hnum_4region=TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_2region=TH1D("hnum_2region"+label+"_"+self.channel,"hnum_2region"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total

        if self.channel=="el":
              tmp_lumi=19531.85;
        else: tmp_lumi=19538.85;

        for i in range(treeIn.GetEntries()):
            if i % 100000 == 0: print "iEntry: ",i
            treeIn.GetEntry(i);

            if i==0:
                tmp_scale_to_lumi=treeIn.LumiWeight*tmp_lumi;
    
            tmp_jet_mass=getattr(treeIn, jet_mass);
            
            ## event in the whole range
            if treeIn.categories==self.categoryID and treeIn.mZZ> rrv_mass_lvj.getMin() and treeIn.mZZ<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() :

                ### weigh MC events
                tmp_event_weight4bias = treeIn.weight*tmp_lumi;

                #### wtagger_eff_reweight
                if label=="_data":
                   tmp_event_weight4bias=1.;
                
                rrv_mass_lvj.setVal(treeIn.mZZ);
                if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
                    rdataset4bias_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );

                    data_category.setLabel("sideband");
                    combData4bias.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4bias);

                if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                    rdataset4bias_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );

                    data_category.setLabel("signal_region");
                    combData4bias.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4bias);
                    hnum_2region.Fill(1,tmp_event_weight4bias);

                if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                    rdataset4bias_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );

                if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max:
                    hnum_4region.Fill(-1,tmp_event_weight4bias);
                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                    hnum_4region.Fill(0,tmp_event_weight4bias);
                if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max:
                    hnum_4region.Fill(1,tmp_event_weight4bias);

                hnum_4region.Fill(2,tmp_event_weight4bias);


        ### scaler to lumi for MC in 4bias datasets
        rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi)
        rrv_scale_to_lumi.Print()
        getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi)

        ### prepare m_lvj dataset to be compared with the fit results
        rrv_number_dataset_AllRange_mlvj = RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(2));
        
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_AllRange_mlvj)

        ### import the dataser
        getattr(self.workspace4bias_,"import")(rdataset4bias_sb_lo_mlvj);
        getattr(self.workspace4bias_,"import")(rdataset4bias_signal_region_mlvj);
        getattr(self.workspace4bias_,"import")(rdataset4bias_sb_hi_mlvj);
        getattr(self.workspace4bias_,"import")(combData4bias);


        ### prepare m_j dataset
        rrv_number_dataset_sb_lo_mj=RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj=RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));
        rrv_number_dataset_sb_hi_mj=RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_sb_lo_mj)
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_signal_region_mj)
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_sb_hi_mj)
                
        #### print everything        
        rdataset4bias_sb_lo_mlvj.Print();
        rdataset4bias_signal_region_mlvj.Print();
        rdataset4bias_sb_hi_mlvj.Print();
        rrv_number_dataset_AllRange_mlvj.Print()
        rrv_number_dataset_sb_lo_mj.Print()
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_sb_hi_mj.Print()

    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist();
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
          if(y == 0):
           hpull.SetPoint(ipoint,x,10)
       
        mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins())));
        medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
        mplot_pull.addObject(medianLine);
        mplot_pull.addPlotable(hpull,"P");
        mplot_pull.SetTitle("");
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
      print "############### draw the banner ########################"

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
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1.):
        print "############### draw the legend ########################"
        if left==-1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetLineColor(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.41+x_offset_low, 0.61+y_offset_low, 0.76+x_offset_high, 0.93+y_offset_high, "", "NDC");
            theLeg.SetName("theLegend");
            if TwoCoulum :
                theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.040);

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";

        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            print objName;
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data()=="STop" : theLeg.AddEntry(theObj, "Single Top","F");
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ/ZZ","F");
                    elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F");
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
                    elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
                    elif TString(objName).Contains("Bulk"):
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M600") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=0.6TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M700") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=0.7TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M800") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M800"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=0.8TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M900") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=0.9TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1000") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1100") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.1TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1200") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.2TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1300") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.3TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1400") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.4TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1500") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.5TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1600") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.6TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1700") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.7TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1800") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1800"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.8TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1900") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=1.9TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2000") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=2TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M3200") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=2.1TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2200") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=2.2TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2300") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=2.3TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2400") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=2.4TeV #tilde{k}=0.2 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2500") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "BulkG M=2.5TeV #tilde{k}=0.2 (#times100)";
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
           theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;

    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, mplot, mplot_pull,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0):# mplot + pull

        print "############### draw the canvas with pull ########################"
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
        Directory = TString(in_directory+self.signal_sample);
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
            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*100);
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1);

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0, frompull=0):

        print "############### draw the canvas without pull ########################"
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);

        if frompull and logy :
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/100)
        elif not frompull and logy :
            in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())
            

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
        
        Directory=TString(in_directory+self.signal_sample);
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
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()*100);
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());
       

    ##### Get Lumi for banner title
    def GetLumi(self):
        if self.channel=="el": return 19.5;
        if self.channel=="mu": return 19.5;


    def biasAnalysis(self):
   
     print"######################## begin the bias analysis ###########################";  
     ## get the signal and fit it
     self.get_mj_and_mlvj_dataset(self.file_signal,"_%s"%(self.signal_sample), "mJJNoKinFit")# to get the shape of m_lvj
     self.fit_mlvj_model_single_MC(self.file_signal,"_%s"%(self.signal_sample),"_signal_region","DoubleCB_v1", 0, 0, 0, 0);

     ## get diboson and fit it
     self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "mJJNoKinFit");
     self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV","_sb_lo","Exp", 0, 0, 1, 0);

     ## get SingleTop and fit it
     self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop", "mJJNoKinFit");
     self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop","_sb_lo","Exp", 0, 0, 1);

     ## get TTbar and fit it
     self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar", "mJJNoKinFit");# to get the shape of m_lvj
     self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar","_sb_lo","Exp");
                 
     ## get Wjets and fit it in the sb
     self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0", "mJJNoKinFit");# to get the shape of m_lvj          
     self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo",options.fgen, 0, 0, 1, 1);

     ## get data in sb and fit it
     self.get_mj_and_mlvj_dataset(self.file_data,"_data", "mJJNoKinFit"); ## global fit of data in the sidand fixing non dominant bkg
     self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo",options.fgen,1);


     ## fix signal and bkg models that are going to be used in the generation
     self.fix_Model("_%s"%self.signal_sample,"_signal_region","_mlvj");
     self.fix_Model("_WJets0","_sb_lo","_mlvj");
     self.fix_Model("_TTbar","_sb_lo","_mlvj");
     self.fix_Model("_STop","_sb_lo","_mlvj");
     self.fix_Model("_VV","_sb_lo","_mlvj") ;
     self.fix_Model("_WJets0","_sb_lo","_mlvj","_from_fitting");

     ### clone the signal shape --> parameter already fixed
     fitted_signal      = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%(self.signal_sample,"signal_region",self.channel,"mlvj"));

     ### make the pdf for signal not in the extend way, cloning the parameter from the fitted value
     constrainslist_signal_wjet = [];
     model_signal = self.make_Pdf("_%s_signal_region_fit"%(self.signal_sample),"DoubleCB_v1","_mlvj",constrainslist_signal_wjet,1);
     model_signal.Print();     
     self.clone_Model(model_signal,"_%s"%self.signal_sample,"_signal_region","_mlvj");
     self.fix_Model("_%s"%self.signal_sample,"_signal_region","_mlvj","_fit",1);
      
     rrv_number_signal_signal_fit_mc = RooRealVar("rrv_number_signal_region_fit_mc","rrv_number_signal_region_fit_mc",0,-self.workspace4bias_.var("rrv_number_"+self.signal_sample+"_signal_region_"+self.channel+"_mlvj").getVal()*300.,self.workspace4bias_.var("rrv_number_"+self.signal_sample+"_signal_region_"+self.channel+"_mlvj").getVal()*300);
     rrv_number_signal_signal_fit_mc.setVal(self.workspace4bias_.var("rrv_number_"+self.signal_sample+"_signal_region_"+self.channel+"_mlvj").getVal());
     rrv_number_signal_signal_fit_mc.setError(self.workspace4bias_.var("rrv_number_"+self.signal_sample+"_signal_region_"+self.channel+"_mlvj").getError());
     rrv_number_signal_signal_fit_mc.Print();
 
     modified_signal_model_mc = RooExtendPdf("model_%s_signal_region_fit_%s_mlvj"%(self.signal_sample,self.channel),
                                             "model_%s_signal_region_fit_%s_mlvj"%(self.signal_sample,self.channel),model_signal,rrv_number_signal_signal_fit_mc);
     modified_signal_model_mc.Print();
     getattr(self.workspace4bias_,"import")(modified_signal_model_mc);

     ############### Make the MC analysis --> make the Entended pdf for the bkg
     if options.skipMC == 0 :
      constrainslist_bkg_wjet = [];
      model_bkg_wjet    = self.make_Model("_WJets0_sb_lo_fit",options.fres,"_mlvj",constrainslist_bkg_wjet,1);
      model_bkg_wjet.Print();


      ##### Total model for MC
      model_Total_mc    = RooAddPdf("model_Total_background_mc","model_Total_background_mc",RooArgList(modified_signal_model_mc,model_bkg_wjet));

      model_Total_mc.Print();
      getattr(self.workspace4bias_,"import")(model_Total_mc);

      ##### generate models  --> the one fixed and already fitted
      generation_model_wjet = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%("WJets0","sb_lo",self.channel,"mlvj"));
      generation_model_wjet.Print();

      self.workspace4bias_.Print();             

      ## variable and RooMC study
      numevents_mc   = self.workspace4bias_.data("rdataset4bias"+"_WJets0"+"_sb_lo"+"_"+self.channel+"_mlvj").sumEntries();
      print"########  numevents mc ",numevents_mc;


      mc_wjet = RooMCStudy(generation_model_wjet,
                           RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                           RooFit.FitModel(model_Total_mc),
                           RooFit.FitOptions(RooFit.Save(kTRUE),RooFit.SumW2Error(kTRUE),RooFit.Minimizer("Minuit2"),RooFit.Extended(kTRUE)),
                           RooFit.Extended(kTRUE),
                           RooFit.Silence());


      mc_wjet.Print();
     
      ## create module for chi2 evaluation
      chi2_mc = RooChi2MCSModule();
      mc_wjet.addModule(chi2_mc);

      ## generate and fit storing the generated distribution for each toy
      mc_wjet.generateAndFit(options.nexp,int(numevents_mc),1);

      generatedData_wjet         = []; ## distribution of generated toy according to bkg only hypothesis
      fittedPdf_wjet             = []; ## fitted pdf signal + bkg
      fitResults_wjet            = [];
      parameterHisto_wjet        = []; ## histo of the parameters of the fitted pdf
      parameterHistoError_wjet   = []; ## errror on the parameters
      parameterHistoPull_wjet    = []; ## pull wrt the generated one
      chi2distribution_wjet      = [];
      nLLdistribution_wjet       = [];
             
      for iToy in range(options.nexp): ## loop on the toy
 
         if not mc_wjet.genData(iToy) : continue ;
         generatedData_wjet.append(mc_wjet.genData(iToy)); ## take the generated dataset ad store them
         fitResults_wjet.append(mc_wjet.fitResult(iToy));         
         fittedPdf_wjet.append(model_Total_mc.Clone("model_Total_mc_toy_%d"%iToy)); ## clone the pdf and take the list to re-buil the pdf shape for plotting reason

         if not mc_wjet.fitResult(iToy) : continue ;
         if fitResults_wjet[len(fitResults_wjet)-1].status() != 0 : continue ;
          
         parset = mc_wjet.fitParams(iToy); ## get the parameters of the fit only those non constant
         if not parset : continue ;
         parlist = RooArgList(parset);
         
         param       = fittedPdf_wjet[len(fittedPdf_wjet)-1].getParameters(generatedData_wjet[len(generatedData_wjet)-1]); ## parameter of the new pdf
         if not param : continue ;
         parameters = RooArgList(param); 

         param_generated  = generation_model_wjet.getParameters(self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%("_WJets0","_sb_lo",self.channel,"_mlvj")));
         ### parameters used in the generation
         
         if not param_generated : continue ;
         parameters_generated = RooArgList(param_generated);

         iparNotConstant = 0;
         iGenerated = 0 ;
         iPull      = 0 ;    

         for ipar in range(parameters.getSize()): ## to clone the pdf value for plotting reasons

              if (parameters.at(ipar).GetName() == parlist.at(ipar).GetName()):
                  parameters.at(ipar).setVal(parlist.at(ipar).getVal());
                  parameters.at(ipar).setError(parlist.at(ipar).getError()) ## copy parameters in order to have the fitted shape
 
             
              if TString(parameters.at(ipar).GetName()).Contains("signal_region") and not TString(parameters.at(ipar).GetName()).Contains("number"): continue ;
             
              if iToy == 0 or len(parameterHisto_wjet)==0 : ## create parameter histo and pull 

               if not TString(parlist.at(ipar).GetName()).Contains("signal_region"):
                 if not TString(parlist.at(ipar).GetName()).Contains("number") :  
                  parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName(),"",100,-math.fabs(parlist.at(ipar).getVal()*2),math.fabs(parlist.at(ipar).getVal()*2)));
                  parameterHistoError_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_error","",100,0.,math.fabs(parlist.at(ipar).getError()*2)));
                  if options.fgen == options.fres :
                   parameterHistoPull_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_pull","",35,-3,3));
                 else : 
                  parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName(),"",50,math.fabs(parlist.at(ipar).getVal())/2,math.fabs(parlist.at(ipar).getVal())*2));
                  parameterHistoError_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_error","",100,0.,math.fabs(parlist.at(ipar).getError()*2)));
                  parameterHistoPull_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_pull","",35,-3,3)); 
               else:
                parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName(),"",100,-25,25));
                parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_fraction","",100,-100,100));                
                parameterHistoError_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_error","",100,0.,math.fabs(parlist.at(ipar).getError()*10)));
                parameterHistoPull_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_pull","",35,-3,3)); 
                 

              if not TString(parlist.at(ipar).GetName()).Contains("signal_region"): ## fill pulls and parameters histo
               parameterHisto_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal());
               parameterHistoError_wjet[iparNotConstant].Fill(parlist.at(ipar).getError());

               if not TString(parlist.at(ipar).GetName()).Contains("number") and options.fgen == options.fres:
                parameterHistoPull_wjet[iPull].Fill((parlist.at(ipar).getVal()-parameters_generated.at(iGenerated).getVal())/parlist.at(ipar).getError());
                iPull = iPull +1;
               elif TString(parlist.at(ipar).GetName()).Contains("number"):
                parameterHistoPull_wjet[iPull].Fill((parlist.at(ipar).getVal()-parlist.find("ngen").getVal())/parlist.at(ipar).getError());
                iPull = iPull +1;                
               iGenerated = iGenerated +1 ;
              else:
                ### fill with the fraction of the fitted signal events (positive or negative) divided by what is given by mc -> signal strentght in a only bkg generation  
                parameterHisto_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal());
                parameterHisto_wjet[iparNotConstant+1].Fill(parlist.at(ipar).getVal()/self.workspace4bias_.var("rrv_number_"+self.signal_sample+"_signal_region_"+self.channel+"_mlvj").getVal());
                parameterHistoError_wjet[iparNotConstant].Fill(parlist.at(ipar).getError());  
                parameterHistoPull_wjet[iPull].Fill(parlist.at(ipar).getVal()/parlist.at(ipar).getError());
                iPull = iPull +1;
                    
              iparNotConstant = iparNotConstant+1;

         ## fill chi2, NNLL
         if len(chi2distribution_wjet)==0 and parlist.find("chi2red"):
          chi2distribution_wjet.append(ROOT.TH1F("chi2distribution_wjet","",50,0.,math.fabs(parlist.find("chi2red").getVal())*4));
         if len(nLLdistribution_wjet)==0 and parlist.find("NLL"):
          nLLdistribution_wjet.append(ROOT.TH1F("nLLdistribution_wjet","",50,math.fabs(parlist.find("NLL").getVal())*0.5,math.fabs(parlist.find("NLL").getVal())*2));

         if parlist.find("chi2red") :
          chi2distribution_wjet[0].Fill(parlist.find("chi2red").getVal());
         if parlist.find("NLL") : 
          nLLdistribution_wjet[0].Fill(parlist.find("NLL").getVal());

      ### Plot in Canvas + Gaussian fit of each histogram 
      canvas_generatedToys_wjet     = []; ## canvas for show generated distribution + fits
      canvas_parameters_wjet        = []; ## canvas parameters
      canvas_parameters_err_wjet    = [];
      canvas_parameters_pull_wjet   = [];
      canvas_chi2_wjet              = [];
      canvas_chi2_wjet_frame        = [];
      canvas_nLL_wjet               = [];

      chi2distribution_wjet_frame   = [];

      ### print the canvas of the single jobs
      for iObj in range(len(generatedData_wjet)):
       storethisPlot = 0;
       if options.nexp <= 10 :
           storethisPlot = 1;
       elif options.nexp <= 100 and iObj%3 == 0:
           storethisPlot = 1;
       elif options.nexp > 100 and iObj%10 == 0:
           storethisPlot = 1;
                   
       parameters = fittedPdf_wjet[iObj].getParameters(generatedData_wjet[iObj]).selectByAttrib("Constant",kFALSE) ;

       if parameters :
        wjet_binned   = generatedData_wjet[iObj].binnedClone();
        ChiSquare = fittedPdf_wjet[iObj].createChi2(wjet_binned,RooFit.Extended(kTRUE),RooFit.SumW2Error(kTRUE));
        if len (chi2distribution_wjet_frame) ==0:
         chi2distribution_wjet_frame.append(ROOT.TH1F("chi2distribution_wjet_frame","",50,(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize())*0.5),(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize())*4)));
         
        chi2distribution_wjet_frame[0].Fill(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize()));

       if options.storeplot and  storethisPlot == 1:
         mplot = self.workspace4bias_.var("rrv_mass_lvj").frame(RooFit.Title("frame_generatedToys_%d"%iObj), RooFit.Bins(self.workspace4bias_.var("rrv_mass_lvj").getBins()));      
         generatedData_wjet[iObj].plotOn(mplot,RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
         draw_error_band_extendPdf(generatedData_wjet[iObj], fittedPdf_wjet[iObj], fitResults_wjet[iObj],mplot,2,"L");
         generatedData_wjet[iObj].plotOn(mplot,RooFit.MarkerSize(1.5),RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0),RooFit.Name(generatedData_wjet[iObj].GetName()+"_curve"));
         fittedPdf_wjet[iObj].plotOn(mplot,RooFit.Name(fittedPdf_wjet[iObj].GetName()+"_curve"));
                     
         canvas_generatedToys_wjet.append(ROOT.TCanvas("canvas_generatedToys_wjet_%d"%iObj,""));
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].cd();
         ROOT.SetOwnership(canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1], False);
         pad1 = ROOT.TPad("pad1_%d"%iObj,"pad1_%d"%iObj,0.,0.24,0.99,1. );
         ROOT.SetOwnership(pad1, False);
         pad1.Draw();
         pad2 = ROOT.TPad("pad2_%d"%iObj,"pad2_%d"%iObj,0.,0.,0.99,0.24 );
         ROOT.SetOwnership(pad2, False);
         pad2.Draw();
         pad1.cd();
         mplot.GetXaxis().SetTitleOffset(1.1);
         mplot.GetYaxis().SetTitleOffset(1.3);
         mplot.GetXaxis().SetTitleSize(0.05);
         mplot.GetYaxis().SetTitleSize(0.05);
         mplot.GetXaxis().SetLabelSize(0.045);
         mplot.GetYaxis().SetLabelSize(0.045);
         mplot.Draw();
         banner = self.banner4Plot() ;
         banner.Draw();
                                
         pad2.cd();
         mplot_pull = self.get_pull(self.workspace4bias_.var("rrv_mass_lvj"), mplot);
         mplot_pull.Draw();
         mplot_pull.GetXaxis().SetLabelSize(0.15);
         mplot_pull.GetYaxis().SetLabelSize(0.15);
         mplot_pull.GetYaxis().SetTitleSize(0.15);
         mplot_pull.GetYaxis().SetNdivisions(205);

         if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
          os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));

         pad1.SetLogy();
         pad1.Update();
         mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*100);
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].Update()

         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s_log.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s_log.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));
        

      ### print plots of the parameters
      for ipar in range(len(parameterHisto_wjet)):
       canvas_parameters_wjet.append(ROOT.TCanvas("canvas_parameters_wjet_%s"%parameterHisto_wjet[ipar].GetName(),""));
       canvas_parameters_wjet[len(canvas_parameters_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_wjet[len(canvas_parameters_wjet)-1], False);
       parameterHisto_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHisto_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHisto_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHisto_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHisto_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHisto_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHisto_wjet[ipar].GetXaxis().SetTitle(parameterHisto_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian%d"%ipar,"gaus",parameterHisto_wjet[ipar].GetXaxis().GetXmin(),parameterHisto_wjet[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHisto_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHisto_wjet[ipar].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
         
       canvas_parameters_wjet[len(canvas_parameters_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_wjet[len(canvas_parameters_wjet)-1].GetName()));
       canvas_parameters_wjet[len(canvas_parameters_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_wjet[len(canvas_parameters_wjet)-1].GetName()));

      ### print the parameters error distribution 
      for ipar in range(len(parameterHistoError_wjet)):
       canvas_parameters_err_wjet.append(ROOT.TCanvas("canvas_parameterHistoError_wjet_%s"%parameterHistoError_wjet[ipar].GetName(),""));
       canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1], False);
       parameterHistoError_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHistoError_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHistoError_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHistoError_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHistoError_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHistoError_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHistoError_wjet[ipar].GetXaxis().SetTitle(parameterHistoError_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian%d"%ipar,"gaus",parameterHistoError_wjet[ipar].GetXaxis().GetXmin(),parameterHistoError_wjet[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHistoError_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHistoError_wjet[ipar].Draw("E");
        
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
         
       canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].GetName()));
       canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].GetName()));
 
      ## print pulls of each parameter
      for ipar in range(len(parameterHistoPull_wjet)):
       canvas_parameters_pull_wjet.append(ROOT.TCanvas("canvas_parameters_pull_wjet_%s"%parameterHistoPull_wjet[ipar].GetName(),""));
       canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1], False);
       parameterHistoPull_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHistoPull_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHistoPull_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHistoPull_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHistoPull_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHistoPull_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHistoPull_wjet[ipar].GetXaxis().SetTitle(parameterHistoPull_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian%d"%ipar,"gaus",-2,2);
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHistoPull_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHistoPull_wjet[ipar].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
         
       canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].GetName()));
       canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].GetName()));
  
      ## print chi2
      if len(chi2distribution_wjet) !=0:
       canvas_chi2_wjet.append(ROOT.TCanvas("canvas_chi2_wjet_%s"%chi2distribution_wjet[0].GetName(),""));
       canvas_chi2_wjet[0].cd();
       ROOT.SetOwnership(canvas_chi2_wjet[0], False);
       chi2distribution_wjet[0].GetXaxis().SetTitleOffset(1.1);
       chi2distribution_wjet[0].GetYaxis().SetTitleOffset(1.3);
       chi2distribution_wjet[0].GetXaxis().SetTitleSize(0.04);
       chi2distribution_wjet[0].GetYaxis().SetTitleSize(0.04);
       chi2distribution_wjet[0].GetXaxis().SetLabelSize(0.035);
       chi2distribution_wjet[0].GetYaxis().SetLabelSize(0.035);
       chi2distribution_wjet[0].GetXaxis().SetTitle(chi2distribution_wjet[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_wjet[0].GetXaxis().GetXmin(),chi2distribution_wjet[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       chi2distribution_wjet[0].Fit(Gaussian,"MSQ");
       chi2distribution_wjet[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
       canvas_chi2_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet[0].GetName()));
       canvas_chi2_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet[0].GetName()));

      ## print chi2 from frame
      if len(chi2distribution_wjet_frame) !=0:
       canvas_chi2_wjet_frame.append(ROOT.TCanvas("canvas_chi2_wjet_frame_%s"%chi2distribution_wjet_frame[0].GetName(),""));
       canvas_chi2_wjet_frame[0].cd();
       ROOT.SetOwnership(canvas_chi2_wjet_frame[0], False);
       chi2distribution_wjet_frame[0].GetXaxis().SetTitleOffset(1.1);
       chi2distribution_wjet_frame[0].GetYaxis().SetTitleOffset(1.3);
       chi2distribution_wjet_frame[0].GetXaxis().SetTitleSize(0.04);
       chi2distribution_wjet_frame[0].GetYaxis().SetTitleSize(0.04);
       chi2distribution_wjet_frame[0].GetXaxis().SetLabelSize(0.035);
       chi2distribution_wjet_frame[0].GetYaxis().SetLabelSize(0.035);
       chi2distribution_wjet_frame[0].GetXaxis().SetTitle(chi2distribution_wjet_frame[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_wjet_frame[0].GetXaxis().GetXmin(),chi2distribution_wjet_frame[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       chi2distribution_wjet_frame[0].Fit(Gaussian,"MSQ");
       chi2distribution_wjet_frame[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
       canvas_chi2_wjet_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet_frame[0].GetName()));
       canvas_chi2_wjet_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet_frame[0].GetName()));

      ## print -log(l)
      if len(nLLdistribution_wjet)!=0:
       canvas_nLL_wjet.append(ROOT.TCanvas("canvas_nLL_wjet_%s"%nLLdistribution_wjet[0].GetName(),""));
       canvas_nLL_wjet[0].cd();
       ROOT.SetOwnership(canvas_nLL_wjet[0], False);
       nLLdistribution_wjet[0].GetXaxis().SetTitleOffset(1.1);
       nLLdistribution_wjet[0].GetYaxis().SetTitleOffset(1.3);
       nLLdistribution_wjet[0].GetXaxis().SetTitleSize(0.04);
       nLLdistribution_wjet[0].GetYaxis().SetTitleSize(0.04);
       nLLdistribution_wjet[0].GetXaxis().SetLabelSize(0.035);
       nLLdistribution_wjet[0].GetYaxis().SetLabelSize(0.035); 
       nLLdistribution_wjet[0].GetXaxis().SetTitle(nLLdistribution_wjet[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",nLLdistribution_wjet[0].GetXaxis().GetXmin(),nLLdistribution_wjet[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       nLLdistribution_wjet[0].Fit(Gaussian,"MSQ");
       nLLdistribution_wjet[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
       canvas_nLL_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_wjet[0].GetName()));
       canvas_nLL_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_wjet[0].GetName()));

     ############### Make the Data analysis --> make the Entended pdf for the bkg
     constrainslist_bkg_data = [];
     model_VV_backgrounds    = self.get_VV_mlvj_Model("_sb_lo");
     model_TTbar_backgrounds = self.get_TTbar_mlvj_Model("_sb_lo");
     model_STop_backgrounds  = self.get_STop_mlvj_Model("_sb_lo");

     model_bkg_data    = self.make_Model("_data_sb_lo_fit",options.fres,"_mlvj",constrainslist_bkg_data,1);
     model_bkg_data.Print();

     ## Add the other bkg component fixed to the total model --> in the extended way
     model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(modified_signal_model_mc,model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds));
                                                                                                       
     model_Total_data.Print();
     getattr(self.workspace4bias_,"import")(model_Total_data);

      ##### generate models  --> the one fixed and already fitted
     generation_model_data = self.workspace4bias_.pdf("model_%s_%s_%s_%s_%s"%("data","WJets0","sb_lo",self.channel,"mlvj"));
     generation_model_data.Print();

     self.workspace4bias_.Print();             

     ## variable and RooMC study
     numevents_data   = self.workspace4bias_.data("rdataset4bias"+"_data"+"_sb_lo"+"_"+self.channel+"_mlvj").sumEntries();
     print"########  numevents data ",numevents_data;

     data_wjet = RooMCStudy(generation_model_data,
                            RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                            RooFit.FitModel(model_Total_data),
                            RooFit.FitOptions(RooFit.Save(kTRUE),RooFit.SumW2Error(kTRUE),RooFit.Minimizer("Minuit2"),RooFit.Extended(kTRUE)),
                            RooFit.Extended(kTRUE),
                            RooFit.Silence());


     data_wjet.Print();
     
     ## create module for chi2 evaluation
     chi2_data = RooChi2MCSModule();
     data_wjet.addModule(chi2_data)

     ## generate and fit storing the generated distribution for each toy
     data_wjet.generateAndFit(options.nexp,int(numevents_data),1);

     generatedData_data         = []; ## distribution of generated toy according to bkg only hypothesis -> RooAbsData
     fittedPdf_data             = []; ## fitted pdf signal + bkg
     fitResults_data            = [];
     parameterHisto_data        = []; ## histo of the parameters of the fitted pdf
     parameterHistoError_data   = []; ## errror on the parameters
     parameterHistoPull_data    = []; ## pull wrt the generated one
     chi2distribution_data      = [];
     nLLdistribution_data       = [];
            
     for iToy in range(options.nexp): ## loop on the toy

         generatedData_data.append(data_wjet.genData(iToy)); ## take the generated dataset ad store them
         fitResults_data.append(data_wjet.fitResult(iToy));       
         fittedPdf_data.append(model_Total_data.Clone("model_Total_data_toy_%d"%iToy)); ## clone the pdf and take the list to re-buil the pdf shape for plotting reason
         if not generatedData_data or not fitResults_data : continue ;
         if fitResults_data[len(fitResults_data)-1].status() != 0 : continue ;
         
         parset  = data_wjet.fitParams(iToy); ## get the parameters of the fit
         if not parset : continue ;
         parlist = RooArgList(parset);

         param       = fittedPdf_data[len(fittedPdf_data)-1].getParameters(generatedData_data[len(generatedData_data)-1]); ## parameter of the new pdf
         if not param : continue ;
         parameters = RooArgList(param); 

         param_generated       = generation_model_data.getParameters(self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%("_Wdata","_sb_lo",self.channel,"_mlvj"))); ### parameters used in the generation
         if not param_generated : continue ;
         parameters_generated = RooArgList(param_generated);

         iparNotConstant = 0;
         iGenerated = 0 ;
         iPull = 0;
             
         for ipar in range(parameters.getSize()): ## to clone the pdf value for plotting reasons

              if (parameters.at(ipar).GetName() == parlist.at(ipar).GetName()):
                  parameters.at(ipar).setVal(parlist.at(ipar).getVal());
                  parameters.at(ipar).setError(parlist.at(ipar).getError()) ## copy parameters in order to have the fitted shape
 
             
              if TString(parameters.at(ipar).GetName()).Contains("signal_region") and not TString(parameters.at(ipar).GetName()).Contains("number"): continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_VV_") : continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_STop_") : continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_TTbar_") : continue ;

              if iToy == 0 or len(parameterHisto_data)==0 :

               if not TString(parlist.at(ipar).GetName()).Contains("signal_region"):
                 if not TString(parlist.at(ipar).GetName()).Contains("number") :  
                  parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data","",100,-math.fabs(parlist.at(ipar).getVal()*2),math.fabs(parlist.at(ipar).getVal()*2)));
                  parameterHistoError_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_error","",100,0.,math.fabs(parlist.at(ipar).getError()*2)));
                  if options.fgen == options.fres :
                   parameterHistoPull_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_pull","",50,-3,3));
                 else : 
                  parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data","",50,math.fabs(parlist.at(ipar).getVal())/2,math.fabs(parlist.at(ipar).getVal())*2));
                  parameterHistoError_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_error","",100,0.,math.fabs(parlist.at(ipar).getError()*2)));
                  parameterHistoPull_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_pull","",50,-3,3)); 
               else:
                parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data","",100,-50,50));
                parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_fraction","",100,-100,100));                
                parameterHistoError_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_error","",100,0.,math.fabs(parlist.at(ipar).getError()*10)));
                parameterHistoPull_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_pull","",50,-3,3)); 
                 
              if not TString(parlist.at(ipar).GetName()).Contains("signal_region"): ## fill pulls and parameters histo
               parameterHistoError_data[iparNotConstant].Fill(parlist.at(ipar).getError());
               
               if not TString(parlist.at(ipar).GetName()).Contains("number") :
                parameterHisto_data[iparNotConstant].Fill(parlist.at(ipar).getVal());
                if  options.fgen == options.fres:                                            
                 parameterHistoPull_data[iPull].Fill((parlist.at(ipar).getVal()-parameters_generated.at(iGenerated).getVal())/parlist.at(ipar).getError());
                 iPull = iPull +1;
               elif TString(parlist.at(ipar).GetName()).Contains("number"):
                parameterHisto_data[iparNotConstant].Fill(parlist.at(ipar).getVal()+self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)).getVal()+self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)).getVal()+self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)).getVal());
                parameterHistoPull_data[iPull].Fill((parlist.at(ipar).getVal()+self.workspace4bias_.var("rrv_number_VV_sb_lo_%s_mlvj"%(self.channel)).getVal()+self.workspace4bias_.var("rrv_number_TTbar_sb_lo_%s_mlvj"%(self.channel)).getVal()+self.workspace4bias_.var("rrv_number_STop_sb_lo_%s_mlvj"%(self.channel)).getVal()-parlist.find("ngen").getVal())/parlist.at(ipar).getError());
                iPull = iPull +1;                
               iGenerated = iGenerated +1 ;
              else:
                ### fill with the fraction of the fitted signal events (positive or negative) divided by what is given by mc -> signal strentght in a only bkg generation  
                parameterHisto_data[iparNotConstant].Fill(parlist.at(ipar).getVal());
                parameterHisto_data[iparNotConstant+1].Fill(parlist.at(ipar).getVal()/self.workspace4bias_.var("rrv_number_"+self.signal_sample+"_signal_region_"+self.channel+"_mlvj").getVal());
                parameterHistoError_data[iparNotConstant].Fill(parlist.at(ipar).getError());  

                parameterHistoPull_data[iPull].Fill(parlist.at(ipar).getVal()/parlist.at(ipar).getError());
                iPull = iPull +1;
                    
              iparNotConstant = iparNotConstant+1;

         ## fill chi2, NNLL
         if len(nLLdistribution_data)==0 and parlist.find("NLL"):
             nLLdistribution_data.append(ROOT.TH1F("nLLdistribution_data","",50,math.fabs(parlist.find("NLL").getVal())*0.5,math.fabs(parlist.find("NLL").getVal())*2));
             
         if len(chi2distribution_data)==0 and parlist.find("chi2red"):
             chi2distribution_data.append(ROOT.TH1F("chi2distribution_data","",50,0.,parlist.find("chi2red").getVal()*5));

         if parlist.find("chi2red") :
          chi2distribution_data[0].Fill(parlist.find("chi2red").getVal());
         if parlist.find("NLL"): 
          nLLdistribution_data[0].Fill(parlist.find("NLL").getVal());

     ### Plot in Canvas + Gaussian fit of each histogram 
     canvas_generatedToys_data     = []; ## canvas for show generated distribution + fits
     canvas_parameters_data        = []; ## canvas parameters
     canvas_parameters_err_data    = [];
     canvas_parameters_pull_data   = [];
     canvas_chi2_data              = [];
     canvas_nLL_data               = [];
     canvas_chi2_data_frame        = [];     

     chi2distribution_data_frame   = [];
     ### print the canvas of the single jobs         
     for iObj in range(len(generatedData_data)):
       storethisPlot = 0;

       if options.nexp <= 10 :
           storethisPlot = 1;
       elif options.nexp <= 100 and iObj%3 == 0:
           storethisPlot = 1;
       elif options.nexp > 100 and iObj%10 == 0:
          storethisPlot = 1;

       parameters = fittedPdf_data[iObj].getParameters(generatedData_data[iObj]).selectByAttrib("Constant",kFALSE) ;
       if parameters :
        data_binned   = generatedData_data[iObj].binnedClone();
        ChiSquare = fittedPdf_data[iObj].createChi2(data_binned,RooFit.Extended(kTRUE),RooFit.SumW2Error(kTRUE));
        if len (chi2distribution_data_frame) ==0:
         chi2distribution_data_frame.append(ROOT.TH1F("chi2distribution_data_frame","",50,(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize())*0.5),(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize())*4)));
         
        chi2distribution_data_frame[0].Fill(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize()));
        
       if options.storeplot and  storethisPlot == 1:

        mplot = self.workspace4bias_.var("rrv_mass_lvj").frame(RooFit.Title("frame_generatedToys_data_%d"%iObj), RooFit.Bins(self.workspace4bias_.var("rrv_mass_lvj").getBins()));      
        generatedData_data[iObj].plotOn(mplot,RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        draw_error_band_extendPdf(generatedData_data[iObj], fittedPdf_data[iObj], fitResults_data[iObj],mplot,2,"L");
        generatedData_data[iObj].plotOn(mplot,RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0),RooFit.Name(generatedData_wjet[iObj].GetName()+"_curve"));
        fittedPdf_data[iObj].plotOn(mplot, RooFit.Name(fittedPdf_data[iObj].GetName()+"_curve"));
            
        canvas_generatedToys_data.append(ROOT.TCanvas("canvas_generatedToys_data_%d"%iObj,""));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].cd();
        ROOT.SetOwnership(canvas_generatedToys_data[len(canvas_generatedToys_data)-1], False);
        pad1 = ROOT.TPad("pad1_%d"%iObj,"pad1_%d"%iObj,0.,0.24,0.99,1. );
        ROOT.SetOwnership(pad1, False);
        pad1.Draw();
        pad2 = ROOT.TPad("pad2_%d"%iObj,"pad2_%d"%iObj,0.,0.,0.99,0.24 );
        ROOT.SetOwnership(pad2, False);
        pad2.Draw();
        pad1.cd();
        mplot = self.workspace4bias_.var("rrv_mass_lvj").frame(RooFit.Title("frame_generatedToys_%d"%iObj), RooFit.Bins(self.workspace4bias_.var("rrv_mass_lvj").getBins()));      
        mplot.GetXaxis().SetTitleOffset(1.1);
        mplot.GetYaxis().SetTitleOffset(1.3);
        mplot.GetXaxis().SetTitleSize(0.05);
        mplot.GetYaxis().SetTitleSize(0.05);
        mplot.GetXaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetLabelSize(0.045);
        mplot.Draw();
        banner = self.banner4Plot() ;
        banner.Draw();
                                
        pad2.cd();
        mplot_pull = self.get_pull(self.workspace4bias_.var("rrv_mass_lvj"), mplot);
        mplot_pull.Draw();
        mplot_pull.GetXaxis().SetLabelSize(0.15);
        mplot_pull.GetYaxis().SetLabelSize(0.15);
        mplot_pull.GetYaxis().SetTitleSize(0.15);
        mplot_pull.GetYaxis().SetNdivisions(205);


        if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));

        pad1.SetLogy();
        pad1.Update();
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*100);
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].Update()

        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s_log.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s_log.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        

     ### print plots of the parameters
     for ipar in range(len(parameterHisto_data)):
      canvas_parameters_data.append(ROOT.TCanvas("canvas_parameters_data_%s"%parameterHisto_data[ipar].GetName(),""));
      canvas_parameters_data[len(canvas_parameters_data)-1].cd();
      ROOT.SetOwnership(canvas_parameters_data[len(canvas_parameters_data)-1], False);
      parameterHisto_data[ipar].GetXaxis().SetTitleOffset(1.1);
      parameterHisto_data[ipar].GetYaxis().SetTitleOffset(1.3);
      parameterHisto_data[ipar].GetXaxis().SetTitleSize(0.04);
      parameterHisto_data[ipar].GetYaxis().SetTitleSize(0.04);
      parameterHisto_data[ipar].GetXaxis().SetLabelSize(0.035);
      parameterHisto_data[ipar].GetYaxis().SetLabelSize(0.035);
      parameterHisto_data[ipar].GetXaxis().SetTitle(parameterHisto_data[ipar].GetName());
      Gaussian = ROOT.TF1("Gaussian%d"%ipar,"gaus",parameterHisto_data[ipar].GetXaxis().GetXmin(),parameterHisto_data[ipar].GetXaxis().GetXmax());
      Gaussian.SetLineColor(kBlue);
      Gaussian.SetLineWidth(2);
      parameterHisto_data[ipar].Fit(Gaussian,"MSQ");
      parameterHisto_data[ipar].Draw("E");
      
      if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

      canvas_parameters_data[len(canvas_parameters_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_data[len(canvas_parameters_data)-1].GetName()));
      canvas_parameters_data[len(canvas_parameters_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_data[len(canvas_parameters_data)-1].GetName()));

     ### print the parameters error distribution 
     for ipar in range(len(parameterHistoError_data)):
      canvas_parameters_err_data.append(ROOT.TCanvas("canvas_parameterHistoError_data_%s"%parameterHistoError_data[ipar].GetName(),""));
      canvas_parameters_err_data[len(canvas_parameters_err_data)-1].cd();
      ROOT.SetOwnership(canvas_parameters_err_data[len(canvas_parameters_err_data)-1], False);
      parameterHistoError_data[ipar].GetXaxis().SetTitleOffset(1.1);
      parameterHistoError_data[ipar].GetYaxis().SetTitleOffset(1.3);
      parameterHistoError_data[ipar].GetXaxis().SetTitleSize(0.04);
      parameterHistoError_data[ipar].GetYaxis().SetTitleSize(0.04);
      parameterHistoError_data[ipar].GetXaxis().SetLabelSize(0.035);
      parameterHistoError_data[ipar].GetYaxis().SetLabelSize(0.035);
      parameterHistoError_data[ipar].GetXaxis().SetTitle(parameterHistoError_data[ipar].GetName());
      Gaussian = ROOT.TF1("Gaussian%d"%ipar,"gaus",parameterHistoError_data[ipar].GetXaxis().GetXmin(),parameterHistoError_data[ipar].GetXaxis().GetXmax());
      Gaussian.SetLineColor(kBlue);
      Gaussian.SetLineWidth(2);
      parameterHistoError_data[ipar].Fit(Gaussian,"MSQ");
      parameterHistoError_data[ipar].Draw("E");
      
      if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

      canvas_parameters_err_data[len(canvas_parameters_err_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_data[len(canvas_parameters_err_data)-1].GetName()));
      canvas_parameters_err_data[len(canvas_parameters_err_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_data[len(canvas_parameters_err_data)-1].GetName()));

     ## print pulls of each parameter
     for ipar in range(len(parameterHistoPull_data)):
      canvas_parameters_pull_data.append(ROOT.TCanvas("canvas_parameters_pull_data_%s"%parameterHistoPull_data[ipar].GetName(),""));
      canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].cd();
      ROOT.SetOwnership(canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1], False);
      parameterHistoPull_data[ipar].GetXaxis().SetTitleOffset(1.1);
      parameterHistoPull_data[ipar].GetYaxis().SetTitleOffset(1.3);
      parameterHistoPull_data[ipar].GetXaxis().SetTitleSize(0.04);
      parameterHistoPull_data[ipar].GetYaxis().SetTitleSize(0.04);
      parameterHistoPull_data[ipar].GetXaxis().SetLabelSize(0.035);
      parameterHistoPull_data[ipar].GetYaxis().SetLabelSize(0.035);
      parameterHistoPull_data[ipar].GetXaxis().SetTitle(parameterHistoPull_data[ipar].GetName());
      Gaussian = ROOT.TF1("Gaussian%d"%ipar,"gaus",-2,2);
      Gaussian.SetLineColor(kBlue);
      Gaussian.SetLineWidth(2);
      parameterHistoPull_data[ipar].Fit(Gaussian,"MSQ");
      parameterHistoPull_data[ipar].Draw("E");
      
      if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

      canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].GetName()));
      canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].GetName()));

     ## print chi2
     if len (chi2distribution_data)!=0: 
      canvas_chi2_data.append(ROOT.TCanvas("canvas_chi2_data_%s"%chi2distribution_data[0].GetName(),""));
      canvas_chi2_data[0].cd();
      ROOT.SetOwnership(canvas_chi2_data[0], False);
      chi2distribution_data[0].GetXaxis().SetTitleOffset(1.1);
      chi2distribution_data[0].GetYaxis().SetTitleOffset(1.3);
      chi2distribution_data[0].GetXaxis().SetTitleSize(0.04);
      chi2distribution_data[0].GetYaxis().SetTitleSize(0.04);
      chi2distribution_data[0].GetXaxis().SetLabelSize(0.035);
      chi2distribution_data[0].GetYaxis().SetLabelSize(0.035);
      chi2distribution_data[0].GetXaxis().SetTitle(chi2distribution_data[0].GetName());
      Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_data[0].GetXaxis().GetXmin(),chi2distribution_data[0].GetXaxis().GetXmax());
      Gaussian.SetLineColor(kBlue);
      Gaussian.SetLineWidth(2);
      chi2distribution_data[0].Fit(Gaussian,"MSQ");
      chi2distribution_data[0].Draw("E");
      
      if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
         
      canvas_chi2_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data[0].GetName()));
      canvas_chi2_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data[0].GetName()));

     ## print -log(l)
     if len(nLLdistribution_data)!=0: 
      canvas_nLL_data.append(ROOT.TCanvas("canvas_nLL_data_%s"%nLLdistribution_data[0].GetName(),""));
      canvas_nLL_data[0].cd();
      ROOT.SetOwnership(canvas_nLL_data[0], False);
      nLLdistribution_data[0].GetXaxis().SetTitleOffset(1.1);
      nLLdistribution_data[0].GetYaxis().SetTitleOffset(1.3);
      nLLdistribution_data[0].GetXaxis().SetTitleSize(0.04);
      nLLdistribution_data[0].GetYaxis().SetTitleSize(0.04);
      nLLdistribution_data[0].GetXaxis().SetLabelSize(0.035);
      nLLdistribution_data[0].GetYaxis().SetLabelSize(0.035);
      nLLdistribution_data[0].GetXaxis().SetTitle(nLLdistribution_data[0].GetName());
      Gaussian = ROOT.TF1("GaussianChi","gaus",nLLdistribution_data[0].GetXaxis().GetXmin(),nLLdistribution_data[0].GetXaxis().GetXmax());
      Gaussian.SetLineColor(kBlue);
      Gaussian.SetLineWidth(2);
      nLLdistribution_data[0].Fit(Gaussian,"MSQ");
      nLLdistribution_data[0].Draw("E");
      
      if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
      canvas_nLL_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_data[0].GetName()));
      canvas_nLL_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_data[0].GetName()));

     ## print chi2 from frame
     if len(chi2distribution_data_frame) !=0:
      canvas_chi2_data_frame.append(ROOT.TCanvas("canvas_chi2_data_frame_%s"%chi2distribution_data_frame[0].GetName(),""));
      canvas_chi2_data_frame[0].cd();
      ROOT.SetOwnership(canvas_chi2_data_frame[0], False);
      chi2distribution_data_frame[0].GetXaxis().SetTitleOffset(1.1);
      chi2distribution_data_frame[0].GetYaxis().SetTitleOffset(1.3);
      chi2distribution_data_frame[0].GetXaxis().SetTitleSize(0.04);
      chi2distribution_data_frame[0].GetYaxis().SetTitleSize(0.04);
      chi2distribution_data_frame[0].GetXaxis().SetLabelSize(0.035);
      chi2distribution_data_frame[0].GetYaxis().SetLabelSize(0.035);
      chi2distribution_data_frame[0].GetXaxis().SetTitle(chi2distribution_data_frame[0].GetName());
      Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_data_frame[0].GetXaxis().GetXmin(),chi2distribution_data_frame[0].GetXaxis().GetXmax());
      Gaussian.SetLineColor(kBlue);
      Gaussian.SetLineWidth(2);
      chi2distribution_data_frame[0].Fit(Gaussian,"MSQ");
      chi2distribution_data_frame[0].Draw("E");
      
      if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
      canvas_chi2_data_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data_frame[0].GetName()));
      canvas_chi2_data_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data_frame[0].GetName()));



#### Main code     
if __name__ == "__main__":

  print "###################### begin the analysis: channel %s, signla name %s, mlvj min %s, mlvj max %s, mj min %s, mj max %s, genfunction %s, fitfunction %s"%(options.channel,sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],options.fgen,options.fres);
  
  fitBiasAnalysis = doBiasStudy_mlvj (options.channel,sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),options.fgen,options.fres)
  fitBiasAnalysis.biasAnalysis();
  


                                  
           
  

