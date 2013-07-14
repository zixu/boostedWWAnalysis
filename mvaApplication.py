#! /usr/bin/env python
import os
import glob
import math
import array

from ROOT import gROOT, gStyle, gSystem, TLatex, RooArgSet, RooFit, RooRealVar
import subprocess
from subprocess import Popen

from sampleWrapperClass import *
from trainingClass      import *
from BoostedWSamples    import * 
from mvaApplication     import *

## Application class

class tmvaApplicator: 

    def __init__(self, listOfVars, weightfile):
        
        self.ListOfTrainingVariables = listOfVars;
        self.nameOfWeightFile = weightfile;
            
        # make a list of arrays
        self.reader = ROOT.TMVA.Reader("!Color:!Silent")
        self.listOfVarArray = [];
        for i in range(len(self.ListOfTrainingVariables)):
            #            curVar = array('f',[0.]);
            self.listOfVarArray.append( array('f',[0.]) );
            varString = self.ListOfTrainingVariables[i]+" := "+self.ListOfTrainingVariables[i];
            self.reader.AddVariable( varString, self.listOfVarArray[i] );
        
        #spectators
        spec1 = array('f',[0.]);
        spec2 = array('f',[0.]);
    #    spec3 = array('f',[0.]);
        self.reader.AddSpectator( "jet_pt_pr", spec1 )
        self.reader.AddSpectator( "jet_mass_pr", spec2 )
    #    reader.AddSpectator( "jet_massdrop_pr", spec3 )        

        self.reader.BookMVA("BDT",self.nameOfWeightFile);

    def eval(self, listOfVarVals):
        
        for i in range(len(listOfVarVals)):
            self.listOfVarArray[i][0] = listOfVarVals[i];
            
        return self.reader.EvaluateMVA("BDT");

class multi_tmvaApplicator:
    def __init__(self,listOfVarVals,weightfile1,weightfile2):
        self.ListOfTrainingVariables = listOfVarVals;
        self.nameOfWeightFile1 =weightfile1;
        self.nameOfWeightFile2 =weightfile2;

        # make a list of arrays
        self.reader1 = ROOT.TMVA.Reader("!Color:!Silent")
        self.reader2 = ROOT.TMVA.Reader("!Color:!Silent")
        self.listOfVarArray = [];
        for i in range(len(self.ListOfTrainingVariables)):
            self.listOfVarArray.append( array('f',[0.]) );
            varString = self.ListOfTrainingVariables[i]+" := "+self.ListOfTrainingVariables[i];
            self.reader1.AddVariable( varString, self.listOfVarArray[i] );
            self.reader2.AddVariable( varString, self.listOfVarArray[i] );
        
        #spectators
        spec1 = array('f',[0.]);
        spec2 = array('f',[0.]);
    #    spec3 = array('f',[0.]);
        self.reader1.AddSpectator( "jet_pt_pr", spec1 )
        self.reader1.AddSpectator( "jet_mass_pr", spec2 )
        self.reader2.AddSpectator( "jet_pt_pr", spec1 )
        self.reader2.AddSpectator( "jet_mass_pr", spec2 )
    #    reader.AddSpectator( "jet_massdrop_pr", spec3 )        

        self.reader1.BookMVA("BDT",self.nameOfWeightFile1);
        self.reader2.BookMVA("BDT",self.nameOfWeightFile2);

    def eval(self, listOfVarVals, Wpt):
        
        for i in range(len(listOfVarVals)):
            self.listOfVarArray[i][0] = listOfVarVals[i];
            
        #print "Wpt=%s"%(Wpt)
        #print "BDT: %s %s"%(self.reader1.EvaluateMVA("BDT"),self.reader2.EvaluateMVA("BDT"))
        if Wpt<275:
            return self.reader1.EvaluateMVA("BDT");
        else: 
            return self.reader2.EvaluateMVA("BDT");
