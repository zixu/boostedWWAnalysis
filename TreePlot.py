#!/usr/bin/env python

########################################
##     
##       Author: Wei Zou
##       
##       Email: weizou.pku@gmail.com
#######################################

from ROOT import *
import ROOT
import os.path

class TreePlot:
      
      def __init__(self):

          self.histogram = {} 
          
          self.bin_ = 0
          self.min_ = 0.
          self.max_ = 0.
          
          self.cutname_ = ""
          self.weightname_ = ""
          
          self.n = 0

      def SetBinMinMax(self,bin,min,max):
          
          self.bin_ = bin
          self.min_ = min
          self.max_ = max
      
      def SetCutWeightName(self,cutname,weightname):
          
          self.cutname_ = cutname
          self.weightname_ = weightname

      def DrawTrees(self,filename,treename,drawname,index,xtitle,ytitle,scalename,lumi,keyname,dict,multiplicitylabel,plotdir):
          
          tfilename = TFile(filename,"READ")
          tfilename.cd()
          obj = tfilename.Get(treename)
          self.n = float(str(hash(gRandom.Gaus(1000,200)))[0:8])
          basefilename = os.path.basename(filename).rstrip(".root")
          if (len(str(index).split(":")) > 1):
             index0 = float(index.split(":")[0])
             index1 = float(index.split(":")[1])
             obj.Draw("%s[%d][%d]>>h%d(%d,%f,%f)"%(drawname,index0,index1,self.n,self.bin_,self.min_,self.max_),"((%s) * (%s))"%(self.weightname_,self.cutname_),"goff")
          else: 
             index = float(index)
             obj.Draw("%s[%d]>>h%d(%d,%f,%f)"%(drawname,index,self.n,self.bin_,self.min_,self.max_),"((%s) * (%s))"%(self.weightname_,self.cutname_),"goff")
          if(hasattr(ROOT,"h%d"%self.n)):
             getattr(ROOT,"h%d"%self.n).SetDirectory(0)
             self.histogram[basefilename] = getattr(ROOT,"h%d"%self.n).Clone(basefilename)
             self.histogram[basefilename].SetDirectory(0)
             self.histogram[basefilename].SetName(basefilename)
             self.histogram[basefilename].Sumw2()
             self.histogram[basefilename].SetTitle(basefilename)
             self.histogram[basefilename].GetXaxis().SetTitle(xtitle)
             self.histogram[basefilename].GetYaxis().SetTitle("Events/%.2f %s"%(self.histogram[basefilename].GetBinWidth(1),ytitle))
             scalefactor = 1.0
             if len(str(scalename).split(":")) > 1:
                keySF = scalename.split(":")[1]
                SFfile = open(scalename.split(":")[0])
                for sfline in SFfile:
                   if sfline.find("#")!=-1: continue
                   if sfline.find(keySF) != -1:
                      scalefactor = float(sfline.split()[1])
                      multiplicitylabel[keyname] = 1.0
                      if len(sfline.split()) > 2:
                         multiplicitylabel[keyname] = float(sfline.split()[2])
                      scalefactor = scalefactor * multiplicitylabel[keyname]
                      break
                SFfile.close()
             else: scalefactor = scalename
             self.histogram[basefilename].Scale(scalefactor * float(lumi))
          else:
             self.histogram[basefilename] = TH1D(basefilename,basefilename,self.bin_,self.min_,self.max_) 
             self.histogram[basefilename].SetDirectory(0)
          
          dict[keyname] = self.histogram[basefilename]
          dict[keyname].SetDirectory(0)
          tfilename.Close()
