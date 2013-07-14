#!/usr/bin/env python

########################################
##     
##       Author: Wei Zou
##       
##       Email: weizou.pku@gmail.com
#######################################

from ROOT import *
import os
import sys
import math
import string
from array import array
from BoostedWSamples import *
from TreePlot import *

class ProducePlot:
      
     def __init__(self,CHANNEL):#CHANNEL = 'mu', or 'el'
          
          self.thstack = THStack("MC","MC")

          self.samples = Samples(CHANNEL)
          self.drawhistograms = {}
          
          self.tmpplot = TreePlot()

          self.bin_ = 0
          self.min_ = 0.
          self.max_ = 0.

          self.cutname_ = ""
          self.weightname_ = ""

          self.plotdir_ = ""
          self.setlogy_ = False

          self.channel = CHANNEL

          self.multiplicitylabel = {} 

          self.c1 = TCanvas()

     def SetBinMinMax(self,bin,min,max):
          
          self.bin_ = int(bin)
          self.min_ = float(min)
          self.max_ = float(max)
      
     def SetCutWeightName(self,cutname,weightname):
          
          self.cutname_  = cutname
          self.weightname_ = weightname

     def SetFilePath(self,filepath):
          
          self.samples.SetFilePath(filepath)
      
     def SetLumiTree(self,lumi,treename):
          
          self.samples.SetLumi(lumi)
          self.samples.SetTreeName(treename)

     def SetPlotDir(self,plotdir):
          
          self.plotdir_ = plotdir
     
     def SetLogy(self,setlogy):
         
         if(setlogy == "False"):
           self.setlogy_ = False
         
         if(setlogy == "True"):
           self.setlogy_ = True

     def SetChannel(self,Channel):

         self.channel = Channel

     def DrawTHStack(self,drawname,index,xtitle,ytitle,scalefactortxtfilename,n,ymin,latexcuttitle):
         
         self.samples.SetFileNames()
         
         print "We are ploting the variable: " + drawname

         print "We are using the cut: " + self.cutname_

         print "We are using the weight: " + self.weightname_
         
         ttbarregionsignal = ""
            
         for tmpfiletype in self.samples.GetFileNames().keys():
             #print tmpfiletype
             tmpcutname = self.cutname_
             if self.plotdir_ == "controlPlots_ttbar" and (str(tmpfiletype).find("H") != -1): tmpcutname = "&&".join(self.cutname_.split("&&")[:-1])
             self.tmpplot.SetBinMinMax(self.bin_,self.min_,self.max_)
             self.tmpplot.SetCutWeightName(tmpcutname,self.weightname_)
             tmphistogram = TH1D(tmpfiletype,tmpfiletype,self.bin_,self.min_,self.max_)
             if (tmpfiletype != "data"):
                 self.tmpplot.DrawTrees(self.samples.GetFileNames()[tmpfiletype],self.samples.GetTreeName(),drawname,index,xtitle,ytitle,"%s:%s"%(scalefactortxtfilename,tmpfiletype),self.samples.GetLumi(),tmpfiletype,self.drawhistograms,self.multiplicitylabel,ttbarregionsignal)
             else:
                 self.tmpplot.DrawTrees(self.samples.GetFileNames()[tmpfiletype],self.samples.GetTreeName(),drawname,index,xtitle,ytitle,1.0,1.0,tmpfiletype,self.drawhistograms,self.multiplicitylabel,ttbarregionsignal)
             self.drawhistograms[tmpfiletype].SetStats(kFALSE)
         
         ################################Data############################
         self.drawhistograms["data"].SetMarkerStyle(20)
         self.drawhistograms["data"].SetMarkerSize(1)
         self.drawhistograms["data"].Sumw2()
         ################################Data############################

         ##########TTbar###########################################
         tt1D = TH1D("tt1D","tt1D",self.bin_,self.min_,self.max_)
         self.drawhistograms["TTbar"].SetFillColor(kRed+1)
         tt1D = self.drawhistograms["TTbar"]
         tt1D.SetFillColor(kRed+1)
         ##########################################################
         
         #############WJet########################################
         wjet1D = TH1D("wjet1D","wjet1D",self.bin_,self.min_,self.max_)
         wjet1D.Sumw2()
         wjet1D.SetStats(kFALSE)

         if(self.samples.GetFileNames().has_key("Wbb") and self.samples.GetFileNames().has_key("Wcc") and self.samples.GetFileNames().has_key("Wlight")):
            wjet1D = self.drawhistograms["Wbb"] + self.drawhistograms["Wcc"] + self.drawhistograms["Wlight"]
         else:
            wjet1D = self.drawhistograms["WJets_Pythia"]
         wjet1D.SetFillColor(kGreen-3)
         #############WJet########################################
         
         ###########QCD############################################
         qcd1D = TH1D("qcd1D","qcd1D",self.bin_,self.min_,self.max_)
         qcd1D.Sumw2()
         qcd1D.SetStats(kFALSE)
         if (self.samples.GetFileNames().has_key("QCDBCtoE3080") and self.samples.GetFileNames().has_key("QCDBCtoE80170") and self.samples.GetFileNames().has_key("QCDEmEn3080") and self.samples.GetFileNames().has_key("QCDEmEn80170")):
            qcd1D = self.drawhistograms["QCDBCtoE3080"] + self.drawhistograms["QCDBCtoE80170"] + self.drawhistograms["QCDEmEn3080"]  + self.drawhistograms["QCDEmEn80170"]
         qcd1D.SetFillColor(kYellow)
         ###########QCD############################################
         
         ############################DiBoson#######################################
         diboson1D = TH1D("diboson1D","diboson1D",self.bin_,self.min_,self.max_)
         diboson1D.Sumw2()
         if (self.samples.GetFileNames().has_key("WW") and self.samples.GetFileNames().has_key("WZ") and self.samples.GetFileNames().has_key("ZZ")):
            diboson1D = self.drawhistograms["WW"]+ self.drawhistograms["WZ"] + self.drawhistograms["ZZ"] 
         else:
            if(self.samples.GetFileNames().has_key("WW") and self.samples.GetFileNames().has_key("WZ")):
               diboson1D = self.drawhistograms["WW"]+ self.drawhistograms["WZ"]
         diboson1D.SetStats(kFALSE)
         diboson1D.SetFillColor(kBlue)
         ############################DiBoson#######################################
         
         #####################################SumSingleTop###########################
         sumsttop1D = TH1D("sumsttop1D","sumsttop1D",self.bin_,self.min_,self.max_)
         sumsttop1D.Sumw2()
         if (self.samples.GetFileNames().has_key("tch") and self.samples.GetFileNames().has_key("tWch") and self.samples.GetFileNames().has_key("sch") and self.samples.GetFileNames().has_key("tch_bar") and self.samples.GetFileNames().has_key("tWch_bar") and self.samples.GetFileNames().has_key("sch_bar")):
            sumsttop1D = self.drawhistograms["tch"] + self.drawhistograms["tWch"] + self.drawhistograms["sch"] + self.drawhistograms["tch_bar"] + self.drawhistograms["tWch_bar"] + self.drawhistograms["sch_bar"]
         sumsttop1D.SetStats(kFALSE)
         sumsttop1D.SetFillColor(kMagenta)
         #####################################SumSingleTop###########################
         
         #####################################ZJet###################################
         zjet1D = TH1D("zjet1D","zjet1D",self.bin_,self.min_,self.max_)
         zjet1D.Sumw2()
         if (self.samples.GetFileNames().has_key("ZJets")):
            zjet1D = self.drawhistograms["ZJets"]
         zjet1D.SetStats(kFALSE)
         zjet1D.SetFillColor(kAzure-3)
         #####################################ZJet###################################

         mcbackground = [self.drawhistograms["TTbar"], wjet1D, qcd1D, diboson1D, sumsttop1D, zjet1D]

         #####################################Data-MC/MC############################
         allmc1D = TH1D("allmc1D","allmc1D",self.bin_,self.min_,self.max_)
         allmc1D.Sumw2()
         allmc1D = tt1D +  wjet1D + diboson1D + sumsttop1D +  zjet1D
         if(qcd1D.GetEntries() > 0):
            allmc1D = tt1D +  wjet1D + diboson1D + sumsttop1D +  zjet1D + qcd1D
         allmc1D.SetStats(kFALSE)

         dataminusmc1D = TH1D("dataminusmc1D","dataminusmc1D",self.bin_,self.min_,self.max_)
         dataminusmc1D.Sumw2()
         dataminusmc1D.Add(self.drawhistograms["data"],allmc1D,1.,-1.)
         dataminusmc1D.SetStats(kFALSE)

         dataovermc1D = TH1D("dataovermc1D","dataovermc1D",self.bin_,self.min_,self.max_)
         dataovermc1D.Sumw2()
         dataovermc1D.Divide(dataminusmc1D,allmc1D,1,1,"B")
         dataminusmc1D.SetStats(kFALSE)
         dataovermc1D.SetXTitle(xtitle)
         dataovermc1D.SetYTitle("(Data - MC) / MC")
         dataovermc1D.SetMarkerStyle(20)
         dataovermc1D.SetMarkerSize(1)
         dataovermc1D.SetTitle("")
         dataovermc1D.SetTitleSize(0.05,"X")
         dataovermc1D.SetTitleSize(0.05,"Y")
         dataovermc1D.GetYaxis().SetRangeUser(-1,1)
         #####################################Data-MC/MC############################
         
         datamcksresult_ = self.drawhistograms["data"].KolmogorovTest(allmc1D)

         #####################################Error_Band###########################
         xvalue = array('d')
         yvalue = array('d')
         xlefterror = array('d')
         xrighterror = array('d')
         ylowerror = array('d')
         yhigherror = array('d')

         xlist = []
         ylist = []
         xleftlist = []
         xrightlist = []
         ylowlist = []
         yhighlist = []

         for ibin in range(1,allmc1D.GetNbinsX() + 1):
             xlist.append(allmc1D.GetBinCenter(ibin))
             ylist.append(allmc1D.GetBinContent(ibin))
             xleftlist.append(0.5 * allmc1D.GetBinWidth(ibin))
             xrightlist.append(0.5 * allmc1D.GetBinWidth(ibin))
             ttbarerror = self.drawhistograms["TTbar"].GetBinContent(ibin) * 0.15
             wjeterror = wjet1D.GetBinContent(ibin) * 0.3
             lumierror = allmc1D.GetBinContent(ibin) * 0.022
             statisticerror = allmc1D.GetBinError(ibin)
             allerror = sqrt(pow(ttbarerror,2) + pow(wjeterror,2) + pow(lumierror,2) + pow(statisticerror,2))
             allerror = statisticerror             
             ylowlist.append(allerror)
             yhighlist.append(allerror)

         xvalue.fromlist(xlist)
         yvalue.fromlist(ylist)
         xlefterror.fromlist(xleftlist)
         xrighterror.fromlist(xrightlist)
         ylowerror.fromlist(ylowlist)
         yhigherror.fromlist(yhighlist)

         mc1Derror = TGraphAsymmErrors(allmc1D.GetNbinsX(),xvalue,yvalue,xlefterror,xrighterror,ylowerror,yhigherror)
         mc1Derror.SetName("MC Uncerntainty")
         mc1Derror.SetFillColor(920+3)
         mc1Derror.SetFillStyle(3008)
         #####################################Error_Band###########################
         
         self.thstack.SetHistogram(self.drawhistograms["TTbar"])
         if(self.drawhistograms["TTbar"].GetEntries() > 0):
           self.thstack.Add(self.drawhistograms["TTbar"])
         if(wjet1D.GetEntries() > 0):
           self.thstack.Add(wjet1D)
         if(sumsttop1D.GetEntries() > 0):
            self.thstack.Add(sumsttop1D)
         if(self.samples.GetFileNames().has_key("ZJets")):
           if(zjet1D.GetEntries() > 0):
              self.thstack.Add(zjet1D)
         if(diboson1D.GetEntries() > 0):
           self.thstack.Add(diboson1D)
         if(qcd1D.GetEntries() > 0):
           self.thstack.Add(qcd1D)
         self.thstack.SetMaximum(float(n) * self.thstack.GetMaximum())
         self.thstack.SetMinimum(float(ymin))
                 
         #L = TLegend(0.66,0.65,0.93,0.93)
         L = TLegend(0.7270115,0.6510417,0.9971264,0.9300595)
         L.SetBorderSize(0)
         L.SetLineStyle(0)
         L.SetTextFont(42)
         L.SetFillStyle(0)
         L.SetMargin(0.12)
         L.SetTextSize(0.025)
         L.SetFillColor(10)
         L.SetBorderSize(0)
         if(self.drawhistograms["data"].GetEntries() > 0):
           L.AddEntry(self.drawhistograms["data"],"Data", "lp")
         if(self.drawhistograms["TTbar"].GetEntries() > 0):
           L.AddEntry(self.drawhistograms["TTbar"],"t#bar{t}", "f")
         if(wjet1D.GetEntries() > 0):
           L.AddEntry(wjet1D,"W#rightarrowl#nu", "f")
         if(self.samples.GetFileNames().has_key("ZJets")):
           if(zjet1D.GetEntries() > 0):
             L.AddEntry(zjet1D,"Z/#gamma^{\*}#rightarrowl^{+}l^{-}", "f")
         if(qcd1D.GetEntries() > 0):
           L.AddEntry(qcd1D,"QCD", "f")
         if(diboson1D.GetEntries() > 0):
           L.AddEntry(diboson1D,"Dibosons", "f")
         if(sumsttop1D.GetEntries() > 0):
           L.AddEntry(sumsttop1D,"Single-Top","f")
         if(self.multiplicitylabel["ggH600"] != 1.0):
           L.AddEntry(self.drawhistograms["ggH600"],"ggHWW 600 #times %d"%self.multiplicitylabel["ggH600"],"f")
         else:
           L.AddEntry(self.drawhistograms["ggH600"],"ggHWW 600","f")

         banner = TLatex(0.25,0.88,("#splitline{CMS Preliminary}{%.1f fb^{-1} at #sqrt{s}=8TeV %s+jets}"%(self.samples.GetLumi(),self.channel)))
         banner.SetNDC()
         banner.SetTextSize(0.035)

         tl = TLatex(0.25,0.85,("#splitline{KS=%.4f}{%s}"%(datamcksresult_,latexcuttitle)))
         tl.SetNDC()
         tl.SetTextSize(0.06)

         self.c1 = TCanvas(drawname + self.cutname_,drawname + self.cutname_,10,10,700,700)
         self.c1.cd()
         pad1 = TPad("pad1","pad1",0.00,0.25,1.00,0.97)
         pad2 = TPad("pad2","pad2",0.00,0.00,1.00,0.25)
         pad1.SetFillColor(0)
         pad2.SetFillColor(0)
         pad1.Draw()
         pad2.Draw()
         pad1.SetTicks(1,1)
         pad2.SetTicks(1,1)
         pad2.SetGridx()
         pad2.SetGridy()

         pad1.cd()
         if(self.setlogy_):
            pad1.SetLogy()
         self.thstack.Draw("hist")
         mc1Derror.Draw("2 same")
         self.drawhistograms["data"].Draw("Esame")
         self.drawhistograms["ggH600"].SetLineStyle(2)
         self.drawhistograms["ggH600"].Draw("histsame")
         L.Draw()
         banner.Draw()
         
         pad2.cd()
         dataovermc1D.Draw()
         tl.Draw()
         
         PWD = os.getcwd() + "/"

         if os.path.isdir(PWD + self.plotdir_):
            print "%s direcotry has aleady been created"%(self.plotdir_)
         else:
            os.system("mkdir %s"%(PWD + self.plotdir_))
         
         tmplatexcuttitle = latexcuttitle
         tmplatexcuttitle = "_".join(tmplatexcuttitle.split(" "))
         tmplatexcuttitle = "".join(tmplatexcuttitle.split("&&"))
         tmplatexcuttitle = "_".join(tmplatexcuttitle.split("__"))
         tmplatexcuttitle = string.replace(tmplatexcuttitle,">","large")
         tmplatexcuttitle = string.replace(tmplatexcuttitle,"=","euqal")
         tmplatexcuttitle = string.replace(tmplatexcuttitle,"<","small")
         tmplatexcuttitle = string.replace(tmplatexcuttitle,"#","")
         
         self.c1.Print( PWD  + self.plotdir_ + "/" + self.channel + "_" + drawname + "_" +  tmplatexcuttitle + ".png")
         self.c1.Print( PWD  + self.plotdir_ + "/" + self.channel + "_" + drawname + "_" +  tmplatexcuttitle + ".eps")
         self.c1.Print( PWD  + self.plotdir_ + "/" + self.channel + "_" + drawname + "_" +  tmplatexcuttitle + ".pdf")
         self.c1.Update()

         #raw_input( 'Press ENTER to continue\n ' )

