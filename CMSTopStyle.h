// -*- C++ -*-

#include "TString.h"

// To use in CINT:
// root [0] #include "CMSTopStyle_v2.cc"
// root [1] CMSTopStyle style           
// root [2] style.TtbarColor
// (int)633

// To use in PyRoot:
// >>> import ROOT
// >>> ROOT.gROOT.ProcessLine ('.L CMSTopStyle.cc+')
// >>> style = ROOT.CMSTopStyle()
// >>> style.TtbarColor
// 633

class CMSTopStyle
{
   public:

      CMSTopStyle();

      // setup everything a la ICHEP
      void setupICHEPv1();

      void setupDefault();

      ////////////
      // Colors //
      ////////////

      // see http://root.cern.ch/root/html/MACRO_TColor_3_wheel.gif
      // for color wheel

      // ttbar
      int TtbarColor;        
      int TtbarOtherColor;   

      // single top                 
      int SingleTopColor;    
      int ST_tWColor;        
      int ST_t_sColor;       

      // W + jets                   
      int WJetsColor;        
      int WLFJetsColor;      
      int WbbJetsColor;      
      int WccJetsColor;      
      int WcJetsColor;       

      // Z + jets / DY              
      int DYZJetsColor;      
      int DYZTauTauJetsColor;
      // we should add Z + flavor too

      // QCD
      int QCDColor;          

      // Other EW
      int DibosonsColor;     
      int GammaJetsColor;    
      int WGammaColor;       


      ///////////
      // Fills //
      ///////////

      // ttbar
      int TtbarFill;        
      int TtbarOtherFill;    

      // single top                 
      int SingleTopFill;    
      int ST_tWFill;        
      int ST_t_sFill;        

      // W + jets                   
      int WJetsFill;        
      int WLFJetsFill;      
      int WbbJetsFill;      
      int WccJetsFill;      
      int WcJetsFill;       

      // Z + jets / DY              
      int DYZJetsFill;      
      int DYZTauTauJetsFill;
      // we should add Z + flavor too

      // QCD
      int QCDFill;          

      // Other EW
      int DibosonsFill;      
      int GammaJetsFill;    
      int WGammaFill;        

      //////////
      // Text //
      //////////

      // ttbar
      TString TtbarText;        
      TString TtbarOtherText;   

      // single top                 
      TString SingleTopText;    
      TString ST_tWText;        
      TString ST_t_sText;       

      // W + jets                   
      TString WJetsText;        
      TString WLFJetsText;      
      TString WbbJetsText;      
      TString WccJetsText;      
      TString WcJetsText;       

      // Z + jets / DY              
      TString DYZJetsText;      
      TString DYZemuJetsText;   
      TString DYZTauTauJetsText;
      // we should add Z + flavor too

      // QCD
      TString QCDText;          

      // Other EW
      TString DibosonsText;     
      TString GammaJetsText;    
      TString WGammaText;       

};

