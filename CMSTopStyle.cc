// -*- C++ -*-

#include "CMSTopStyle.h"

// To use in CINT:
// root [0] #include "CMSTopStyle.h"
// root [1] CMSTopStyle::kTtbarColor
// (const int)633
// root [2] using namespace CMSTopStyle
// root [3] kTtbarColor
// (const int)633

// To use in PyRoot:
// >>> import ROOT
// >>> ROOT.gROOT.ProcessLine('#include "CMSTopStyle.h"')
// >>> ROOT.CMSTopStyle.kTtbarColor
// 633
// >>> style = ROOT.CMSTopStyle
// >>> style.kTtbarColor
// 633

////////////
// Colors //
////////////

// see http://root.cern.ch/root/html/MACRO_TColor_3_wheel.gif
// for color wheel

CMSTopStyle::CMSTopStyle()
{
   setupDefault();
}

void
CMSTopStyle::setupDefault()
{
   // ttbar
   TtbarColor         = kRed + 1;
   TtbarOtherColor    = kRed - 7;    // Try to avoid
                               
   // single top                 
   SingleTopColor     = kOrange;
   ST_tWColor         = kOrange;
   ST_t_sColor        = kOrange + 2; // Try to avoid
                               
   // W + jets                   
   WJetsColor         = kGreen - 9;
   WLFJetsColor       = kGreen - 9;
   WbbJetsColor       = kGreen - 3;
   WccJetsColor       = kGreen - 5;
   WcJetsColor        = kGreen - 7;
                               
   // Z + jets / DY              
   DYZJetsColor       = kAzure - 2;
   DYZTauTauJetsColor = kAzure + 8;
   // we should add Z + flavor too
 
   // QCD
   QCDColor           = kYellow;
 
   // Other EW
   DibosonsColor      = kMagenta;  // White is a bad color...
   GammaJetsColor     = kGray + 2;
   WGammaColor        = kGray + 3;
 
 
   ///////////
   // Fills //
   ///////////
 
   // ttbar
   TtbarFill         = 1001;
   TtbarOtherFill    = 1001; 
                               
   // single top                 
   SingleTopFill     = 1001;
   ST_tWFill         = 1001;
   ST_t_sFill        = 1001; 
                               
   // W + jets                   
   WJetsFill         = 1001;
   WLFJetsFill       = 1001;
   WbbJetsFill       = 1001;
   WccJetsFill       = 1001;
   WcJetsFill        = 1001;
                               
   // Z + jets / DY              
   DYZJetsFill       = 1001;
   DYZTauTauJetsFill = 1001;
   // we should add Z + flavor too
 
   // QCD
   QCDFill           = 1001;
 
   // Other EW
   DibosonsFill      = 1001; 
   GammaJetsFill     = 1001;
   WGammaFill        = 1001;
 
   //////////
   // Text //
   //////////
 
   // ttbar
   TtbarText         = "t#bar{t}";
   TtbarOtherText    = "t#bar{t} other"; 
                               
   // single top                 
   SingleTopText     = "Single-Top";
   ST_tWText         = "tW";
   ST_t_sText        = "s+t channels"; 
                               
   // W + jets                   
   WJetsText         = "W#rightarrowl#nu";
   WLFJetsText       = "W#rightarrowl#nu+light jets";
   WbbJetsText       = "W#rightarrowl#nu+b#bar{b}";
   WccJetsText       = "W#rightarrowl#nu+c#bar{c}";
   WcJetsText        = "W#rightarrowl#nu+c";
                               
   // Z + jets / DY              
   DYZJetsText       = "Z/#gamma*#rightarrowl^{+}l^{-}";
   DYZemuJetsText    = "Z/#gamma*#rightarrowl^{+}l^{-} (l=e,#mu)";
   DYZTauTauJetsText = "Z/#gamma*#rightarrow#tau^{+}#tau^{-}"; 
   // we should add Z + flavor too
 
   // QCD
   QCDText           = "QCD";
 
   // Other EW
   DibosonsText      = "Dibosons"; 
   GammaJetsText     = "#gamma";
   WGammaText        = "W#gamma";
 
}

void 
CMSTopStyle::setupICHEPv1()
{
   // single top                 
   SingleTopColor     = kMagenta;
   ST_tWColor         = kMagenta;
   ST_t_sColor        = kMagenta + 2; // Try to avoid
                               
   // W + jets                   
   WJetsColor         = kGreen - 3;
   WLFJetsColor       = kGreen - 3;
   WbbJetsColor       = kGreen - 9;
   WccJetsColor       = kGreen - 7;
   WcJetsColor        = kGreen - 5;
                               
   // Other EW
   DibosonsColor      = kWhite;
   GammaJetsColor     = kOrange - 3;
   WGammaColor        = kGray;
}
