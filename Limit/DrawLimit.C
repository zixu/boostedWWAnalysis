// tdrGrid: Turns the grid lines on (true) or off (false)

// void tdrGrid(bool gridOn) {
//   tdrStyle->SetPadGridX(gridOn);
//   tdrStyle->SetPadGridY(gridOn);
// }

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.18);
  tdrStyle->SetPadRightMargin(0.06);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.030, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.1);
  tdrStyle->SetTitleYOffset(1.3);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.03, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}



Bool_t same(Float_t a1, Float_t a2){
    if(TMath::Abs(a1-a2) < 0.001)return 1;
    else return 0;
}


void DrawLimit(char*channel, char* model="unbin", double el_lumi=19.3, double mu_lumi=19.3, bool showobs=1)
{
	setTDRStyle();
	TFile inFile(Form("higgisCombin_%s_%s.root",channel,model)) ;
    TTree* tree_limit=inFile.Get("limit");

    Double_t tmp_mh;
    Float_t  tmp_quantileExpected;
    Double_t tmp_limit;

    limit->SetBranchAddress("limit",&tmp_limit);
    limit->SetBranchAddress("mh",&tmp_mh);
    limit->SetBranchAddress("quantileExpected",&tmp_quantileExpected);

	const double nmass = tree_limit.GetEntries()/6.;
	double Mass[nmass] ;
	double Zero[nmass];
    for(Int_t i=0;i<nmass;i++)Zero[i]=0.;

	double Observed[nmass];
	double Median[nmass];
	double N1S [nmass];
	double P1S [nmass];
	double N2S [nmass];
	double P2S [nmass];


	int masspoint=-1;
    for(Int_t i=0;i<nmass*6;i++){
        masspoint=i/6;
        tree_limit->GetEntry(i);

        if( same(tmp_quantileExpected,0.025)) N2S[masspoint]=tmp_limit;
        if( same(tmp_quantileExpected,0.16 )) N1S[masspoint]=tmp_limit;
        if( same(tmp_quantileExpected,0.5  )) Median[masspoint]=tmp_limit;
        if( same(tmp_quantileExpected,0.84 )) P1S[masspoint]=tmp_limit;
        if( same(tmp_quantileExpected,0.975)) P2S[masspoint]=tmp_limit;
        if( same(tmp_quantileExpected,-1.0 )) Observed[masspoint]=tmp_limit;
        Mass[masspoint]=tmp_mh;
    }

	for (int i=0 ;i<nmass ;i++ ) {
		N1S [i] = Median[i] - N1S [i];
		P1S [i] = P1S [i]   - Median[i];
		N2S [i] = Median[i] - N2S [i];
		P2S [i] = P2S [i]   - Median[i];
	}

	TGraph *likelihd_limit_d = new TGraph(nmass,Mass,Observed);                                                                                              
	likelihd_limit_d->SetLineColor(kBlack);                                                                                                                  
	likelihd_limit_d->SetLineWidth(2);                                                                                                                       
	likelihd_limit_d->SetLineStyle(1);
    likelihd_limit_d->SetMarkerStyle(20);

	TGraph *likelihd_limit_c = new TGraph(nmass,Mass,Median);                                                                                                
	likelihd_limit_c->SetLineColor(kBlack);                                                                                                                  
	likelihd_limit_c->SetLineWidth(2);                                                                                                                       
	likelihd_limit_c->SetLineStyle(2);

	TGraphAsymmErrors *likelihd_limit_1sigma = new TGraphAsymmErrors(nmass,Mass,Median,Zero,Zero,N1S,P1S);                                                                                                  
	likelihd_limit_1sigma->SetFillColor(kGreen);

	TGraphAsymmErrors *likelihd_limit_2sigma = new TGraphAsymmErrors(nmass,Mass,Median,Zero,Zero,N2S,P2S);                                                                                                  
	likelihd_limit_2sigma->SetFillColor(kYellow);   

	TMultiGraph *likelihd_limit = new TMultiGraph("exclusionlimit_p",";m_{H} [GeV]; 95% CL limit on #sigma/#sigma_{SM}");  
	likelihd_limit->Add(likelihd_limit_2sigma,"E3");                                                                                                         
	likelihd_limit->Add(likelihd_limit_1sigma,"E3");                                                                                                         
	likelihd_limit->Add(likelihd_limit_c, "L");                                                                                                              
    if(showobs){ likelihd_limit->Add(likelihd_limit_d, "LP"); }

	TCanvas *c1 = new TCanvas("Limits_Canvas","Limits Canvas",600,600);
	gStyle->SetOptStat(0);
	//c1->SetLogy () ;
	//c1->SetGridy(1);
	//c1->SetGrid(1);

	Double_t y_max_h2=TMath::Max(P2S[nmass-1]+Median[nmass-1]+1, Observed[nmass-1]+1);
	y_max_h2=10;
	cout<<P2S[nmass-1]+Median[nmass-1]+1<<","<<Observed[nmass-1]+1<<endl;

	//TH2F *h2=new TH2F("h2",";m_{H} (GeV/c^{2});95% CL limit on #sigma/#sigma_{SM}",100,Mass[0]-30,Mass[nmass-1]+30,20,0, y_max_h2);
	TH2F *h2=new TH2F("h2",";Higgs boson mass (GeV/c^{2});95% CL limit on #sigma/#sigma_{SM}",100,Mass[0]-1 ,Mass[nmass-1]+1 ,20,0, y_max_h2);
	h2->Draw();

	likelihd_limit->Draw("");	

	double min_x = h2->GetXaxis ()->GetXmin () ;
	double max_x = h2->GetXaxis ()->GetXmax () ;
	double min_y = h2->GetYaxis ()->GetXmin () ; 
	double max_y = h2->GetYaxis ()->GetXmax () ; 

	TLine *line = new TLine(min_x,1,max_x,1);
	line->SetLineColor(kRed); line->SetLineWidth(2);
	line->Draw();

	//draw grid on top of limits
	TH1D* postGrid=new TH1D("postGrid","postGrid",1,600,1000);
	postGrid->GetYaxis()->SetRangeUser(0.,10);
	postGrid->Draw("AXISSAME");
	//postGrid->Draw("AXIGSAME");
	line->Draw();

	TLegend * theLeg = new TLegend(0.25, 0.57, 0.62, 0.87, "", "NDC");
	theLeg.SetName("theLegend");
	theLeg.SetBorderSize(0);
	theLeg.SetLineColor(0);
	theLeg.SetFillColor(0);
	//theLeg.SetFillStyle(1001);
	theLeg.SetFillStyle(0);
	theLeg.SetLineWidth(0);
	theLeg.SetLineStyle(0);
	theLeg.SetTextFont(42);
	theLeg.SetTextSize(.030);

	theLeg->AddEntry(likelihd_limit_d, "95% C.L.Observed Limit","l");
	theLeg->AddEntry(likelihd_limit_c, "95% C.L.Expected Limit","l");
	theLeg->AddEntry(likelihd_limit_1sigma, "#pm1 #sigma Expected Limit","f");
	theLeg->AddEntry(likelihd_limit_2sigma, "#pm2 #sigma Expected Limit","f");
	theLeg->AddEntry(line, "SM Expected","l");
	theLeg->Draw();

	//banner = TLatex(0.18,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s}=8TeV %s+jets"%(self.GetLumi(),self.channel)));
	float tmp_lumi=0.;
	double tmp_banner_offset=0.;
	TString tmp_channel;
	if(channel=="mu"){ tmp_lumi=mu_lumi; tmp_channel="W#rightarrow #mu #nu";}
	if(channel=="em"){ tmp_lumi=(mu_lumi+el_lumi)/2.; tmp_channel="e+#mu";tmp_banner_offset=0.08;}
	if(channel=="el"){ tmp_lumi=el_lumi; tmp_channel="W#rightarrow e #nu";}
	//banner = TLatex(0.18,0.96,Form("CMS Preliminary, %0.1f fb^{-1} at #sqrt{s} = 8 TeV, %s", (tmp_lumi), tmp_channel.Data() ) );
	banner = TLatex(0.33+tmp_banner_offset,0.96,Form("CMS Preliminary, %0.1f fb^{-1} at #sqrt{s} = 8 TeV, %s", (tmp_lumi), tmp_channel.Data() ) );
	banner.SetNDC(); banner.SetTextSize(0.028);
	banner.Draw();

	//c1->SetLogy(0);
	c1->Update();
	c1->Print(Form("CL_%s_%s.pdf",channel,model),"pdf");
	c1->Print(Form("CL_%s_%s.png",channel,model),"png");
}
