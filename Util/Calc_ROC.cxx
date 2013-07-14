#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <vector> 
#include <iostream>
#include <fstream>


using namespace std;

class Cut_Eff_Calculator{
	public:
		vector<TString> cutname_vect;
		Int_t number_cut;
		Double_t cut_min;
		Double_t cut_max;
		Double_t cut_width;
		Double_t cut_halfwidth;
		vector<double> cut_vect;
		vector<double> cut_eff_vect;
		vector<double> cut_eff_err_vect;
		Double_t eventnumber_before_cut;
		Double_t eventnumber_after_cut;

		Cut_Eff_Calculator(TString cutname_in, Int_t number_cut_in, Double_t cut_min_in, Double_t cut_max_in){
			cutname_vect.push_back(cutname_in);
			number_cut=number_cut_in;
			cut_min=cut_min_in;
			cut_max=cut_max_in;
			cut_width=(cut_max-cut_min)/number_cut;
			cut_halfwidth=cut_width/2.;
			eventnumber_before_cut=0.;
			eventnumber_after_cut =0.;
		};
		Cut_Eff_Calculator( Int_t number_cut_in ){
			number_cut=number_cut_in;
			eventnumber_before_cut=0.;
			eventnumber_after_cut =0.;
		};
		~Cut_Eff_Calculator(){};
};

void Calc_ROC(char* channel){

	cout<<Form("../trainingtrees_19_Jan22/trainingtrees_%s/ofile_pseudodata.root",channel)<<endl;
	TString mcfile = Form("../trainingtrees_%s/ofile_pseudodata.root",channel);
	TString datafile = Form("../trainingtrees_%s/ofile_data.root",channel);
	TChain * chain_data = new TChain("otree"); chain_data->Add(datafile);
	TChain * chain_mc   = new TChain("otree"); chain_mc->Add(  mcfile);

	TCanvas * c1 =new TCanvas();

	//Cut_Eff_Calculator cec_tau2tau1_data("jet_tau2tau1", 20,0. , 1.);
	//Cut_Eff_Calculator cec_tau2tau1_mc  ("jet_tau2tau1", 20,0. , 1.);
	Cut_Eff_Calculator cec_tau2tau1_data("tau2tau1_ttbar", 20,0. , 1.);
	Cut_Eff_Calculator cec_tau2tau1_mc  ("tau2tau1_ttbar", 20,0. , 1.);

	double efficiency_benchmark = 0.50;

    double self_nPV_min=0;
    double self_nPV_max=100;
    double self_deltaPhi_METj_cut=2.0;
    double self_pfMET_cut=50;
    double self_lpt_cut=30;
    if(TString(channel).Contains("el")){
        self_pfMET_cut=70; 
        self_lpt_cut=35;
        }


	//TString precut =Form(" ungroomed_jet_pt > 200.  &&  jet_mass_pr >= 70 &&  jet_mass_pr<= 100 &&  nbjets >=1 &&  mass_lvj >= 400 &&  mass_lvj<=1400 &&   nPV >=%g &&  nPV<=%g &&  deltaphi_METca8jet>%g &&  mvaMET>%g ",self_nPV_min ,self_nPV_max,self_deltaPhi_METj_cut,self_pfMET_cut);
	TString precut =Form( "jet_mass_pr>=70 && jet_mass_pr<=100 && mass_lvj>0 && isttbar > 0 && l_pt >= %g && pfMET > %g && jet_pt_ttbar > 200 && v_pt > 200", self_lpt_cut, self_pfMET_cut );

	TString weight_data="1";
	TString weight_mc="totalEventWeight*interference_Weight_H600";

	TH1F * temp_hist = new TH1F ("temp_hist","temp_hist",100,0,99999);

	chain_data->Draw("jet_mass_pr>>temp_hist",Form("(%s)*(%s)",weight_data.Data(),precut.Data()));
	double dataall = (double)temp_hist->Integral(); temp_hist->Reset();

	chain_mc->Draw("jet_mass_pr>>temp_hist",Form("(%s)*(%s)",weight_mc.Data(),precut.Data()));
	double mcall = (double)temp_hist->Integral(); temp_hist->Reset();

	for(double cut_iter =cec_tau2tau1_data.cut_min; cut_iter<cec_tau2tau1_data.cut_max; cut_iter=cut_iter+cec_tau2tau1_data.cut_width)
	{	
		TString cut =precut;

		cut+=Form(" && %s<", cec_tau2tau1_data.cutname_vect.at(0).Data()); cut+=cut_iter;

		chain_data->Draw("jet_mass_pr>>temp_hist",Form("(%s)*(%s)",weight_data.Data(),cut.Data()));
		double datapass = (double)temp_hist->Integral();
		double dataeff_ = datapass/dataall;
		double dataefferror_ = dataeff_ * sqrt(1./datapass+1./dataall);
        if(datapass==0 || dataall ==0){ dataeff_=0;dataefferror_=0; }
		temp_hist->Reset();
		cout<<dataall<<","<<datapass<<","<<dataeff_<<","<<dataefferror_<<"  || "<<Form("(%s)*(%s)",weight_data.Data(),cut.Data())<<endl;

		cec_tau2tau1_data.cut_vect.push_back(cut_iter);
		cec_tau2tau1_data.cut_eff_vect.push_back(dataeff_);
		cec_tau2tau1_data.cut_eff_err_vect.push_back(dataefferror_);
	}


	for(double cut_iter =cec_tau2tau1_mc.cut_min; cut_iter<cec_tau2tau1_mc.cut_max; cut_iter=cut_iter+cec_tau2tau1_mc.cut_width)
	{	
		TString cut =precut;

		cut+=Form(" && %s<", cec_tau2tau1_mc.cutname_vect.at(0).Data()); cut+=cut_iter;

		chain_mc  ->Draw("jet_mass_pr>>temp_hist",Form("(%s)*(%s)",weight_mc.Data(),cut.Data()));
		double mcpass = (double)temp_hist->Integral();
        //cout<<"mcall="<<mcall<<","<<mcpass<<endl;
		double mceff_ = mcpass/mcall;
		double mcefferror_ = mceff_ * sqrt(1./mcpass+1./mcall);
        if(mcpass==0 || mcall ==0){ mceff_=0;mcefferror_=0; }
		temp_hist->Reset();
		//cout<<mcall<<","<<mcpass<<","<<mceff_<<","<<mcefferror_<<"  || "<<Form("(%s)*(%s)",weight_mc.Data(),cut.Data())<<endl;

		cec_tau2tau1_mc.cut_vect.push_back(cut_iter);
		cec_tau2tau1_mc.cut_eff_vect.push_back(mceff_);
		cec_tau2tau1_mc.cut_eff_err_vect.push_back(mcefferror_);
	}

	c1->SetGridx(1);
	c1->SetGridy(1);
	Double_t *x =new Double_t[cec_tau2tau1_data.number_cut];
	Double_t *xe=new Double_t[cec_tau2tau1_data.number_cut];
	Double_t *y =new Double_t[cec_tau2tau1_data.number_cut];
	Double_t *ye=new Double_t[cec_tau2tau1_data.number_cut];
	
	Double_t *mc_x =new Double_t[cec_tau2tau1_mc.number_cut];
	Double_t *mc_xe=new Double_t[cec_tau2tau1_mc.number_cut];
	Double_t *mc_y =new Double_t[cec_tau2tau1_mc.number_cut];
	Double_t *mc_ye=new Double_t[cec_tau2tau1_mc.number_cut];
	for(Int_t i=0;i<cec_tau2tau1_data.number_cut;i++){
		x[i]=cec_tau2tau1_data.cut_vect.at(i);
		xe[i]=0;
		y[i]=cec_tau2tau1_data.cut_eff_vect.at(i);
		ye[i]=cec_tau2tau1_data.cut_eff_err_vect.at(i);
	}
	TGraphErrors *gr = new TGraphErrors(cec_tau2tau1_data.number_cut,x,y,xe,ye);
	//TGraph *gr = new TGraph(n,x,y);

	for(Int_t i=0;i<cec_tau2tau1_mc.number_cut;i++){
		mc_x[i]=cec_tau2tau1_mc.cut_vect.at(i);
		mc_xe[i]=0;
		mc_y[i]=cec_tau2tau1_mc.cut_eff_vect.at(i);
		mc_ye[i]=cec_tau2tau1_mc.cut_eff_err_vect.at(i);
        //cout<<mc_x[i]<<","<<mc_xe[i]<<","<<mc_y[i]<<","<<mc_ye[i]<<endl;
	}
	TGraphErrors *gr_mc = new TGraphErrors(cec_tau2tau1_mc.number_cut,mc_x,mc_y,mc_xe,mc_ye);
	
	TString grname = "Cut on tau2/tau1 ";
	gr->SetTitle(grname);
	gr->GetXaxis()->SetTitle("tau2/tau1");
	gr->GetYaxis()->SetTitle("efficiency");
	gr->GetXaxis()->SetRangeUser(0,1.2);
	gr->GetYaxis()->SetRangeUser(0,1.2);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.5);
	gr->SetLineColor(2);
    gr->SetLineWidth(2);
    gr->SetLineColor(2); gr->SetMarkerColor(1); gr->SetFillColor(2);
	gr->Draw("AC");
    gr_mc->SetLineColor(4); gr_mc->SetMarkerColor(1); gr_mc->SetFillColor(4);
    gr_mc->SetLineWidth(2);
    gr_mc->Draw("C");
    TLegend *lg=new TLegend(0.1,0.7,0.48,0.9);
    if(TString(channel).Contains("el")){ 
       lg->AddEntry(gr,"2012 Data: 19.2/fb el","le");
    }else{ 
       lg->AddEntry(gr,"2012 Data: 19.3/fb mu","le");
    }
    lg->AddEntry(gr_mc,"MC:TTbar+SingleTop+VV+WJets","le");
    lg->Draw();
	//c1->Print("data_eff_vs_tau2tau1.eps","eps");
	c1->Print(Form("totalMC_eff_vs_tau2tau1_%s.eps",channel),"eps");
	c1->Print(Form("totalMC_eff_vs_tau2tau1_%s.png",channel),"png");


	/*ofstream outFile("600_qjet_tau21.txt");
	outFile<<"Looking for point : ("<<signal_efficiency_benchmark<<","<<bkg_reject_benchmark<<")"<<endl;
	for(int k =0; k<n; k++)
	{
		if(fabs(signal_efficiency_benchmark-x[k])<0.05&&fabs(bkg_reject_benchmark-y[k])<0.05)
		{
			outFile<<"point "<<k<<": "<<endl;
			outFile<<"precut:  "<<precut<<endl;
			if(cutmassdrop)outFile<<"massdrop < "<<massdropcut_vect.at(k)<<endl;
			if(cutqjet)outFile<<"qjet < "<<qjetcut_vect.at(k)<<endl;
			if(cuttau21)outFile<<"tau2/tau1 < "<<tau21cut_vect.at(k)<<endl;
			outFile<<"signal efficiency = "<<x[k]<<" +- "<<xe[k]<<endl;
			outFile<<"background rejection = "<<y[k]<<" +- "<<ye[k]<<endl<<endl;
		}
	}*/
}
