/*
 * =====================================================================================
 *
 *       Filename:  test.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/26/2012 03:50:37 PM CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * =====================================================================================
 */

{
	using namespace RooFit;
	gROOT.ProcessLine(".L Util.cxx++");
	gROOT.ProcessLine(".L PdfDiagonalizer.cc+");
	RooRealVar x("x","x",-10,10,"GeV/c^{2}") ;
	//x.setBins(100);
	TTree* tree = makeTTree() ;
	RooDataSet ds("ds","ds",RooArgSet(x),Import(*tree)) ;
	RooRealVar mean("mean","mean",0,-10,10) ;
	RooRealVar sigma("sigma","sigma",3,0.1,10) ;
	RooRealVar number("number","number",100,0,10000);
	number.setError(10);
	RooGaussian gauss("gauss","gauss",x,mean,sigma) ;
	RooExtendPdf gauss_ex("gauss_ex","gauss_ex",gauss, number);
	//RooFitResult *rfres=gauss_ex.fitTo(ds, Save(1)) ;
	RooFitResult *rfres=gauss.fitTo(ds, Save(1)) ;
	rfres.Print("v");
	rfres.covarianceMatrix().Print();

	RooPlot*mplot=x.frame();
    //draw_error_band(ds,gauss_ex,number,rfres,mplot);

	ds.plotOn(mplot);
	gauss.plotOn(mplot,RooFit::VisualizeError(*rfres,1,kFALSE),DrawOption("F"),RooFit::FillColor(kOrange));
	gauss.plotOn(mplot);
	draw_error_band(ds,gauss,number,rfres,mplot);
	//gauss.plotOn(mplot);
	ds.plotOn(mplot);


	TCanvas *c1=new TCanvas("c1","c1");
	mplot->Draw();

	cout<<endl<<"Next Section******************************+++++++++++++++++++++++++++++++++++++++++++++++++++++*"<<endl<<endl;
	RooWorkspace wssig("wssig");
	PdfDiagonalizer DecorSig("DecorSig",&wssig,*rfres);
	RooAbsPdf* sigfDec=DecorSig.diagonalize(gauss);
	wssig.Print();

	RooWorkspace wsDecor("wsDecor");

	sigfDec->Print("v");
	wsDecor.Print();

	wsDecor.import(*sigfDec);
	RooFitResult *rfresDec=sigfDec.fitTo(ds, Save(1)) ;
	//rfresDec.Print("v");
	//rfresDec.covarianceMatrix().Print();


	RooPlot*mplot2=x.frame();
	ds.plotOn(mplot2);
	sigfDec.plotOn(mplot2,RooFit::VisualizeError(*rfresDec,1,kFALSE),DrawOption("F"),RooFit::FillColor(kOrange));
	//draw_error_band_Decor( wsDecor,number,mplot2);
	//vector<string> vect_para_names;
	//vect_para_names.push_back("DecorSig_eig0");
	//vect_para_names.push_back("DecorSig_eig1");
	RooArgList paras(*(wsDecor.var("DecorSig_eig0")),*(wsDecor.var("DecorSig_eig1")));
	draw_error_band_Decor("gauss_DecorSig","x", paras, wsDecor,number,mplot2);
	sigfDec.plotOn(mplot2);
	ds.plotOn(mplot2);

	TCanvas *c2=new TCanvas("c2","c2");
	mplot2->Draw();
	


}
