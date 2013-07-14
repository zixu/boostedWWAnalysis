/*
 * =====================================================================================
 *
 *       Filename:  test2Gaussian.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/22/2010 11:07:09 AM CST
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
	gROOT->ProcessLine(".L RooErfExpPdf.cxx+");
	gROOT->ProcessLine(".L RooErfExp_Gaus_Pdf.cxx+");
	gROOT->ProcessLine(".L RooErfExp_2Gaus_Pdf.cxx+");
	
/*	RooRealVar x("x","x",90,260);
	x.setBins(100000);
	RooRealVar mean("mean","mean",100);
	RooRealVar mean0("mean0","mean0",0);
	RooRealVar sigma("sigma","sigma",0.3);
	RooRealVar c0("c0","c0",0.1,0,10);
	RooRealVar c1("c1","c1",1.5,-10,10);
	RooRealVar c2("c2","c2",1.75,0,10);
	Roo2Gaussian g("g","g",x,mean0,sigma,c0,c1,c2);
	xframe=x.frame();
	//g.plotOn(xframe);
	RooRealVar width("width","width",0.1);
	RooBreitWigner sourse("sourse","sourse",x,mean,width);
	RooFFTConvPdf sig("sig","sig",x,sourse,g);
	sig.plotOn(xframe);
*/
	RooRealVar x("x","x",30,200);
	//RooRealVar x("x","x",150,180);
	x.setBins(10000);
	RooRealVar mean("mean","mean",80);
	RooRealVar mean2("mean2","mean2",170);
	RooRealVar sigma("sigma","sigma",10);
	RooRealVar sigma2("sigma2","sigma2",20);
	RooRealVar high("high","high",1);
	RooRealVar high2("high2","high2",1);
	RooRealVar c("c","c",-0.05);
	RooRealVar off("off","off",100,10,210);
	RooRealVar width("width","width",30);
	RooErfExpPdf g_erfexp("g_erfexp","g_erfexp",x,c,off,width);
	RooErfExp_Gaus_Pdf g_erfexp_gaus("g_erfexp_gaus","g_erfexp_gaus",x,c,off,width,mean,sigma,high);
	RooErfExp_2Gaus_Pdf g_erfexp_2gaus("g_erfexp_2gaus","g_erfexp_2gaus",x,c,off,width,mean,sigma,high,mean2,sigma2,high2);
	//RooRealVar aera("aera","aera",3,0,10);
	//RooExtendPdf g_area("g_area","g_area",g,aera);
	xframe=x.frame();
	//g_area.plotOn(xframe);
    g_erfexp.plotOn(xframe);
    g_erfexp_gaus.plotOn(xframe);
    g_erfexp_2gaus.plotOn(xframe);
	xframe.Draw();
}

