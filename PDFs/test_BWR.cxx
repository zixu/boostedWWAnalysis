/*
 * =====================================================================================
 *
 *       Filename:  test.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/11/2013 02:59:28 PM CST
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
	gROOT->ProcessLine(".L RooBWRunPdf.cxx+");
	gROOT->ProcessLine(".L RooSeymourPdf.cxx+");

	RooRealVar x("x","x",200,2200);
	RooRealVar mean("mean","mean",1000);
	RooRealVar width("width","width",600,0,1000);
	RooBWRunPdf bw("bw","bw",x,mean,width);
	//RooSeymourPdf bw("bw","bw",x,mean,width);
	RooPlot* xframe=x.frame();
	bw.plotOn(xframe);
	xframe.Draw();
}
