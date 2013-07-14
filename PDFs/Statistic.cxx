/*
 * =====================================================================================
 *
 *       Filename:  Statistic.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/05/2013 11:32:54 PM CST
 *       Revision:  none
 *       Compiler:  gcc, root
 *
 *         Author:  Zijun Xu, xuzijun123@gmail.com
 *        Company:  School of Physics, Peking Univ.
 *
 * =====================================================================================
 */
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

// from Kalanand Code
double BinP(int N, double p, int x1, int x2) {
   double q=p/(1-p); 
   int k=0; 
   double v = 1; 
   double s=0; 
   double tot=0.0;
    while(k<=N) {
       tot=tot+v;
       if(k>=x1 & k<=x2) { s=s+v; }
       if(tot>1e30){s=s/1e30; tot=tot/1e30; v=v/1e30;}
       k=k+1; 
       v=v*q*(N+1-k)/k;
    }
    return s/tot;
}

void ClopperPearsonLimits(double numerator, double denominator, 
double lowerLimit, double upperLimit, const double CL_low=1.0, 
const double CL_high=1.0) 
{  
//Confidence intervals are in the units of \sigma.

   double ratio = numerator/denominator;
   
// first get the lower limit
   if(numerator==0)   lowerLimit = 0.0; 
   else { 
      double v=ratio/2; 
      double vsL=0; 
      double vsH=ratio; 
      double p=CL_low/100;
      while((vsH-vsL)>1e-5) { 
         if(BinP(denominator,v,numerator,denominator)>p) 
         { vsH=v; v=(vsL+v)/2; } 
         else { vsL=v; v=(v+vsH)/2; } 
      }
      lowerLimit = v; 
   }
   
// now get the upper limit
   if(numerator==denominator) upperLimit = 1.0;
   else { 
      double v=(1+ratio)/2; 
      double vsL=ratio; 
      double vsH=1; 
      double p=CL_high/100;
      while((vsH-vsL)>1e-5) { 
         if(BinP(denominator,v,0,numerator)<p) { vsH=v; v=(vsL+v)/2; } 
         else { vsL=v; v=(v+vsH)/2; } 
      }
      upperLimit = v;
   }
   cout<<"ratio=("<< lowerLimit<<" , "<<ratio<<" , "<<upperLimit<<")\n";
   cout<<"error=("<<ratio-lowerLimit<<" , "<<upperLimit-ratio<<") ="<<(upperLimit-lowerLimit)/2<<endl;
}
