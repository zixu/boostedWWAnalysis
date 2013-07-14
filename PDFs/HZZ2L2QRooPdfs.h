#ifndef HZZ2L2QROOPDFS
#define HZZ2L2QROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"


double ModMoyal(double x, double mu, double s1, double s2);

class RooModMoyal : 
	public RooAbsPdf {
		public:
			RooModMoyal();
			RooModMoyal(const char *name, const char *title,
						RooAbsReal& _x,
						RooAbsReal& _mu,
						RooAbsReal& _s1,
						RooAbsReal& _s2
					   );
			RooModMoyal(const RooModMoyal& other, const char* name=0) ;
			virtual TObject* clone(const char* newname) const { return new RooModMoyal(*this,newname); }
			inline virtual ~RooModMoyal() { }

		protected:

			RooRealProxy x ;
			RooRealProxy mu;
			RooRealProxy s1;
			RooRealProxy s2;

			Double_t evaluate() const ;

		private:

			ClassDef(RooModMoyal,1)
	};



class RooAlpha : 
	public RooAbsPdf {
		public:
			RooAlpha();
			RooAlpha(const char *name, const char *title,
						RooAbsReal& _x,
						RooAbsReal& _mu,
						RooAbsReal& _s1,
						RooAbsReal& _s2,
						RooAbsReal& _mua,
						RooAbsReal& _s1a,
						RooAbsReal& _s2a
					);
			RooAlpha(const RooAlpha& other, const char* name=0) ;
			virtual TObject* clone(const char* newname) const { return new RooAlpha(*this,newname); }
			inline virtual ~RooAlpha() { }

		protected:

			RooRealProxy x ;
			RooRealProxy mu;
			RooRealProxy s1;
			RooRealProxy s2;
			RooRealProxy mua;
			RooRealProxy s1a;
			RooRealProxy s2a;

			Double_t evaluate() const ;

		private:

			ClassDef(RooAlpha,1)
	};



class RooCB : public RooAbsPdf {
	public:
		RooCB();
		RooCB(const char *name, const char *title,
					RooAbsReal& _x,
					RooAbsReal& _mean,
					RooAbsReal& _width,
					RooAbsReal& _alpha,
					RooAbsReal& _n,
					RooAbsReal& _theta
			 );
		RooCB(const RooCB& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooCB(*this,newname); }
		inline virtual ~RooCB() { }

	protected:

		RooRealProxy x ;
		RooRealProxy mean;
		RooRealProxy width;
		RooRealProxy alpha;
		RooRealProxy n;
		RooRealProxy theta;

		Double_t evaluate() const ;

	private:

		ClassDef(RooCB,1)
};


class RooDoubleCB : public RooAbsPdf {
	public:
		RooDoubleCB();
		RooDoubleCB(const char *name, const char *title,
					RooAbsReal& _x,
					RooAbsReal& _mean,
					RooAbsReal& _width,
					RooAbsReal& _alpha1,
					RooAbsReal& _n1,
					RooAbsReal& _alpha2,
					RooAbsReal& _n2
				   );
		RooDoubleCB(const RooDoubleCB& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooDoubleCB(*this,newname); }
		inline virtual ~RooDoubleCB() { }

	protected:

		RooRealProxy x ;
		RooRealProxy mean;
		RooRealProxy width;
		RooRealProxy alpha1;
		RooRealProxy n1;
		RooRealProxy alpha2;
		RooRealProxy n2;

		Double_t evaluate() const ;

	private:

		ClassDef(RooDoubleCB,1)
};

class RooDoubleCBAnalytic : public RooAbsPdf {
	public:
		RooDoubleCBAnalytic();
		RooDoubleCBAnalytic(const char *name, const char *title,
					RooAbsReal& _x,
					RooAbsReal& _mean,
					RooAbsReal& _width,
					RooAbsReal& _alpha1,
					RooAbsReal& _n1,
					RooAbsReal& _alpha2,
					RooAbsReal& _n2
					);
		RooDoubleCBAnalytic(const RooDoubleCBAnalytic& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooDoubleCBAnalytic(*this,newname); }
		inline virtual ~RooDoubleCBAnalytic() { }
		Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
		Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

	protected:

		RooRealProxy x ;
		RooRealProxy mean;
		RooRealProxy width;
		RooRealProxy alpha1;
		RooRealProxy n1;
		RooRealProxy alpha2;
		RooRealProxy n2;

		Double_t evaluate() const ;

	private:

		ClassDef(RooDoubleCBAnalytic,1)
};

class RooFermi : public RooAbsPdf {
	public:
		RooFermi();
		RooFermi(const char *name, const char *title,
					RooAbsReal& _x,
					RooAbsReal& _cutOff,
					RooAbsReal& _beta
				);
		RooFermi(const RooFermi& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooFermi(*this,newname); }
		inline virtual ~RooFermi() { }

	protected:

		RooRealProxy x ;
		RooRealProxy cutOff ;
		RooRealProxy beta ;

		Double_t evaluate() const ;

	private:

		ClassDef(RooFermi,1) 
};

class RooRelBW : public RooAbsPdf {
	public:
		RooRelBW();
		RooRelBW(const char *name, const char *title,
					RooAbsReal& _x,
					RooAbsReal& _mean,
					RooAbsReal& _width,
					RooAbsReal& _n
				);
		RooRelBW(const RooRelBW& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooRelBW(*this,newname); }
		inline virtual ~RooRelBW() { }

	protected:

		RooRealProxy x ;
		RooRealProxy mean ;
		RooRealProxy width ;
		RooRealProxy n ;

		Double_t evaluate() const ;

	private:

		ClassDef(RooRelBW,1)
};


class Triangle : public RooAbsPdf {
	public:
		Triangle();
		Triangle(const char *name, const char *title,                
					RooAbsReal& _m,
					RooAbsReal& _start,
					RooAbsReal& _turn,
					RooAbsReal& _stop
				);	

		Triangle(const Triangle& other, const char* name = 0);
		virtual TObject* clone(const char* newname) const { 
			return new Triangle(*this,newname); }

		inline virtual ~Triangle() { }
		Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
		Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

	protected:

		RooRealProxy m;
		RooRealProxy start;
		RooRealProxy turn;
		RooRealProxy stop;

		Double_t evaluate() const;

	private:

		ClassDef(Triangle,1)
};


class RooLevelledExp : public RooAbsPdf {
	public:
		RooLevelledExp();
		RooLevelledExp(const char *name, const char *title,
					RooAbsReal& _x,
					RooAbsReal& _sigma,
					RooAbsReal& _alpha,
					RooAbsReal& _m,
					//RooAbsReal& _k,
					RooAbsReal& _theta
					);

		RooLevelledExp(const RooLevelledExp& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooLevelledExp(*this,newname); }
		inline virtual ~RooLevelledExp() { }

	protected:

		RooRealProxy x ;
		RooRealProxy sigma;
		RooRealProxy alpha;
		RooRealProxy m;
		// RooRealProxy k;
		RooRealProxy theta;


		Double_t evaluate() const ;

	private:

		ClassDef(RooLevelledExp,1)
};

class RooBkg : public RooAbsPdf {
	public:
		RooBkg();
		RooBkg(const char *name, const char *title,
					RooAbsReal& _ns,
					RooAbsReal& _nsb,
					RooAbsPdf* _pdfSigD,
					RooAbsPdf* _pdfSigMC,
					RooAbsPdf* _pdfSbMC
			  );

		RooBkg(const RooBkg& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooBkg(*this,newname); }
		inline virtual ~RooBkg() { }

	protected:

		RooRealProxy ns ;
		RooRealProxy nsb ;
		RooAbsPdf* pdfSig;
		RooAbsPdf* pdfSigMC;
		RooAbsPdf* pdfSbMC;


		Double_t evaluate() const ;

	private:

		ClassDef(RooBkg,1)
};


#endif
