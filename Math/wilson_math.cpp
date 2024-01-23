#include <cmath>
#include "Math.h"

constexpr double pi = 3.141592654;
/*--------------------------------------------------------------------*/

double Bplus(double x, double y)
{
	return y/(x-y)*(log(y)/(y-1.)-log(x)/(x-1.));
}

/*---------------------------------------------------------------------*/

double D2(double x, double y)
{
	if(fabs(x-y)<1.e-5)
	{
		if(fabs(x-1.)<1.e-5) return -0.5;	
		return (1.-x+log(x))/(1.-x)/(1.-x);
	}

	return (D3(x)-D3(y))/(x-y);
}

/*---------------------------------------------------------------------*/

double D3(double x)
{
	if(fabs(x)<1.e-5) return 0.;
	if(fabs(x-1.)<1.e-5) return -1.;
	
	return x*log(x)/(1.-x);
}

/*--------------------------------------------------------------------*/

double h10(double x)
{
	if(fabs(1.-x)<1.e-5) return h10(0.9999);
	
	return (3.*x*x-2.*x)/3./pow(x-1.,4.)*log(x) 
	-(8.*x*x+5.*x-7.)/18./pow(x-1.,3.);
}

/*--------------------------------------------------------------------*/

double h20(double x)
{
	if(fabs(1.-x)<1.e-5) return h20(0.9999);
	
	return (-6.*x*x+4.*x)/3./pow(x-1.,3.)*log(x) 
	+(7.*x-5.)/3./pow(x-1.,2.);
}

/*--------------------------------------------------------------------*/

double h30(double x)
{
	if(fabs(1.-x)<1.e-5) return h30(0.9999);
	
	return (-6.*x*x*x+9.*x*x-2)/9./pow(x-1.,4.)*log(x) 
	+(52.*x*x-101.*x+43.)/54./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h40(double x)
{
	if(fabs(1.-x)<1.e-5) return h40(0.9999);
	
	return -log(x)/3./pow(x-1.,4.)+(2.*x*x-7.*x+11.)/18./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h50(double x)
{
	if(fabs(1.-x)<1.e-5) return h50(0.9999);
	
	return -x/pow(x-1.,4.)*log(x)+(-x*x+5.*x+2.)/6./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h60(double x)
{
	if(fabs(1.-x)<1.e-5) return h60(0.9999);

	return 2.*x/pow(x-1.,3.)*log(x)-(x+1.)/pow(x-1.,2.);
}

/*----------------------------------------------------------------------*/

double f20(double x)
{
	if(fabs(1.-x)<1.e-5) return f20(0.9999);

	return -x/(x-1.)*(1.-1./(x-1.)*log(x));
}

/*----------------------------------------------------------------------*/

double f30(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f30(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f30(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f30(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f30(y*0.9998,y);

	return x*log(x)/(x-1.)/(x-y)+y*log(y)/(y-1.)/(y-x);
}

/*----------------------------------------------------------------------*/

double f40(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f40(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f40(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f40(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f40(y*0.9998,y);
	
	return x*x*log(x)/(x-1.)/(x-y)+y*y*log(y)/(y-1.)/(y-x);
}

/*----------------------------------------------------------------------*/

double f50(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f50(0.9996,0.9998,1.0002);

	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f50(0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f50(0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f50(x,0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f50(0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f50(x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f50(x,y,0.9999);
 
	if(fabs(1.-x/y)<1.e-5) return f50(y*0.9998,y,z);
	if(fabs(1.-x/y)<1.e-5) return f50(y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f50(x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f50(x,y,x*0.9998);

	return x*x*log(x)/(x-1.)/(x-y)/(x-z)+y*y*log(y)/(y-1.)/(y-x)/(y-z)
	+z*z*log(z)/(z-1.)/(z-x)/(z-y);
}

/*----------------------------------------------------------------------*/

double f60(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f60(0.9996,0.9998,1.0002);

	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f60(0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f60(0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f60(x,0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f60(0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f60(x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f60(x,y,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f60(y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f60(x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f60(x,y,x*0.9998); 
	
	return x*log(x)/(x-1.)/(x-y)/(x-z)+y*log(y)/(y-1.)/(y-x)/(y-z)
	+z*log(z)/(z-1.)/(z-x)/(z-y);
}

/*----------------------------------------------------------------------*/

double f70(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f70(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f70(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f70(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f70(y*0.9998,y);
	
	return x*log(x)/(x-1.)/(x-y)+x*log(y)/(y-1.)/(y-x);
}

/*----------------------------------------------------------------------*/

double f80(double x)
{
	if(fabs(x-1.)<1.e-5) return 1.;
	
	return x*log(x)/(x-1.);
}

/*----------------------------------------------------------------------*/

double f90(double w, double x, double y, double z)
{
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9996,0.9998,1.0002,1.0004);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f90(0.9996,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(w,0.9996,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9996,x,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9996,0.9998,y,1.0002);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)) return f90(0.9998,1.0002,y,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f90(0.9998,x,1.0002,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9998,x,y,1.0002);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f90(w,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(w,0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(w,x,0.9998,1.0002);
	
	if(fabs(1.-w)<1.e-5) return f90(0.9999,x,y,z);
	if(fabs(1.-x)<1.e-5) return f90(w,0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f90(w,x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f90(w,x,y,0.9999);

	if(fabs(1.-w/x)<1.e-5) return f90(x*0.9998,x,y,z);
	if(fabs(1.-w/y)<1.e-5) return f90(y*0.9998,x,y,z);
	if(fabs(1.-w/z)<1.e-5) return f90(z*0.9998,x,y,z);
	if(fabs(1.-x/y)<1.e-5) return f90(w,y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f90(w,x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f90(w,x,y,x*0.9998); 


	return w*w*log(w)/(w-1.)/(w-x)/(w-y)/(w-z) +x*x*log(x)/(x-1.)/(x-w)/(x-y)/(x-z) 
	+y*y*log(y)/(y-1.)/(y-x)/(y-w)/(y-z) +z*z*log(z)/(z-1.)/(z-x)/(z-y)/(z-w);
}

/*----------------------------------------------------------------------*/

double f100(double w, double x, double y, double z)
{
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9996,0.9998,1.0002,1.0004);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f100(0.9996,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(w,0.9996,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9996,x,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9996,0.9998,y,1.0002);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)) return f100(0.9998,1.0002,y,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f100(0.9998,x,1.0002,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9998,x,y,1.0002);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f100(w,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(w,0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(w,x,0.9998,1.0002);

	if(fabs(1.-w)<1.e-5) return f100(0.9999,x,y,z);
	if(fabs(1.-x)<1.e-5) return f100(w,0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f100(w,x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f100(w,x,y,0.9999);

	if(fabs(1.-w/x)<1.e-5) return f100(x*0.9998,x,y,z);
	if(fabs(1.-w/y)<1.e-5) return f100(y*0.9998,x,y,z);
	if(fabs(1.-w/z)<1.e-5) return f100(z*0.9998,x,y,z);
	if(fabs(1.-x/y)<1.e-5) return f100(w,y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f100(w,x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f100(w,x,y,x*0.9998); 
	
	return w*log(w)/(w-1.)/(w-x)/(w-y)/(w-z) +x*log(x)/(x-1.)/(x-w)/(x-y)/(x-z) 
	+y*log(y)/(y-1.)/(y-x)/(y-w)/(y-z) +z*log(z)/(z-1.)/(z-x)/(z-y)/(z-w);
}

/*----------------------------------------------------------------------*/

double f110(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f110(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f110(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f110(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f110(y*0.9998,y);
	
	return x*log(x)/(x-y)+x*log(y)/(y-x);
}

/*----------------------------------------------------------------------*/

double h11(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h11(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h11(0.99,y);
	if(fabs(1.-y)<1.e-3) return h11(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h11(y*0.98,y);

	return ((-48.*x*x*x-104.*x*x+64.*x)*Li2(1.-1./x)
	+(-378.*x*x*x-1566.*x*x+850.*x+86.)/9./(x-1.)*log(x)
	+(2060.*x*x*x+3798.*x*x-2664.*x-170.)/27.
	+((12.*x*x*x-124.*x*x+64.*x)/(x-1.)*log(x)+(-56.*x*x*x+258.*x*x+24.*x-82.)/3.)*y)/9./pow(x-1.,4.);
}

/*----------------------------------------------------------------------*/

double h21(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h21(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h21(0.99,y);
	if(fabs(1.-y)<1.e-3) return h21(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h21(y*0.98,y);

	return ((224.*x*x-96.*x)*Li2(1.-1./x)
	+(-24.*x*x*x+352.*x*x-128.*x-32.)*log(x)/(x-1.)
	+(-340.*x*x+132.*x+40.)
	+((-24.*x*x*x+176.*x*x-80.*x)*log(x)/(x-1.)+(-28.*x*x-108.*x+64.))*y)/9./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h31(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h31(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h31(0.99,y);
	if(fabs(1.-y)<1.e-3) return h31(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h31(y*0.98,y);

	return (32.*x*x*x+120.*x*x-384.*x+128.)*Li2(1.-1./x)/81./pow(x-1.,4.)
	+(-108.*x*x*x*x+1058.*x*x*x-898.*x*x-1098.*x+710.)/81./pow(x-1.,5.)*log(x)
	+(-304.*x*x*x-13686.*x*x+29076.*x-12062.)/729./pow(x-1.,4.)
	+((540.*x*x*x-972.*x*x+232.*x+56.)/81./pow(x-1.,5.)*log(x)
	+(-664.*x*x*x+54.*x*x+1944.*x-902.)/243./pow(x-1.,4.))*y;
}

/*----------------------------------------------------------------------*/

double h41(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h41(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h41(0.99,y);
	if(fabs(1.-y)<1.e-3) return h41(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h41(y*0.98,y);

	return (-562.*x*x*x+1101.*x*x-420.*x+101.)/54./pow(x-1.,4.)*Li2(1.-1./x)
	+(-562.*x*x*x+1604.*x*x-799.*x+429.)/54./pow(x-1.,5.)*log(x)
	+(17470.*x*x*x-47217.*x*x+31098.*x-13447.)/972./pow(x-1.,4.)
	+((89.*x+55.)*log(x)/27./pow(x-1.,5.)+(38.*x*x*x-135.*x*x+54.*x-821.)/162./pow(x-1.,4.))*y;
}

/*----------------------------------------------------------------------*/

double h51(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h51(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h51(0.99,y);
	if(fabs(1.-y)<1.e-3) return h51(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h51(y*0.98,y);

	return ((9.*x*x*x+46.*x*x+49.*x)*Li2(1.-1./x)/2.
	+(81.*x*x*x+594.*x*x+1270.*x+71.)*log(x)/18./(x-1.)
	+(-923.*x*x*x-3042.*x*x-6921.*x-1210.)/108.
	+((10.*x*x+38.*x)/(x-1.)*log(x)+(-7.*x*x*x+30.*x*x-141.*x-26.)/3.)*y)/3./pow(x-1.,4.);
}

/*----------------------------------------------------------------------*/

double h61(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h61(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h61(0.99,y);
	if(fabs(1.-y)<1.e-3) return h61(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h61(y*0.98,y);

	return ((-32.*x*x-24.*x)*Li2(1.-1./x)
	+(-52.*x*x-109.*x-7.)*log(x)/(x-1.)
	+(95.*x*x+180.*x+61.)/2.
	+((-20.*x*x-52.*x)/(x-1.)*log(x)+(-2.*x*x+60.*x+14.))*y)/3./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h71(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h71(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h71(0.99,y);
	if(fabs(1.-y)<1.e-3) return h71(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h71(y*0.98,y);
	
	return (-20.*x*x*x+60.*x*x-60.*x-20.)*Li2(1.-1./x)/27./pow(x-1.,4.)
	+(-60.*x*x+240.*x+4.)/81./pow(x-1.,4.)*log(x)
	+(132.*x*x-382.*x+186.)/81./pow(x-1.,3.)
	+(20.*log(x)/27./pow(x-1.,4.)+(-20.*x*x+70.*x-110.)/81./pow(x-1.,3.))*y;
}

/*----------------------------------------------------------------------*/

double f31(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f31(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f31(0.99,y);
	if(fabs(1.-y)<1.e-3) return f31(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f31(y*0.98,y);
	
	return -28.*y/3./(y-1.)/(x-y)+2.*x*(11.*x+3.*y)/3./(x-1.)/(x-y)/(x-y)*log(x)
	+2.*y*(25.*x-11.*x*y-11.*y-3.*y*y)/3./(y-1.)/(y-1.)/(x-y)/(x-y)*log(y)
	+4.*(1.+y)/(x-1.)/(y-1)*Li2(1.-1./y)+4.*(x+y)/(x-1.)/(x-y)*Li2(1.-x/y);
}

/*----------------------------------------------------------------------*/

double f41(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f41(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f41(0.99,y);
	if(fabs(1.-y)<1.e-3) return f41(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f41(y*0.98,y);
	
	return (59.*x*(1.-y)-y*(59.-3.*y))/6./(y-1.)/(x-y)
	+4.*x*(7.*x*x-3.*x*y+3.*y*y)/3./(x-1.)/(x-y)/(x-y)*log(x)+2.*log(y)*log(y)
	+4.*y*y*(18.*x-11.*x*y-11.*y+4.*y*y)/3./(y-1.)/(y-1.)/(x-y)/(x-y)*log(y)
	+4.*(1.+y*y)/(x-1.)/(y-1)*Li2(1.-1./y)+4.*(x*x+y*y)/(x-1.)/(x-y)*Li2(1.-x/y);
}

/*----------------------------------------------------------------------*/

double f51(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f51(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f51(0.99,y);
	if(fabs(1.-y)<1.e-3) return f51(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f51(y*0.98,y);

	return (-83.-27.*x*(y-1.)+27.*y)/6./(x-1.)/(y-1.)
	
	-4.*x*(1.+x*(12.+y)-y-6.*x*x)/3./(x-1.)/(x-1.)/(x-y)*log(x)
	+2.*(1.+6.*x*x*(y-1.)-3.*x*x*x*(y-1.)+x*(3.*y-4.))/3./(x-1.)/(x-1.)/(y-1.)/(x-y)*log(x)*log(x)
	-4.*y*(3.*x*x*(y-1.)+x*y*(3.-2.*y)+y*y*(y-2.))/3./(x-1.)/(y-1)/(x-y)/(x-y)*Li2(1.-x/y)
	-4.*(1.-3.*x-x*x*(3.-6.*y)-x*x*x)/3./(x-1.)/(y-1)/(x-y)*Li2(1.-1./x)
	
	-4.*y*(1.+y*(12.+x)-x-6.*y*y)/3./(y-1.)/(y-1.)/(y-x)*log(y)
	+2.*(1.+6.*y*y*(x-1.)-3.*y*y*y*(x-1.)+y*(3.*x-4.))/3./(y-1.)/(y-1.)/(x-1.)/(y-x)*log(y)*log(y)
	-4.*x*(3.*y*y*(x-1.)+y*x*(3.-2.*x)+x*x*(x-2.))/3./(y-1.)/(x-1)/(y-x)/(y-x)*Li2(1.-y/x)
	-4.*(1.-3.*y-y*y*(3.-6.*x)-y*y*y)/3./(y-1.)/(x-1)/(y-x)*Li2(1.-1./y)
	+4.*log(x)*(f40(x,y)+(f40(1.0001*x,y)-f40(0.9999*x,y))/0.0002+(f40(x,1.0001*y)-f40(x,0.9999*y))/0.0002);
}

/*----------------------------------------------------------------------*/

double f81(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f81(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f81(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f81(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f81(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f81(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f81(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f81(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f81(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f81(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f81(x,y,x*0.98);

	return -28.*y*y/3./(y-1.)/(x-y)/(y-z)
	+4.*x*(7.*x*x-3.*x*y+3.*y*y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+4.*z*(7.*z*z-3.*z*y+3.*y*y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	-4.*y*y*(x*(4.*y*y+18.*z-11.*y*(1.+z))+y*(3.*y*y-11.*z+4.*y*(1.+z)))/3./(y-1.)/(y-1.)/(x-y)/(x-y)/(y-z)/(y-z)*log(y)
	-4.*(1.+y*y)/(x-1.)/(y-1.)/(z-1.)*Li2(1.-1./y)
	+4.*(x*x+y*y)/(x-1.)/(x-y)/(x-z)*Li2(1.-x/y)
	+4.*(z*z+y*y)/(z-1.)/(z-y)/(z-x)*Li2(1.-z/y);
}

/*----------------------------------------------------------------------*/

double f91(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f91(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f91(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f91(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f91(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f91(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f91(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f91(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f91(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f91(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f91(x,y,x*0.98);
		
	return -28.*y/3./(y-1.)/(x-y)/(y-z)
	+2.*x*(11.*x+3.*y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+2.*z*(11.*z+3.*y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	+2.*y*(x*(3.*y*y-25.*z+11.*y*(1.+x))+y*(-17.*y*y+11.*z+3.*y*(1.+z)))/3./(y-1.)/(y-1.)/(x-y)/(x-y)/(y-z)/(y-z)*log(y)
	-4.*(1.+y)/(x-1.)/(y-1.)/(z-1.)*Li2(1.-1./y)
	+4.*(x+y)/(x-1.)/(x-y)/(x-z)*Li2(1.-x/y)
	+4.*(z+y)/(z-1.)/(z-y)/(z-x)*Li2(1.-z/y);
}

/*----------------------------------------------------------------------*/

double f111(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f111(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f111(0.99,y);
	if(fabs(1.-y)<1.e-3) return f111(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f111(y*0.98,y);
	
	return 4.*x*(8.*y+(x-1.)*(x-y)*pi*pi)/3./y/(x-1.)/(x-y)
	-8.*x*(x*x-7.*y+3.*x*(1.+y))/3./(x-y)/(x-y)/(x-1.)/(x-1.)*log(x)
	-8.*x*(3.*x-7.*y)/3./(x-y)/(x-y)/(y-1.)*log(y)
	-8.*x/(y-1.)*Li2(1.-1./x)
	+8.*x/y/(y-1.)*Li2(1.-y/x);
}

/*----------------------------------------------------------------------*/

double f121(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f121(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f121(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f121(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f121(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f121(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f121(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f121(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f121(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f121(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f121(x,y,x*0.98);
	
	return -28.*y*y/3./(x-y)/(y-1.)/(y-z)
	+4.*x*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+4.*z*z*(6.*z+y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	-4.*y*y*(x*(6.*y*y+20.*z-13.*y*(1.+z))+y*(y*y-13.*z+6.*y*(1.+z)))/3./(x-y)/(x-y)/(y-1.)/(y-1.)/(y-z)/(y-z)*log(y);
	
}

/*----------------------------------------------------------------------*/

double f131(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f131(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f131(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f131(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f131(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f131(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f131(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f131(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f131(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f131(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f131(x,y,x*0.98);
	
	return -28.*y/3./(x-y)/(y-1.)/(y-z)
	+4.*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+4.*z*(6.*z+y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	+4.*y*(x*(y*y-13.*z+6.*y*(1.+z))+y*(y-8.*y*y+6.*z+y*z))/3./(x-y)/(x-y)/(y-1.)/(y-1.)/(y-z)/(y-z)*log(y);
}

/*----------------------------------------------------------------------*/

double f141(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f141(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f141(0.99,y);
	if(fabs(1.-y)<1.e-3) return f141(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f141(y*0.98,y);
	
	return 32.*x*x/3./(x-1.)/(x-y)
	-8.*x*x*(7.*x*(1.+y)-11.*y-3.*x*x)/3./(x-1.)/(x-1.)/(x-y)/(x-y)*log(x)
	-8.*x*y*(3.*x-7.*y)/3./(x-y)/(x-y)/(y-1.)*log(y)
	-8.*x/(y-1.)*Li2(1.-1./x)
	+8.*x/(y-1.)*Li2(1.-y/x);
}

/*----------------------------------------------------------------------*/

double f151(double x)
{
	if(fabs(1.-x)<1.e-3) return f151(0.99);

	return (1.-3.*x)/(x-1.)+2.*x/(x-1.)/(x-1.)*log(x)+2.*x/(x-1.)*Li2(1.-1./x);
}


/*----------------------------------------------------------------------*/

double f161(double x)
{
	if(fabs(1.-x)<1.e-3) return f161(0.99);

	return 28./3./(x-1.)-4.*x*(13.-6.*x)/3./(x-1.)/(x-1.)*log(x);
}

/*----------------------------------------------------------------------*/

double f171(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f171(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f171(0.99,y);
	if(fabs(1.-y)<1.e-3) return f171(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f171(y*0.98,y);
	
	return -28./3./(x-1.)/(y-1.)
	+4.*y*(10.-3.*y)/3./(x-y)/(y-1.)/(y-1.)*log(y)
	-4.*y/(x-y)/(y-1.)/(y-1.)*log(y)*log(y)
	+(4.*(13.*x-6.*x*x-3.*y-7.*x*y+3.*x*x*y)/3./(x-1.)/(x-1.)/(x-y)/(y-1.)+4.*y*log(y)/(x-y)/(y-1.)/(y-1.))*log(x);
}

/*----------------------------------------------------------------------*/

double f181(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f181(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f181(0.99,y);
	if(fabs(1.-y)<1.e-3) return f181(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f181(y*0.98,y);
	
	return -28.*y/3./(x-y)/(y-1.)
	+4.*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)*log(x)
	-4.*y*(y*(6.+y)-x*(13.-6.*y))/3./(x-y)/(x-y)/(y-1.)/(y-1.)*log(y);
}

/*----------------------------------------------------------------------*/

double f191(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f191(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f191(0.99,y);
	if(fabs(1.-y)<1.e-3) return f191(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f191(y*0.98,y);
	
	return -28.*(x*(y-1.)+y)/3./(x-y)/(y-1.)
	+4.*x*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)*log(x)
	+4.*y*y*(x*(20.-13.*y)-y*(13.-6.*y))/3./(x-y)/(x-y)/(y-1.)/(y-1.)*log(y);
}

/*----------------------------------------------------------------------*/

double q11(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q11(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q11(0.99,y);
	if(fabs(1.-y)<1.e-3) return q11(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q11(y*0.98,y);

	return 4./3./(x-y)*(x*x*log(x)/pow(x-1.,4.)-y*y*log(y)/pow(y-1.,4.)) +(4.*x*x*y*y+10.*x*y*y-2.*y*y+10.*x*x*y-44.*x*y+10.*y-2.*x*x+10.*x+4.)/9./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double q21(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q21(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q21(0.99,y);
	if(fabs(1.-y)<1.e-3) return q21(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q21(y*0.98,y);
	
	return 4./3./(x-y)*(x*log(x)/pow(x-1.,4.)-y*log(y)/pow(y-1.,4.)) +(-2.*x*x*y*y+10.*x*y*y+4.*y*y+10.*x*x*y-20.*x*y-14.*y+4.*x*x-14.*x+22.)/9./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double q31(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q31(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q31(0.99,y);
	if(fabs(1.-y)<1.e-3) return q31(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q31(y*0.98,y);

	return 8./3./(x-y)*(-x*x*log(x)/pow(x-1.,3.)+y*y*log(y)/pow(y-1.,3.)) +(-12.*x*y+4.*y+4.*x+4.)/3./pow(x-1.,2.)/pow(y-1.,2.);
}

/*----------------------------------------------------------------------*/

double q41(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q41(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q41(0.99,y);
	if(fabs(1.-y)<1.e-3) return q41(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q41(y*0.98,y);

	return 8./3./(x-y)*(-x*log(x)/pow(x-1.,3.)+y*log(y)/pow(y-1.,3.)) +(-4.*x*y-4.*y-4.*x+12.)/3./pow(x-1.,2.)/pow(y-1.,2.);
}

/*----------------------------------------------------------------------*/

double q51(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q51(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q51(0.99,y);
	if(fabs(1.-y)<1.e-3) return q51(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q51(y*0.98,y);
	
	return 4./27./(x-y)*((6.*x*x*x-9.*x*x+2.)*log(x)/pow(x-1.,4.)-(6.*y*y*y-9.*y*y+2.)*log(y)/pow(y-1.,4.))
	+(104.*x*x*y*y-202.*x*y*y+86.*y*y-202.*x*x*y+380.*x*y-154.*y+86.*x*x-154.*x+56.)/81./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double q61(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q61(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q61(0.99,y);
	if(fabs(1.-y)<1.e-3) return q61(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q61(y*0.98,y);

	return 4./9./(x-y)*(log(x)/pow(x-1.,4.)-log(y)/pow(y-1.,4.))
	+(4.*x*x*y*y-14.*x*y*y+22.*y*y-14.*x*x*y+52.*x*y-62.*y+22.*x*x-62.*x+52.)/27./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double A0t(double x)
{
	return (-46.+159.*x-153.*x*x+22.*x*x*x)/36./pow(1.-x,3.)+x*x*(-3.*x+2.)/2./pow(1.-x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double A1t(double x, double l)
{
	return (32.*pow(x,4.)+244.*pow(x,3.)-160.*x*x+16.*x)/9./pow(1.-x,4.)*Li2(1.-1./x)
	+(-774.*x*x*x*x-2826.*x*x*x+1994.*x*x-130.*x+8.)/81./pow(1.-x,5.)*log(x)
	+(-94.*x*x*x*x-18665.*x*x*x+20682.*x*x-9113.*x+2006.)/243./pow(1.-x,4.)
	+((-12.*x*x*x*x-92.*x*x*x+56.*x*x)/3./pow(1.-x,5.)*log(x)+(-68.*x*x*x*x-202.*x*x*x-804.*x*x+794.*x-152.)/27./pow(1.-x,4.))*l;
}

/*--------------------------------------------------------------------*/

double B0t(double x)
{
	return x/4./pow(1.-x,2.)*log(x)+1./4./(1.-x);
}

/*----------------------------------------------------------------------*/

double C0t(double x)
{
	return (3.*x+2.)*x/8./pow(1.-x,2.)*log(x)+(-x+6.)*x/8./(1.-x);
}

/*----------------------------------------------------------------------*/

double D0t(double x)
{
	return (-3.*pow(x,4.)+30.*pow(x,3.)-54.*x*x+32.*x-8.)/18./pow(1.-x,4.)*log(x)+(-47.*pow(x,3.)+237.*x*x-312.*x+104.)/108./pow(1.-x,3.);
}

/*----------------------------------------------------------------------*/

double B1t(double x, double l)
{
	return -2.*x/(1.-x)/(1.-x)*Li2(1.-1./x) +(-x+17.)*x/3./pow(1.-x,3.)*log(x) +(13.*x+3)/3./(1.-x)/(1.-x) +((2.*x+2)*x/pow(1.-x,3.)*log(x)+4.*x/(1.-x)/(1.-x))*l;
}

/*----------------------------------------------------------------------*/

double C1t(double x, double l)
{
	return (-x*x-4.)*x/(1.-x)/(1.-x)*Li2(1.-1./x) +(3.*x*x+14.*x+23.)*x/3./pow(1.-x,3.)*log(x) +(4.*x*x+7.*x+29.)*x/3./(1.-x)/(1.-x) +((8.*x+2.)*x/pow(1.-x,3.)*log(x)+(x*x+x+8.)*x/(1.-x)/(1.-x))*l;
}

/*----------------------------------------------------------------------*/

double D1t(double x, double l)
{
	return (380.*pow(x,4.)-1352.*pow(x,3.)+1656.*x*x-784.*x+256.)/81./pow(1.-x,4.)*Li2(1.-1./x) +(304.*pow(x,4.)+1716.*pow(x,3.)-4644.*x*x+2768.*x-720.)/81./pow(1.-x,5.)*log(x) +(-6175.*pow(x,4.)+41608.*pow(x,3.)-66723.*x*x+33106.*x-7000.)/729./pow(1.-x,4.) +((648.*pow(x,4.)-720.*pow(x,3.)-232.*x*x-160.*x+32.)/81./pow(1.-x,5.)*log(x)
	+(-352.*pow(x,4.)+4912.*pow(x,3.)-8280.*x*x+3304.*x-880.)/243./pow(1.-x,4.))*l;	
}

/*----------------------------------------------------------------------*/

double F0t(double x)
{
	return  (5.*x*x*x-9.*x*x+30.*x-8.)/12./pow(1.-x,3.)+3.*x*x/2./pow(1.-x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F1t(double x,double l)
{
	return (4.*pow(x,4.)-40.*pow(x,3.)-41.*x*x-x)/3./pow(1.-x,4.)*Li2(1.-1./x)
	+(-144.*x*x*x*x+3177.*x*x*x+3661.*x*x+250.*x-32.)/108./pow(1.-x,5.)*log(x)
	+(-247.*x*x*x*x+11890.*x*x*x+31779.*x*x-2966.*x+1016.)/648./pow(1.-x,4.)
	+((17.*x*x*x+31.*x*x)/pow(1.-x,5.)*log(x)+(-35.*x*x*x*x+170.*x*x*x+447.*x*x+338.*x-56.)/18./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double E0t(double x)
{
	return (-9.*x*x+16.*x-4.)/6./pow(1.-x,4.)*log(x)+(-7.*x*x*x-21.*x*x+42.*x+4.)/36./pow(1.-x,3.);
}

/*----------------------------------------------------------------------*/

double G1t(double x, double l)
{
	return (10.*pow(x,4.)-100.*pow(x,3.)+30.*x*x+160.*x-40.)/27./pow(1.-x,4.)*Li2(1.-1./x)
	+(30.*pow(x,3.)-42.*x*x-332.*x+68.)/81./pow(1.-x,4.)*log(x)
	+(-6.*pow(x,3.)-293.*x*x+161.*x+42.)/81./pow(1.-x,3.)
	+((90.*x*x-160.*x+40.)/27./pow(1.-x,4.)*log(x)+(35.*pow(x,3.)+105.*x*x-210.*x-20.)/81./pow(1.-x,3.))*l;
}

/*----------------------------------------------------------------------*/

double E1t(double x, double l)
{
	return (515.*pow(x,4.)-614.*pow(x,3.)-81.*x*x-190.*x+40.)/54./pow(1.-x,4.)*Li2(1.-1./x)
	+(-1030.*pow(x,4.)+435.*pow(x,3.)+1373.*x*x+1950.*x-424.)/108./pow(1.-x,5.)*log(x)
	+(-29467.*pow(x,4.)+45604.*pow(x,3.)-30237.*x*x+66532.*x-10960.)/1944./pow(1.-x,4.)
	+((-1125.*pow(x,3.)+1685.*x*x+380.*x-76.)/54./pow(1.-x,5.)*log(x)
	+(133.*pow(x,4.)-2758.*pow(x,3.)-2061.*x*x+11522.*x-1652.)/324./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double T(double x)
{
	return -(16.*x+8.)*sqrt(4.*x-1.)*Cl2(2.*asin(0.5/sqrt(x)))+(16.*x+20./3.)*log(x)+32.*x+112./9.;
}

/*----------------------------------------------------------------------*/

double F7_1(double x)
{
	return x*(7.-5.*x-8.*x*x)/24./pow(x-1.,3.)+x*x*(3.*x-2.)/4./pow(x-1.,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F7_2(double x)
{
	return x*(3.-5.*x)/12./pow(x-1.,2.)+x*(3.*x-2.)/6./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_1(double x)
{
	return x*(2.+5.*x-x*x)/8./pow(x-1.,3.)-3.*x*x/4./pow(x-1.,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_2(double x)
{
	return x*(3.-x)/4./pow(x-1.,2.)-x/2./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double H2(double x, double y)
{
	return D2(x,y);
}

/*----------------------------------------------------------------------*/

double B(double m1, double m2, double Q)
{
	double x=pow(m2/m1,2.);

	if(fabs(x-1.)<1.e-5) return -0.5*log(m2*m2/Q/Q);

	return 0.5*(0.5+1./(1.-x)+log(x)/pow(1.-x,2.)-log(m2*m2/Q/Q));
}

/*----------------------------------------------------------------------*/

double G7H(double x, double lu, double ld)
{
	return lu*(ld*4.*x*(4.*(-3.+7.*x-2.*x*x)*Li2(1.-1./x)+(8.-14.*x-3.*x*x)*log(x)*log(x)/(x-1.)
	+2.*(-3.-x+12.*x*x-2.*x*x*x)*log(x)/(x-1.)+3.*(7.-13.*x+2.*x*x))/9./pow((x-1.),3.)
	+lu*2.*x*(x*(18.-37.*x+8.*x*x)*Li2(1.-1./x)+x*(-14.+23.*x+3.*x*x)*log(x)*log(x)/(x-1.)+
	(-50.+251.*x-174.*x*x-192.*x*x*x+21.*x*x*x*x)*log(x)/9./(x-1.)+(797.-5436.*x+7569.*x*x-1202.*x*x*x)/108.)/pow((x-1.),4.)/9.);
}

/*----------------------------------------------------------------------*/
	
double Delta7H(double x, double lu, double ld)
{
	return 2.*x/9./pow((x-1.),4.)*lu*
	(ld*((x-1.)*(21.-47.*x+8.*x*x)+2.*(-8.+14.*x+3.*x*x)*log(x))
	+lu*((-31.-18.*x+135.*x*x-14.*x*x*x)/6.+x*(14.-23.*x-3.*x*x)*log(x)/(x-1.)));
}

/*----------------------------------------------------------------------*/

double G8H(double x, double lu, double ld)
{
	return lu*(ld*x*(0.5*(-36.+25.*x-17.*x*x)*Li2(1.-1./x)+(19.+17.*x)*log(x)*log(x)/(x-1.)
	+0.25*(-3.-187.*x+12.*x*x-14.*x*x*x)*log(x)/(x-1.)+3.*(143.-44.*x+29.*x*x)/8.)/3./pow((x-1.),3.)
	+ lu*x*(x*(30.-17.*x+13.*x*x)*Li2(1.-1./x)-x*(31.+17.*x)*log(x)*log(x)/(x-1.)+
	(-226.+817.*x+1353.*x*x+318.*x*x*x+42.*x*x*x*x)*log(x)/36./(x-1.)+(1130.-18153.*x+7650.*x*x-4451.*x*x*x)/216.)/pow((x-1.),4.)/6.);
}

/*----------------------------------------------------------------------*/

double Delta8H(double x, double lu, double ld)
{
	return (x/6./pow((x-1.),4.))*lu*(ld*((x-1.)*(81.-16.*x+7.*x*x)-2.*(19.+17.*x)*log(x))
	+lu*((-38.-261.*x+18.*x*x-7.*x*x*x)/6.+x*(31.+17.*x)*log(x)/(x-1.)));	
}

/*----------------------------------------------------------------------*/

double EH(double x, double lu)
{
	return lu*lu*x*((x-1.)*(16.-29.*x+7.*x*x)+6.*(3.*x-2.)*log(x))/36./pow((x-1.),4.);
}

/*----------------------------------------------------------------------*/

double G4H(double x, double lu)
{
	return lu*lu*((515.*x*x*x-906.*x*x+99.*x+182.)*x*Li2(1.-1./x)/54./pow(x-1.,4.)
	+(1030.*x*x*x-2763.*x*x-15.*x+980.)*x*log(x)/108./pow(x-1.,5.)
	+(-29467.*x*x*x+68142.*x*x-6717.*x-18134.)*x/1944./pow(x-1.,4.));
}

/*----------------------------------------------------------------------*/

double Delta4H(double x, double lu)
{
	return -lu*lu*((-375.*x*x-95.*x+182.)*x*log(x)/54./pow(x-1.,5.)
	+(133.*x*x*x-108.*x*x+4023.*x-2320.)*x/324./pow(x-1.,4.));
}

/*----------------------------------------------------------------------*/

double G3H(double x, double lu)
{
	return lu*lu*((10.*x*x*x+30.*x-20.)*x*Li2(1.-1./x)/27./pow(x-1.,4.)
	+(30.*x*x-66.*x-56.)*x*log(x)/81./pow(x-1.,4.)
	+(6.*x*x-187.*x+213.)*x/81./pow(x-1.,3.));
}

/*----------------------------------------------------------------------*/

double Delta3H(double x, double lu)
{
	return -lu*lu*((-30.*x+20.)*x*log(x)/27./pow(x-1.,4.)
	+(-35.*x*x+145.*x-80.)*x/81./pow(x-1.,3.));
}

/*----------------------------------------------------------------------*/

double C9llH0(double x, double y, double lu)
{
	return x/y/8.*lu*lu*(-log(y)/(y-1.)+1.)*y*y/(y-1.);
}

/*----------------------------------------------------------------------*/

double D9H0(double x, double lu)
{
	return lu*lu*((-3.*x*x*x+6.*x-4.)*x/18./pow(x-1.,4.)*log(x)+(47.*x*x-79.*x+38.)*x/108./pow(x-1.,3.));	
}

/*----------------------------------------------------------------------*/

double C9llH1(double x, double y, double lu, double L)
{
	if(fabs(y-1.)<1.e-5) return C9llH1(x,0.9999,lu,L);
	
 	return x/y/8.*lu*lu*((-8.*y*y*y+16.*y*y)/pow(y-1.,2.)*Li2(1.-1./y)
	+(-24.*y*y*y+88.*y*y)/3./pow(y-1.,3.)*log(y)
	+(32.*y*y*y-96.*y*y)/3./(y-1.)/(y-1.)
	+(16.*y*y*log(y)/pow(y-1.,3.)+(8.*y*y*y-24.*y*y)/(y-1.)/(y-1.))*L);
}

/*----------------------------------------------------------------------*/

double D9H1(double x, double lu, double L)
{
	if(fabs(x-1.)<1.e-5) return D9H1(0.9999,lu,L);
	
	return lu*lu*((380.*x*x*x-528.*x*x+72.*x+128.)*x/81.*pow(1.-x,4.)*Li2(1.-1./x)
	+(596.*x*x*x-672.*x*x+64.*x+204.)*x*log(x)/81./pow(x-1.,5.)
	+(-6175.*x*x*x+9138.*x*x-3927.*x-764.)*x/729./pow(x-1.,4.)
	+((432.*x*x*x-456.*x*x+40.*x+128.)*x*log(x)/81./pow(x-1.,5.)
	+(-352.*x*x*x-972.*x*x+1944.*x-1052.)*x/243./pow(x-1.,4.))*L);	
}

/*----------------------------------------------------------------------*/

double C7t2mt(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.4) return 12.06+12.93*z+3.013*z*log(z)+96.71*z*z+52.73*z*z*log(z)+147.9*pow(z,3.)
	+187.7*pow(z,3.)*log(z)-144.9*pow(z,4.)+236.1*pow(z,4.)*log(z);
	
	else return 11.74+0.3642*w+0.1155*w*w-0.003145*pow(w,3.)-0.03263*pow(w,4.)-0.03528*pow(w,5.)
	-0.03076*pow(w,6.)-0.02504*pow(w,7.)-0.01985*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C7c2MW(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.4) return 1.525-0.1165*z+0.01975*z*log(z)+0.06283*z*z+0.005349*z*z*log(z)
	+0.01005*pow(z*log(z),2.)-0.04202*pow(z,3.)+0.01535*pow(z,3.)*log(z)-0.00329*z*pow(z*log(z),2.)
	+0.002372*pow(z,4.)-0.0007910*pow(z,4.)*log(z);
		
	else return 1.432+0.06709*w+0.01257*w*w+0.004710*pow(w,3.)+0.002373*pow(w,4.)
	+0.001406*pow(w,5.)+0.0009216*pow(w,6.)+0.00064730*pow(w,7.)+0.0004779*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C8t2mt(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.35) return -0.8954-7.043*z-98.34*z*z-46.21*z*z*log(z)-127.1*pow(z,3.)
	-181.6*pow(z,3.)*log(z)+535.8*pow(z,4.)-76.76*pow(z,4.)*log(z);
	
	else return -0.6141-0.8975*w-0.03492*w*w+0.06791*pow(w,3.)+0.07966*pow(w,4.)
	+0.07226*pow(w,5.)+0.06132*pow(w,6.)+0.05096*pow(w,7.)+0.04216*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C8c2MW(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.35) return -1.870+0.1010*z-0.1218*z*log(z)+0.1045*z*z-0.03748*z*z*log(z)
	+0.01151*pow(z*log(z),2.)-0.01023*pow(z,3.)+0.004342*pow(z,3.)*log(z)+0.0003031*z*pow(z*log(z),2.)
	-0.001537*pow(z,4.)+0.0007532*pow(z,4.)*log(z);
	
	else return -1.676-0.1179*w-0.02926*w*w-0.01297*pow(w,3.)-0.007296*pow(w,4.)
	-0.004672*pow(w,5.)-0.003248*pow(w,6.)-0.002389*pow(w,7.)-0.001831*pow(w,8.);
}
