#include <cmath>
#include "Math.h"


double F1SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F1SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F1SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F1SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F1SP(xH*0.9998,xH);

	return 1./4./(xH-xt)*(xt*log(xt)/(xt-1.)-xH*log(xH)/(xH-1.));
}


/*--------------------------------------------------------------------*/

double F2SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F2SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F2SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F2SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F2SP(xH*0.9998,xH);

	return 1./8./(xH-xt)*(xH/(xH-1.)+xt*xt*log(xt)/(xt-1.)/(xH-xt)-xH*(xH*xt+xH-2.*xt)/(xH-1.)/(xH-1.)/(xH-xt)*log(xH));
}

/*--------------------------------------------------------------------*/

double F3SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F3SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F3SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F3SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F3SP(xH*0.9998,xH);

	return 1./8./(xH-xt)*((xH-xt)/(xH-1.)/(xt-1.)+xt*(xt-2.)*log(xt)/(xt-1.)/(xt-1.)-xH*(xH-2.)/(xH-1.)/(xH-1.)*log(xH));
}

/*--------------------------------------------------------------------*/

double F4SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F4SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F4SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F4SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F4SP(xH*0.9998,xH);

	return xt/(xH-xt)*(1.-xH/(xH-xt)*log(xH/xt));
}

/*--------------------------------------------------------------------*/

double F5SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F5SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F5SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F5SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F5SP(xH*0.9998,xH);

	return xt/2./(xH-xt)/(xH-xt)*((xH+xt)/2.-xH*xt/(xH-xt)*log(xH/xt));
}

/*--------------------------------------------------------------------*/

double F6SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F6SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F6SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F6SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F6SP(xH*0.9998,xH);

	return 1./2./(xH-xt)*(-xH+xt+xH*log(xH)-xt*log(xt));
}

/*--------------------------------------------------------------------*/

double F7SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F7SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F7SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F7SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F7SP(xH*0.9998,xH);

	return 1./2./(xH-xt)*(xt-xH*xt/(xH-xt)*log(xH/xt));
}

/*--------------------------------------------------------------------*/

double F8SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F8SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F8SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F8SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F8SP(xH*0.9998,xH);

	return 1./2./(xH-xt)*(xH-xH*xH*log(xH)/(xH-xt)+xt*(2.*xH-xt)*log(xt)/(xH-xt));
}

/*--------------------------------------------------------------------*/

double F9SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F9SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F9SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F9SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F9SP(xH*0.9998,xH);

	return 1./4./(xH-xt)/(xH-xt)*(xt*(3.*xH-xt)/2.-xH*xH*xt/(xH-xt)*log(xH/xt));
}

/*--------------------------------------------------------------------*/

double F10SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F10SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F10SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F10SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F10SP(xH*0.9998,xH);

	return 1./4./(xH-xt)/(xH-xt)*(xt*(xH-3.*xt)/2.-xH*xt*(xH-2.*xt)/(xH-xt)*log(xH/xt));
}

/*--------------------------------------------------------------------*/

double F11SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F11SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F11SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F11SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F11SP(xH*0.9998,xH);
	
	double xt2=xt*xt;
	double xt3=xt*xt2;

	return 1./2./(xH-xt)*(xt*(xt2-3.*xH*xt+9.*xH-5.*xt-2.)/4./(xt-1.)/(xt-1.)+xH*(xH*xt-3.*xH+2.*xt)/2./(xH-1.)/(xH-xt)*log(xH)
	+(xH*xH*(-2.*xt3+6.*xt2-9.*xt+2.)+3.*xH*xt2*(xt2-2.*xt+3.)-xt2*(2.*xt3-2.*xt2+3.*xt+1.))*log(xt)/2./pow(xt-1.,3.)/(xH-xt));
}

/*--------------------------------------------------------------------*/

double F12SP(double xt, double xH)
{
	if((fabs(1.-xt)<1.e-5)&&(fabs(1.-xH)<1.e-5)) return F12SP(0.9998,1.0002);

	if(fabs(1.-xt)<1.e-5) return F12SP(0.9999,xH);
	if(fabs(1.-xH)<1.e-5) return F12SP(xt,0.9999);

	if(fabs(1.-xt/xH)<1.e-5) return F12SP(xH*0.9998,xH);
	
	double xt2=xt*xt;
	double xt3=xt*xt2;

	return 1./2./(xH-xt)*((xt2+xt-8.)*(xH-xt)/4./(xt-1.)/(xt-1.)-xH*(xH+2.)*log(xH)/2./(xH-1.)
	+(xH*(xt3-3.*xt2+3.*xt+2.)+3.*(xt-2.)*xt2)*log(xt)/2./pow(xt-1.,3.));
}

/*--------------------------------------------------------------------*/

double CSc_2HDM(double xH, double xt, double lu, double ld, double ll)
{
		if((fabs(1.-xH)<1.e-5)&&(fabs(1.-xt)<1.e-5)) return CSc_2HDM(0.9998,1.0002,lu,ld,ll);

		if(fabs(1.-xH)<1.e-5) return CSc_2HDM(0.9999,xt,lu,ld,ll);
		if(fabs(1.-xt)<1.e-5) return CSc_2HDM(xH,0.9999,lu,ld,ll);

		if(fabs(1.-xH/xt)<1.e-5) return CSc_2HDM(xt*0.9998,xt,lu,ld,ll);
		
		return xt/8./(xH-xt)*(2.*ld*ll*(1./(xH-1.)*log(xH)-1./(xt-1.)*log(xt))
		+lu*ll*(1./(xH-1.)+xH/(xH-xt)/(xt-1.)*log(xt)-xH*(2.*xH-xt-1.)/(xH-xt)/(xH-1.)/(xH-1.)*log(xH))
		-ll*lu*((xt-xH)/(xH-1.)/(xt-1.)+xt/(xt-1.)/(xt-1.)*log(xt)-xH/(xH-1.)/(xH-1.)*log(xH)));
}

/*--------------------------------------------------------------------*/

double CPc_2HDM(double xH, double xt, double lu, double ld, double ll, double sw2)
{
		if((fabs(1.-xH)<1.e-5)&&(fabs(1.-xt)<1.e-5)) return CPc_2HDM(0.9998,1.0002,lu,ld,ll,sw2);

		if(fabs(1.-xH)<1.e-5) return CPc_2HDM(0.9999,xt,lu,ld,ll,sw2);
		if(fabs(1.-xt)<1.e-5) return CPc_2HDM(xH,0.9999,lu,ld,ll,sw2);

		if(fabs(1.-xH/xt)<1.e-5) return CPc_2HDM(xt*0.9998,xt,lu,ld,ll,sw2);
		
		return -xt/8./(xH-xt)*(2.*ld*ll*(1./(xH-1.)*log(xH)-1./(xt-1.)*log(xt))
		+lu*ll*(1./(xH-1.)+xH/(xH-xt)/(xt-1.)*log(xt)-xH*(2.*xH-xt-1.)/(xH-xt)/(xH-1.)/(xH-1.)*log(xH))
		+ll*lu*((xt-xH)/(xH-1.)/(xt-1.)+xt/(xt-1.)/(xt-1.)*log(xt)-xH/(xH-1.)/(xH-1.)*log(xH)))
		+xt/4./(xH-xt)/(xH-xt)*(-ld*lu*(-(xt+xH)/2.+xt*xH/(xH-xt)*log(xH/xt))
		+lu*lu/6./(xH-xt)*((xH*xH-8.*xH*xt-17.*xt*xt)/6.+xt*xt*(3.*xH+xt)/(xH-xt)*log(xH/xt)))
		+sw2*xt/6./(xH-xt)/(xH-xt)*(-ld*lu*((5.*xt-3.*xH)/2.+xH*(2.*xH-3.*xt)/(xH-xt)*log(xH/xt))
		+lu*lu/6./(xH-xt)*((4.*xH*xH*xH-12.*xH*xH*xt+xH*xt*xt+3.*xt*xt*xt)/(xH-xt)*log(xH/xt)-(17.*xH*xH-64.*xH*xt+71.*xt*xt)/6.))
		+(1.-sw2)*lu*lu*xt*xt/4./(xH-xt)/(xH-xt)*(xH*log(xH/xt)+xt-xH);
}
