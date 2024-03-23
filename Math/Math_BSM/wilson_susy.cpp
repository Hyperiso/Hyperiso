#include <cmath>
#include "Math.h"
double C7H2lulumt(double r)
{
	double ubar=1.-r;
	double u=1.-1./r;
	
	if(r<0.16) return -8.27754*r-480.917*r*r-1158.21*pow(r,3.)-1491.65*pow(r,4.)-823.03*pow(r,5.)+4.31733*r*log(r)-396.082*r*r*log(r)-1292.39*pow(r,3.)*log(r)-2540.33*pow(r,4.)*log(r)-3362.49*pow(r,5.)*log(r)+0.922497*r*log(r)*log(r)-112.353*pow(r*log(r),2.)-348.239*pow(r,3.)*log(r)*log(r)-541.372*pow(r,4.)*log(r)*log(r)-412.426*pow(r,5.)*log(r)*log(r)-20.7325*r*r*pow(log(r),3.)-34.5021*pow(r*log(r),3.)-23.2551*pow(r,4.)*pow(log(r),3.)+42.3018*pow(r,5.)*pow(log(r),3.);

	else if(r<=1.) return 1.2835-0.715834*ubar-0.303943*pow(ubar,2.)-0.154911*pow(ubar,3.)-0.0862524*pow(ubar,4.)-0.0502032*pow(ubar,5.)-0.0296977*pow(ubar,6.)-0.0173999*pow(ubar,7.)-0.00975217*pow(ubar,8.)-0.00487721*pow(ubar,9.)-0.0017209*pow(ubar,10.)+0.000337773*pow(ubar,11.)+0.00167884*pow(ubar,12.)+0.00254213*pow(ubar,13.)+0.00308287*pow(ubar,14.)+0.0034037*pow(ubar,15.)+0.00357366*pow(ubar,16.);
	
	else if (r<10.5) return 1.2835+0.715834*u+0.411891*pow(u,2.)+0.262859*pow(u,3.)+0.182486*pow(u,4.)+0.134723*pow(u,5.)+0.104026*pow(u,6.)+0.0830584*pow(u,7.)+0.0680429*pow(u,8.)+0.0568813*pow(u,9.)+0.0483331*pow(u,10.)+0.0416253*pow(u,11.)+0.0362548*pow(u,12.)+0.0318817*pow(u,13.)+0.0282691*pow(u,14.)+0.0252475*pow(u,15.)+0.0226928*pow(u,16.);
	
	else return 3.96988-274.233/pow(r,5.)+24.4124/pow(r,4.)+79.1456/pow(r,3.)+47.0925/pow(r,2.)+15.3541/r-(72.1251*log(r))/pow(r,5.)-(168.261*log(r))/pow(r,4.)-(103.794*log(r))/pow(r,3.)-(38.1211*log(r))/pow(r,2.)-(8.75279*log(r))/r;
}

/*----------------------------------------------------------------------*/

double C7H2ldlumt(double r)
{
	double ubar=1.-r;
	double u=1.-1./r;
	
	if(r<0.16) return -572.224*r-524.114*pow(r,2.)+166.709*pow(r,3.)+1479.91*pow(r,4.)+2828.07*pow(r,5.)-453.485*r*log(r)-870.291*pow(r,2.)*log(r)-826.157*pow(r,3.)*log(r)+169.874*pow(r,4.)*log(r)+1985.76*pow(r,5.)*log(r)-123.519*r*pow(log(r),2.)-195.667*pow(r,2.)*pow(log(r),2.)-46.6111*pow(r,3.)*pow(log(r),2.)+323.191*pow(r,4.)*pow(log(r),2.)+469.396*pow(r,5.)*pow(log(r),2.)-20.9383*r*pow(log(r),3.)-8.88889*pow(r,2.)*pow(log(r),3.)+19.7284*pow(r,3.)*pow(log(r),3.)+36.0768*pow(r,4.)*pow(log(r),3.)-66.631*pow(r,5.)*pow(log(r),3.);

	else if(r<=1.) return 12.8225+1.66323*ubar+0.77799*pow(ubar,2.)+0.375487*pow(ubar,3.)+0.158083*pow(ubar,4.)+0.0302107*pow(ubar,5.)-0.0486763*pow(ubar,6.)-0.0986447*pow(ubar,7.)-0.130636*pow(ubar,8.)-0.15103*pow(ubar,9.)-0.163732*pow(ubar,10.)-0.171227*pow(ubar,11.)-0.175145*pow(ubar,12.)-0.17658*pow(ubar,13.)-0.176283*pow(ubar,14.)-0.174777*pow(ubar,15.)-0.172432*pow(ubar,16.);
	
	else if (r<10.5) return 12.8225-1.66323*u-0.885238*pow(u,2.)-0.482735*pow(u,3.)-0.297636*pow(u,4.)-0.202068*pow(u,5.)-0.147045*pow(u,6.)-0.112501*pow(u,7.)-0.0893095*pow(u,8.)-0.0729095*pow(u,9.)-0.0608295*pow(u,10.)-0.0516382*pow(u,11.)-0.0444588*pow(u,12.)-0.038729*pow(u,13.)-0.0340731*pow(u,14.)-0.030232*pow(u,15.)-0.0270216*pow(u,16.);
	
	else return 8.08756+194.326/pow(r,5.)-24.9665/pow(r,4.)-78.8992/pow(r,3.)-49.3185/pow(r,2.)-12.9103/r+(101.082*log(r))/pow(r,5.)+(168.442*log(r))/pow(r,4.)+(106.216*log(r))/pow(r,3.)+(38.4258*log(r))/pow(r,2.)+(9.75728*log(r))/r;	
}


double C7H2(double yt, double lu, double ld, double L)
{
	if(fabs(1.-yt)<1.e-5) return C7H2(0.9999,lu,ld,L);

	return lu*lu*(C7H2lulumt(yt)+L*(-yt*(67930.*pow(yt,4.)-470095.*pow(yt,3.)+1358478.*yt*yt-700243.*yt+54970.)/2187./pow(yt-1.,5.)+yt*(10422.*pow(yt,4.)-84390.*pow(yt,3.)+322801.*yt*yt-146588.*yt+1435.)/729./pow(yt-1.,6.)*log(yt)+2.*yt*yt*(260.*pow(yt,3.)-1515.*yt*yt+3757.*yt-1446.)/27./pow(yt-1.,5.)*Li2(1.-1./yt))
	+L*L*(yt*(-518.*pow(yt,4.)+3665.*pow(yt,3.)-17397.*yt*yt+3767.*yt+1843.)/162./pow(yt-1.,5.)+yt*yt*(-63.*pow(yt,3.)+532.*yt*yt+2089.*yt-1118.)/27./pow(yt-1.,6.)*log(yt))
	)
	+ld*lu*(C7H2ldlumt(yt)+L*(yt*(3790.*pow(yt,3.)-22511.*yt*yt+53614.*yt-21069.)/81./pow(yt-1.,4.)+2.*yt*(-1266.*pow(yt,3.)+7642.*yt*yt-21467.*yt+8179.)/81./pow(yt-1.,5.)*log(yt)-8.*yt*(139.*pow(yt,3.)-612.*yt*yt+1103.*yt-342.)/27./pow(yt-1.,4.)*Li2(1.-1./yt))
	+L*L*(yt*(284.*pow(yt,3.)-1435.*yt*yt+4304.*yt-1425.)/27./pow(yt-1.,4.)+2.*yt*(63.*pow(yt,3.)-397.*yt*yt-970.*yt+440.)/27./pow(yt-1.,5.)*log(yt))
	);
}

double C8H2lulumt(double r)
{
	double ubar=1.-r;
	double u=1.-1./r;
	
	if(r<0.16) return 0.743716*r-805.468*pow(r,2.)-3357.02*pow(r,3.)-9016.34*pow(r,4.)-19606.1*pow(r,5.)+3.23783*r*log(r)-602.665*pow(r,2.)*log(r)-3077.38*pow(r,3.)*log(r)-10102.2*pow(r,4.)*log(r)-26090.2*pow(r,5.)*log(r)+0.690844*r*pow(log(r),2.)-169.105*pow(r,2.)*pow(log(r),2.)-779.588*pow(r,3.)*pow(log(r),2.)-2243.8*pow(r,4.)*pow(log(r),2.)-5251.1*pow(r,5.)*pow(log(r),2.)-22.9807*pow(r,2.)*pow(log(r),3.)-66.3202*pow(r,3.)*pow(log(r),3.)-143.424*pow(r,4.)*pow(log(r),3.)-226.68*pow(r,5.)*pow(log(r),3.);
	else if(r<=1.) return 1.18809-0.407808*ubar-0.20763*pow(ubar,2.)-0.126464*pow(ubar,3.)-0.0856991*pow(ubar,4.)-0.0620372*pow(ubar,5.)-0.0468892*pow(ubar,6.)-0.0365182*pow(ubar,7.)-0.0290741*pow(ubar,8.)-0.0235439*pow(ubar,9.)-0.0193272*pow(ubar,10.)-0.0160457*pow(ubar,11.)-0.0134495*pow(ubar,12.)-0.0113673*pow(ubar,13.)-0.00967771*pow(ubar,14.)-0.00829289*pow(ubar,15.)-0.00714781*pow(ubar,16.);
	else if (r<10.5) return 1.18809+0.407808*u+0.200178*pow(u,2.)+0.119012*pow(u,3.)+0.0786112*pow(u,4.)+0.0553134*pow(u,5.)+0.0406051*pow(u,6.)+0.0307493*pow(u,7.)+0.023859*pow(u,8.)+0.0188842*pow(u,9.)+0.0151985*pow(u,10.)+0.0124088*pow(u,11.)+0.0102586*pow(u,12.)+0.00857494*pow(u,13.)+0.00723815*pow(u,14.)+0.00616365*pow(u,15.)+0.00529038*pow(u,16.);	
	else return 0.27835+826.151/pow(r,5.)+96.3468/pow(r,4.)-66.3939/pow(r,3.)-39.7588/pow(r,2.)-5.21409/r-(300.663*log(r))/pow(r,5.)+(91.89*log(r))/pow(r,4.)+(78.5823*log(r))/pow(r,3.)+(20.0187*log(r))/pow(r,2.);
}

/*----------------------------------------------------------------------*/

double C8H2ldlumt(double r)
{
	double ubar=1.-r;
	double u=1.-1./r;
	
	if(r<0.16) return -929.846*r-2942.88*pow(r,2.)-6480.66*pow(r,3.)-11683.8*pow(r,4.)-15961.9*pow(r,5.)-658.39*r*log(r)-2769.99*pow(r,2.)*log(r)-7906.35*pow(r,3.)*log(r)-17770.*pow(r,4.)*log(r)-29962.3*pow(r,5.)*log(r)-174.653*r*pow(log(r),2.)-612.625*pow(r,2.)*pow(log(r),2.)-1438.79*pow(r,3.)*pow(log(r),2.)-2776.64*pow(r,4.)*pow(log(r),2.)-2626.73*pow(r,5.)*pow(log(r),2.)-19.7963*r*pow(log(r),3.)-31.8333*pow(r,2.)*pow(log(r),3.)-40.6759*pow(r,3.)*pow(log(r),3.)+54.6584*pow(r,4.)*pow(log(r),3.)+1002.93*pow(r,5.)*pow(log(r),3.);
	else if(r<=1.) return -0.610999+1.09548*ubar+0.649151*pow(ubar,2.)+0.459582*pow(ubar,3.)+0.356928*pow(ubar,4.)+0.291027*pow(ubar,5.)+0.243786*pow(ubar,6.)+0.207532*pow(ubar,7.)+0.178511*pow(ubar,8.)+0.154649*pow(ubar,9.)+0.134673*pow(ubar,10.)+0.117736*pow(ubar,11.)+0.103236*pow(ubar,12.)+0.0907291*pow(ubar,13.)+0.0798733*pow(ubar,14.)+0.0704007*pow(ubar,15.)+0.0620968*pow(ubar,16.);
	else if (r<10.5) return -0.610999-1.09548*u-0.446333*pow(u,2.)-0.256764*pow(u,3.)-0.169848*pow(u,4.)-0.119684*pow(u,5.)-0.0876105*pow(u,6.)-0.0659531*pow(u,7.)-0.0507877*pow(u,8.)-0.0398674*pow(u,9.)-0.0318194*pow(u,10.)-0.0257687*pow(u,11.)-0.0211391*pow(u,12.)-0.0175414*pow(u,13.)-0.0147063*pow(u,14.)-0.0124441*pow(u,15.)-0.0106184*pow(u,16.);
	else return -3.17375-1002.87/pow(r,5.)-205.702/pow(r,4.)+62.2604/pow(r,3.)+63.7396/pow(r,2.)+10.8891/r+(476.897*log(r))/pow(r,5.)-(71.6153*log(r))/pow(r,4.)-(110.665*log(r))/pow(r,3.)-(35.4207*log(r))/pow(r,2.);
}

/*----------------------------------------------------------------------*/

double C8H2(double yt, double lu, double ld, double L)
{
	if(fabs(1.-yt)<1.e-5) return C8H2(0.9999,lu,ld,L);

	return lu*lu*(C8H2lulumt(yt)+L*(yt*(51948.*pow(yt,4.)-233781.*pow(yt,3.)+48634.*yt*yt-698693.*yt+2452.)/1944./pow(yt-1.,6.)*log(yt)-yt*(522347.*pow(yt,4.)-2423255.*pow(yt,3.)+2706021.*yt*yt-5930609.*yt+148856.)/11664./pow(yt-1.,5.)+yt*yt*(481.*pow(yt,3.)-1950.*yt*yt+1523.*yt-2550.)/18./pow(yt-1.,5.)*Li2(1.-1./yt))
	+L*L*(yt*(-259.*pow(yt,4.)+1117.*pow(yt,3.)+2925.*yt*yt+28411.*yt+2366.)/216./pow(yt-1.,5.)-yt*yt*(139.*yt*yt+2938.*yt+2683.)/36./pow(yt-1.,6.)*log(yt))
	)
	+ld*lu*(C8H2ldlumt(yt)+L*(yt*(1463.*pow(yt,3.)-5794.*yt*yt+5543.*yt-15036.)/27./pow(yt-1.,4.)+yt*(-1887.*pow(yt,3.)+7115.*yt*yt+2519.*yt+19901.)/54./pow(yt-1.,5.)*log(yt)+yt*(-629.*pow(yt,3.)+2178.*yt*yt-1729.*yt+2196.)/18./pow(yt-1.,4.)*Li2(1.-1./yt))
	+L*L*(yt*(259.*pow(yt,3.)-947.*yt*yt-251.*yt-5973.)/36./pow(yt-1.,4.)+yt*(139.*yt*yt+2134.*yt+1183.)/18./pow(yt-1.,5.)*log(yt))
	);
}